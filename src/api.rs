// References: https://www.youtube.com/watch?v=0zOg8_B71gE&t=22s
// https://deterministic.space/elegant-apis-in-rust.html
// https://rust-lang.github.io/api-guidelines/
// https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust

use core::fmt::Debug;
use nd::ShapeBuilder;
use ndarray as nd;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    ops::Range,
    path::{Path, PathBuf},
};

use derive_builder::Builder;
// !!! might want to use this instead use typed_builder::TypedBuilder;

use crate::{
    counts, read_no_alloc, BedError, BedErrorPlus, Missing, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
    CB_HEADER_USIZE,
};

// https://crates.io/crates/typed-builder
// (or https://docs.rs/derive_builder/latest/derive_builder/)
// Somehow ndarray can do this: 	Array::zeros((3, 4, 5).f())
//       see https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html
#[derive(Builder)]
#[builder(build_fn(skip))]
pub struct Bed {
    // https://stackoverflow.com/questions/32730714/what-is-the-right-way-to-store-an-immutable-path-in-a-struct
    pub path: PathBuf, // !!!cmk later always clone?

    #[builder(default = "true")]
    count_a1: bool,

    #[builder(default = "None", setter(strip_option))]
    iid_count: Option<usize>,

    #[builder(default = "None", setter(strip_option))]
    sid_count: Option<usize>,

    #[builder(default = "None", setter(strip_option))]
    iid: Option<nd::Array1<String>>,

    #[builder(default = "None", setter(strip_option))]
    sid: Option<nd::Array1<String>>,

    #[builder(default = "None", setter(strip_option))]
    chromosome: Option<nd::Array1<String>>,
}

impl Bed {
    pub fn builder<P: AsRef<Path>>(path: P) -> BedBuilder {
        BedBuilder::new(path)
    }

    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedErrorPlus> {
        Bed::builder(path).build()
    }
}

impl BedBuilder {
    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        Self {
            path: Some(PathBuf::from(path.as_ref())),
            count_a1: Some(true),
            iid_count: None,
            sid_count: None,
            iid: None,
            sid: None,
            chromosome: None,
        }
    }

    pub fn build(&self) -> Result<Bed, BedErrorPlus> {
        let bed = Bed {
            path: match self.path {
                Some(ref value) => Clone::clone(value),
                None => {
                    return Err(BedErrorPlus::UninitializedFieldError(
                        ::derive_builder::UninitializedFieldError::from("path"),
                    ))
                }
            },
            count_a1: match self.count_a1 {
                Some(ref value) => Clone::clone(value),
                None => true,
            },
            iid_count: match self.iid_count {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            sid_count: match self.sid_count {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            iid: match self.iid {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            sid: match self.sid {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            chromosome: match self.chromosome {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
        };

        // !!!cmk later similar code elsewhere
        let mut buf_reader = BufReader::new(File::open(&bed.path)?);
        let mut bytes_vector: Vec<u8> = vec![0; CB_HEADER_USIZE];
        buf_reader.read_exact(&mut bytes_vector)?;
        if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
            return Err(BedError::IllFormed(PathBuf::from(&bed.path).display().to_string()).into());
        }

        Result::Ok(bed)
    }
}

impl Bed {
    // !!!cmk later
    // properties: Mapping[str, List[Any]] = {},
    // num_threads: Optional[int] = None,
    // skip_format_check: bool = False,
    // fam_filepath: Union[str, Path] = None,
    // bim_filepath: Union[str, Path] = None,

    // !!!cmk later is this how you do lazy accessors?
    // !!!cmk later should this be "try_get_..." or just "iid_count" or as is
    pub fn get_iid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let (iid_count, sid_count) = counts(&self.path)?;
            self.iid_count = Some(iid_count);
            self.sid_count = Some(sid_count);
            Ok(iid_count)
        }
    }
    pub fn get_sid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let (iid_count, sid_count) = counts(&self.path)?;
            self.iid_count = Some(iid_count);
            self.sid_count = Some(sid_count);
            Ok(sid_count)
        }
    }

    /// !!!cmk later don't re-read for every column
    fn read_fam_or_bim(
        &self,
        suffix: &str,
        field_index: usize,
    ) -> Result<nd::Array1<String>, BedErrorPlus> {
        // !!!cmk later allow fam file to be specified.
        let path_buf = Path::new(&self.path).with_extension(suffix);
        let file = if let Ok(file) = File::open(&path_buf) {
            file
        } else {
            return Err(BedError::CannotOpenFamOrBim(path_buf.display().to_string()).into());
        };

        // !!!cmk later use the correct delimiters (here is the Python spec:)
        // count = self._counts[suffix]
        // delimiter = _delimiters[suffix]
        // if delimiter in {r"\s+"}:
        //     delimiter = None
        //     delim_whitespace = True
        // else:
        //     delim_whitespace = False

        let reader = BufReader::new(file);
        let mut line_vec = Vec::<String>::new();
        for line in reader.lines() {
            let line = line?;
            let field = line.split_whitespace().nth(field_index);
            if let Some(field) = field {
                line_vec.push(field.to_string());
            } else {
                // !!!cmk later update to more specific error message
                return Err(BedError::CannotOpenFamOrBim(path_buf.display().to_string()).into());
            }
        }

        let array = match nd::Array1::from_shape_vec(line_vec.len(), line_vec) {
            Ok(array) => array,
            Err(_error) => {
                // !!!cmk later update to more specific error message
                return Err(BedError::CannotOpenFamOrBim(path_buf.display().to_string()).into());
                // return Err(error).into();
            }
        };

        Ok(array)
    }

    pub fn get_iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if let Some(ref iid) = self.iid {
            Ok(iid)
        } else {
            let iid = self.read_fam_or_bim("fam", 1)?;
            self.iid = Some(iid);
            // This unwrap is safe because we just created 'Some'
            Ok(self.iid.as_ref().unwrap())
        }
    }

    pub fn get_sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if let Some(ref sid) = self.sid {
            Ok(sid)
        } else {
            let sid = self.read_fam_or_bim("bim", 1)?;
            self.sid = Some(sid);
            // This unwrap is safe because we just created 'Some'
            Ok(self.sid.as_ref().unwrap())
        }
    }

    pub fn get_chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if let Some(ref chromosome) = self.chromosome {
            Ok(chromosome)
        } else {
            let chromosome = self.read_fam_or_bim("bim", 0)?;
            self.chromosome = Some(chromosome);
            // This unwrap is safe because we just created 'Some'
            Ok(self.chromosome.as_ref().unwrap())
        }
    }

    pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing>(
        &mut self,
    ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_options = ReadOptions::<TOut>::builder().build()?;
        self.read_with_options(read_options)
    }

    pub(crate) fn read_with_options<
        TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing,
    >(
        &mut self,
        read_options: ReadOptions<TOut>,
    ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let iid_count = self.get_iid_count()?;
        let sid_count = &self.get_sid_count()?;

        // !!!cmk later do something with read_options.num_threads

        let iid_index = to_vec(read_options.iid_index, iid_count);
        let sid_index = to_vec(read_options.sid_index, *sid_count);

        let shape = ShapeBuilder::set_f(
            (iid_index.len(), sid_index.len()),
            read_options.output_is_orderf,
        );
        let mut val = nd::Array2::<TOut>::default(shape);

        read_no_alloc(
            &self.path,
            iid_count,
            *sid_count,
            self.count_a1,
            &iid_index,
            &sid_index,
            read_options.missing_value,
            &mut val.view_mut(),
        )?;

        Ok(val)
    }
}

// !!!cmk later test every case
// !!!cmk 0 implement this with conversion impl
pub fn to_vec(index: Index, count: usize) -> Vec<usize> {
    match index {
        Index::None => (0..count).collect(),
        Index::Vec(vec) => vec,
        Index::NDArrayBool(nd_array_bool) => {
            // !!!cmk later check that bool_index.len() == iid_count
            // !!!cmk 0 use enumerate() instead of zip
            (0..count)
                .zip(nd_array_bool)
                .filter(|(_, b)| *b)
                .map(|(i, _)| i)
                .collect()
        }
        // !!!cmk later can we implement this without two allocations?
        Index::NDSliceInfo(nd_slice_info) => {
            let full_array: nd::Array1<usize> = (0..count).collect();
            let array = full_array.slice(nd_slice_info);
            array.to_vec()
        }
        Index::Range(range) => range.collect(),
        Index::NDArray(nd_array) => nd_array.to_vec(),
        Index::One(one) => vec![one],
        Index::VecBool(vec_bool) => {
            // !!!cmk later similar code elsewhere
            // !!!cmk later check that vec_bool.len() == iid_count
            // !!!cmk 0 use enumerate() instead of zip
            (0..count)
                .zip(vec_bool)
                .filter(|(_, b)| *b)
                .map(|(i, _)| i)
                .collect()
        }
    }
}

pub(crate) type SliceInfo1 =
    nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>;

// !!!cmk later add docs to type typedbuilder stuff: https://docs.rs/typed-builder/latest/typed_builder/derive.TypedBuilder.html#customisation-with-attributes
#[derive(Debug, Clone)]
pub enum Index {
    None,
    One(usize),
    Vec(Vec<usize>),
    NDArray(nd::Array1<usize>),
    VecBool(Vec<bool>),
    NDArrayBool(nd::Array1<bool>),
    NDSliceInfo(SliceInfo1),
    Range(Range<usize>),
}

// !!!cmk later see if what ref conversions. See https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust
impl From<SliceInfo1> for Index {
    fn from(slice_info: SliceInfo1) -> Index {
        Index::NDSliceInfo(slice_info)
    }
}

impl From<Range<usize>> for Index {
    fn from(range: Range<usize>) -> Index {
        Index::Range(range)
    }
}

impl From<&[usize]> for Index {
    fn from(slice: &[usize]) -> Index {
        Index::Vec(slice.to_vec())
    }
}

impl From<usize> for Index {
    fn from(one: usize) -> Index {
        Index::One(one)
    }
}

impl From<Vec<usize>> for Index {
    fn from(vec: Vec<usize>) -> Index {
        Index::Vec(vec)
    }
}
impl From<nd::Array1<usize>> for Index {
    fn from(nd_array: nd::Array1<usize>) -> Index {
        Index::NDArray(nd_array)
    }
}

impl From<nd::Array1<bool>> for Index {
    fn from(nd_array_bool: nd::Array1<bool>) -> Index {
        Index::NDArrayBool(nd_array_bool)
    }
}

impl From<Vec<bool>> for Index {
    fn from(vec_bool: Vec<bool>) -> Index {
        Index::VecBool(vec_bool)
    }
}

impl From<()> for Index {
    fn from(_: ()) -> Index {
        Index::None
    }
}

// See https://nullderef.com/blog/rust-parameters/
#[derive(Builder)]
pub struct ReadOptions<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing> {
    #[builder(default = "TOut::missing()")]
    missing_value: TOut,

    #[builder(default = "Index::None")]
    iid_index: Index,

    #[builder(default = "Index::None")]
    sid_index: Index,

    #[builder(default = "true")]
    output_is_orderf: bool, // !!!cmk later use enum or .f()

    #[builder(default, setter(strip_option))]
    pub num_threads: Option<usize>,
    // num_threads=None,
}

impl<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing + Clone> ReadOptions<TOut> {
    pub fn builder() -> ReadOptionsBuilder<TOut> {
        ReadOptionsBuilder::default()
    }
}

impl<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing + Clone>
    ReadOptionsBuilder<TOut>
{
    pub fn read(&self, bed: &mut Bed) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_option = self.build()?;
        bed.read_with_options(read_option)
    }
}

// !!!cmk ask could a macro likes be nice?
// #[macro_export]
// macro_rules! read {
//     ($bed:expr) => {
//         $bed.read()
//     };
//     ($bed:expr, $option: expr) => {
//         ReadOptions::builder().$option.read($bed)
//     };
// }
