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
    // !!!cmk later or file_name or a Path,
    // !!!cmk later confirm that this is required by derive_builder
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
            // !!!cmk 0 understand this as_ref
            // !!!cmk 0 rename filename to path
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
    ) -> Result<impl Iterator<Item = String>, BedErrorPlus> {
        // !!!cmk later allow fam file to be specified.
        let path_buf = Path::new(&self.path).with_extension(suffix);
        // !!!cmk later here and elsewhere if there are only two arms, use 'if let' maybe?
        let file = match File::open(&path_buf) {
            Err(_) => {
                return Err(BedError::CannotOpenFamOrBim(path_buf.display().to_string()).into());
            }
            Ok(file) => file,
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
        let iter = reader.lines().map(move |line| {
            // !!!cmk later unwrap
            let line = line.unwrap();
            // !!!cmk later unwrap
            let field = line.split_whitespace().nth(field_index).unwrap();
            field.to_string()
        });

        Ok(iter)
    }

    // !!!cmk later should not have any unwraps in this whole file
    pub fn get_iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if self.iid.is_some() {
            // !!!cmk later unwrap
            Ok(self.iid.as_ref().unwrap())
        } else {
            let iid = self.read_fam_or_bim("fam", 1)?.collect();
            self.iid = Some(iid);
            // !!!cmk later unwrap
            Ok(self.iid.as_ref().unwrap())
        }
    }

    pub fn get_sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if self.sid.is_some() {
            // !!!cmk later unwrap
            Ok(self.sid.as_ref().unwrap())
        } else {
            let sid = self.read_fam_or_bim("bim", 1)?.collect();
            self.sid = Some(sid);
            // !!!cmk later unwrap
            Ok(self.sid.as_ref().unwrap())
        }
    }

    pub fn get_chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        if self.chromosome.is_some() {
            // !!!cmk later unwrap
            Ok(self.chromosome.as_ref().unwrap())
        } else {
            let chromosome = self.read_fam_or_bim("bim", 0)?.collect();
            self.chromosome = Some(chromosome);
            // !!!cmk later unwrap
            Ok(self.chromosome.as_ref().unwrap())
        }
    }

    pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing>(
        &self,
    ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_options = ReadOptions::<TOut>::builder().build()?;
        self.read_with_options(read_options)
    }

    pub(crate) fn read_with_options<
        TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing,
    >(
        &self,
        read_options: ReadOptions<TOut>,
    ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        // !!!cmk later this should be lazy in Bed object, not here
        let (iid_count, sid_count) = counts(&self.path)?;

        // !!!cmk later do something with read_options.num_threads

        let iid_index = to_vec(read_options.iid_index, iid_count);
        let sid_index = to_vec(read_options.sid_index, sid_count);

        let shape = ShapeBuilder::set_f(
            (iid_index.len(), sid_index.len()),
            read_options.output_is_orderf,
        );
        let mut val = nd::Array2::<TOut>::default(shape);

        read_no_alloc(
            &self.path,
            iid_count,
            sid_count,
            self.count_a1,
            &iid_index,
            &sid_index,
            read_options.missing_value,
            &mut val.view_mut(),
        )?;

        Ok(val)
    }
}

// !!!cmk later implement this with conversion impl
pub fn to_vec(index: Index, count: usize) -> Vec<usize> {
    match index {
        Index::None => (0..count).collect(),
        Index::Vec(iid_index) => iid_index,
        // !!!cmk ask how to support slices
        // !!!cmk ask Index::RustSlice(iid_index) => iid_index.to_vec(),
        Index::Bool(bool_index) => {
            // !!!cmk later check that bool_index.len() == iid_count
            // !!!cmk later use enumerate() instead of zip
            (0..count)
                .zip(bool_index)
                .filter(|(_, b)| *b)
                .map(|(i, _)| i)
                .collect()
        }
        // !!!cmk later can we implement this without two allocations?
        Index::NDSlice(slice_index) => {
            let full_array: nd::Array1<usize> = (0..count).collect();
            let array = full_array.slice(slice_index);
            array.to_vec()
        }
    }
}

pub(crate) type SliceInfo1 =
    nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>;

// !!!cmk later add docs to type typedbuilder stuff: https://docs.rs/typed-builder/latest/typed_builder/derive.TypedBuilder.html#customisation-with-attributes
#[derive(Debug, Clone)]
pub enum Index {
    None,
    Vec(Vec<usize>),
    Bool(nd::Array1<bool>),
    NDSlice(SliceInfo1),
    // !!! cmk0 what about supporting ranges?
}

// !!!cmk later see if what ref conversions. See https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust
impl From<SliceInfo1> for Index {
    fn from(slice_info: SliceInfo1) -> Index {
        Index::NDSlice(slice_info)
    }
}

impl From<Vec<usize>> for Index {
    fn from(full: Vec<usize>) -> Index {
        Index::Vec(full)
    }
}

impl From<nd::Array1<bool>> for Index {
    fn from(bool_array: nd::Array1<bool>) -> Index {
        Index::Bool(bool_array)
    }
}

impl From<()> for Index {
    fn from(_: ()) -> Index {
        Index::None
    }
}

// !!!cmk  later "Arg" is unlikely to be a good name
// See https://nullderef.com/blog/rust-parameters/
// !!!cmk later note that ndarray can do this: a.slice(s![1..4;2, ..;-1])
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
    pub fn read(&self, bed: &Bed) -> Result<nd::Array2<TOut>, BedErrorPlus> {
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
