use core::fmt::Debug;
use nd::ShapeBuilder;
use ndarray as nd;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use derive_builder::Builder;
// !!! might want to use this instead use typed_builder::TypedBuilder;

use crate::{
    counts, read_no_alloc, BedError, BedErrorPlus, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
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
    pub filename: String, // !!!cmk always clone?

    #[builder(default = "true")]
    count_a1: bool,

    #[builder(default = "None", setter(strip_option))]
    iid_count: Option<usize>,

    #[builder(default = "None", setter(strip_option))]
    sid_count: Option<usize>,

    #[builder(default = "None", setter(strip_option))]
    iid: Option<Vec<String>>,

    #[builder(default = "None", setter(strip_option))]
    sid: Option<Vec<String>>,

    #[builder(default = "None", setter(strip_option))]
    chromosome: Option<Vec<String>>,
}

impl BedBuilder {
    pub fn build(&self) -> Result<Bed, BedErrorPlus> {
        let bed = Bed {
            filename: match self.filename {
                Some(ref value) => Clone::clone(value),
                None => {
                    panic!("filename is required");
                    // !!!cmk later put this code back
                    // return Result::Err(Into::into(
                    //     ::derive_builder::UninitializedFieldError::from("filename"),
                    // ))
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
        let filename = bed.filename.to_string();
        let mut buf_reader = BufReader::new(File::open(&filename).unwrap());
        let mut bytes_vector: Vec<u8> = vec![0; CB_HEADER_USIZE];
        buf_reader.read_exact(&mut bytes_vector).unwrap();
        if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
            return Err(BedError::IllFormed(filename.to_string()).into());
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
    pub fn get_iid_count(&mut self) -> usize {
        if let Some(iid_count) = self.iid_count {
            iid_count
        } else {
            let (iid_count, sid_count) = counts(&self.filename).unwrap();
            self.iid_count = Some(iid_count);
            self.sid_count = Some(sid_count);
            iid_count
        }
    }
    pub fn get_sid_count(&mut self) -> usize {
        if let Some(sid_count) = self.sid_count {
            sid_count
        } else {
            let (iid_count, sid_count) = counts(&self.filename).unwrap();
            self.iid_count = Some(iid_count);
            self.sid_count = Some(sid_count);
            sid_count
        }
    }

    /// !!!cmk later don't re-read for every column
    fn read_fam_or_bim(
        &self,
        suffix: &str,
        field_index: usize,
    ) -> Result<impl Iterator<Item = String>, BedErrorPlus> {
        // !!!cmk later allow fam file to be specified.
        let path_buf = Path::new(&self.filename).with_extension(suffix);
        // !!!cmk later here and elsewhere if there are only two arms, use 'if let' maybe?
        let file = match File::open(&path_buf) {
            Err(_) => {
                let string_path = path_buf.to_string_lossy().to_string();
                return Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                    string_path,
                )));
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
            let line = line.unwrap();
            let field = line.split_whitespace().nth(field_index).unwrap();
            field.to_string()
        });

        Ok(iter)
    }

    // !!!cmk later should not have any unwraps in this whole file
    pub fn get_iid(&mut self) -> &Vec<String> {
        if self.iid.is_some() {
            self.iid.as_ref().unwrap()
        } else {
            let iid = self.read_fam_or_bim("fam", 1).unwrap().collect();
            self.iid = Some(iid);
            self.iid.as_ref().unwrap()
        }
    }

    pub fn get_sid(&mut self) -> &Vec<String> {
        if self.sid.is_some() {
            self.sid.as_ref().unwrap()
        } else {
            let sid = self.read_fam_or_bim("bim", 1).unwrap().collect();
            self.sid = Some(sid);
            self.sid.as_ref().unwrap()
        }
    }

    pub fn get_chromosome(&mut self) -> &Vec<String> {
        if self.chromosome.is_some() {
            self.chromosome.as_ref().unwrap()
        } else {
            let chromosome = self.read_fam_or_bim("bim", 0).unwrap().collect();
            self.chromosome = Some(chromosome);
            self.chromosome.as_ref().unwrap()
        }
    }

    pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send>(
        &self,
        read_arg: ReadArg<TOut>,
    ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        // !!!cmk later this should be lazy in Bed object, not here
        let (iid_count, sid_count) = counts(&self.filename)?;

        // !!!cmk later do something with read_arg.num_threads

        let iid_index = to_vec(read_arg.iid_index, iid_count);
        let sid_index = to_vec(read_arg.sid_index, sid_count);

        let shape = ShapeBuilder::set_f(
            (iid_index.len(), sid_index.len()),
            read_arg.output_is_orderf,
        );
        let mut val = nd::Array2::<TOut>::default(shape);

        read_no_alloc(
            &self.filename,
            iid_count,
            sid_count,
            self.count_a1,
            &iid_index,
            &sid_index,
            read_arg.missing_value,
            &mut val.view_mut(),
        )?;

        Ok(val)
    }
}

// !!!cmk 0 avoid saying Index:: so much
pub fn to_vec(index: Index, count: usize) -> Vec<usize> {
    match index {
        // !!!cmk 0 faster to add slice support to Index???
        Index::None => (0..count).collect(),
        Index::Full(iid_index) => iid_index,
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
        Index::Slice(slice_index) => {
            let full_array: nd::Array1<usize> = (0..count).collect();
            let array = full_array.slice(slice_index);
            array.to_vec()
        }
    }
}

// !!!cmk later add docs to type typedbuilder stuff: https://docs.rs/typed-builder/latest/typed_builder/derive.TypedBuilder.html#customisation-with-attributes
#[derive(Debug, Clone)]
pub enum Index {
    None,
    Full(Vec<usize>),
    Bool(Vec<bool>),
    Slice(nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>),
    // !!! cmk0 what about supporting ranges?
}

// !!!cmk  later "Arg" is unlikely to be a good name
// See https://nullderef.com/blog/rust-parameters/
// !!!cmk 0 note that ndarray can do this: a.slice(s![1..4;2, ..;-1])
#[derive(Builder)]
pub struct ReadArg<TOut: Copy + Default + From<i8> + Debug + Sync + Send> {
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
