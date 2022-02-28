// References: https://www.youtube.com/watch?v=0zOg8_B71gE&t=22s
// https://deterministic.space/elegant-apis-in-rust.html
// https://rust-lang.github.io/api-guidelines/
// https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust

// !!!cmk later write doc tests (see https://deterministic.space/elegant-apis-in-rust.html#what-makes-an-api-elegant)
// !!!cmk later To enforce that every public API item is documented, use #![deny(missing_docs)].
// !!!cmk later conventions for formatting Rust documentation https://deterministic.space/machine-readable-inline-markdown-code-cocumentation.html

// !!!cmk later document and add issue that File(s) are not held, incorrectly allowing for the file to be changed between reads.

use core::fmt::Debug;
use nd::ShapeBuilder;
use ndarray as nd;
use std::{
    env,
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Allele {
    A1,
    A2,
}

#[derive(Debug, Clone)]
pub enum OptionOrSkip<T> {
    None,
    Some(T),
    Skip,
}

impl<T> OptionOrSkip<T> {
    pub const fn as_ref(&self) -> OptionOrSkip<&T> {
        match *self {
            OptionOrSkip::Some(ref x) => OptionOrSkip::Some(x),
            OptionOrSkip::None => OptionOrSkip::None,
            OptionOrSkip::Skip => OptionOrSkip::Skip,
        }
    }
    pub fn unwrap(self) -> T {
        match self {
            OptionOrSkip::Some(val) => val,
            OptionOrSkip::None => panic!("called `OptionOrSkip::::unwrap()` on a `None` value"),
            OptionOrSkip::Skip => panic!("called `OptionOrSkip::::unwrap()` on a `Skip` value"),
        }
    }
}

// https://crates.io/crates/typed-builder
// (or https://docs.rs/derive_builder/latest/derive_builder/)
// Somehow ndarray can do this: 	Array::zeros((3, 4, 5).f())
//       see https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html
#[derive(Clone, Debug, Builder)]
#[builder(build_fn(private, name = "build_no_file_check", error = "BedErrorPlus"))]
pub struct Bed {
    // https://stackoverflow.com/questions/32730714/what-is-the-right-way-to-store-an-immutable-path-in-a-struct
    // don't emit a setter, but keep the field declaration on the builder
    #[builder(setter(custom))]
    pub path: PathBuf, // !!!cmk later always clone?

    // !!!cmk 0 play with having an enum and extra set methods
    #[builder(default = "Allele::A1")]
    allele_to_count: Allele,

    // !!!cmk 0 play with having an enum and extra set methods
    #[builder(default = "true")]
    do_format_check: bool,

    #[builder(default, setter(strip_option))]
    iid_count: Option<usize>,

    #[builder(default, setter(strip_option))]
    sid_count: Option<usize>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    iid: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    sid: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    chromosome: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    allele_1: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    allele_2: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    sex: OptionOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    fid: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    father: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    mother: OptionOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    cm_position: OptionOrSkip<nd::Array1<f32>>,

    #[builder(setter(custom))]
    #[builder(default = "OptionOrSkip::None")]
    bp_position: OptionOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,
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
            path: Some(path.as_ref().into()),
            allele_to_count: None,
            do_format_check: None,
            iid_count: None,
            sid_count: None,
            fid: None,
            iid: None,
            father: None,
            mother: None,
            sid: None,
            chromosome: None,
            allele_1: None,
            allele_2: None,
            sex: None,
            cm_position: None,
            bp_position: None,
            fam_path: None,
            bim_path: None,
        }
    }

    // !!!cmk later play with aliasing this as bed_reader::Result<T>
    pub fn build(&self) -> Result<Bed, BedErrorPlus> {
        let bed = self.build_no_file_check()?;

        if bed.do_format_check {
            // !!!cmk later similar code elsewhere
            let mut buf_reader = BufReader::new(File::open(&bed.path)?);
            let mut bytes_vector: Vec<u8> = vec![0; CB_HEADER_USIZE];
            buf_reader.read_exact(&mut bytes_vector)?;
            if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
                return Err(
                    BedError::IllFormed(PathBuf::from(&bed.path).display().to_string()).into(),
                );
            }
        }

        Ok(bed)
    }

    // !!!cmk 0 is it confusing to use 'skip' for both format_check and metadata?
    pub fn skip_format_check(mut self) -> Self {
        self.do_format_check = Some(false);
        self
    }

    pub fn iid_skip(mut self) -> Self {
        self.iid = Some(OptionOrSkip::Skip);
        self
    }

    pub fn sid_skip(mut self) -> Self {
        self.sid = Some(OptionOrSkip::Skip);
        self
    }

    pub fn chromosome_skip(mut self) -> Self {
        self.chromosome = Some(OptionOrSkip::Skip);
        self
    }

    pub fn allele_1_skip(mut self) -> Self {
        self.allele_1 = Some(OptionOrSkip::Skip);
        self
    }

    pub fn allele_2_skip(mut self) -> Self {
        self.allele_2 = Some(OptionOrSkip::Skip);
        self
    }

    pub fn fid_skip(mut self) -> Self {
        self.fid = Some(OptionOrSkip::Skip);
        self
    }
    pub fn father_skip(mut self) -> Self {
        self.father = Some(OptionOrSkip::Skip);
        self
    }
    pub fn mother_skip(mut self) -> Self {
        self.mother = Some(OptionOrSkip::Skip);
        self
    }

    pub fn sex_skip(mut self) -> Self {
        self.sex = Some(OptionOrSkip::Skip);
        self
    }

    pub fn cm_position_skip(mut self) -> Self {
        self.cm_position = Some(OptionOrSkip::Skip);
        self
    }

    pub fn bp_position_skip(mut self) -> Self {
        self.bp_position = Some(OptionOrSkip::Skip);
        self
    }

    // https://stackoverflow.com/questions/38183551/concisely-initializing-a-vector-of-strings
    // https://stackoverflow.com/questions/65250496/how-to-convert-intoiteratoritem-asrefstr-to-iteratoritem-str-in-rust
    pub fn iid<I, T>(mut self, iid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            iid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.iid = Some(OptionOrSkip::Some(new_string_array));
        self
    }

    pub fn sid<I, T>(mut self, sid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            sid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.sid = Some(OptionOrSkip::Some(new_string_array));
        self
    }

    // !!!cmk 0 can we set new fathers, etc and skip them?
    pub fn chromosome<I, T>(mut self, chromosome: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> = chromosome
            .into_iter()
            .map(|s| s.as_ref().to_string())
            .collect();
        self.chromosome = Some(OptionOrSkip::Some(new_string_array));
        self
    }

    pub fn fam_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.fam_path = Some(Some(path.as_ref().into()));
        self
    }

    pub fn bim_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.bim_path = Some(Some(path.as_ref().into()));
        self
    }
}

enum FamOrBim {
    Fam,
    Bim,
}

impl Bed {
    // !!!cmk later
    // fam_filepath: Union[str, Path] = None,
    // bim_filepath: Union[str, Path] = None,

    // !!!cmk later is this how you do lazy accessors?
    // !!!cmk later should this be "try_get_..." or just "iid_count" or as is
    pub fn iid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let (iid_count, sid_count) = counts(&self.path)?;
            self.iid_count = Some(iid_count);
            self.sid_count = Some(sid_count);
            Ok(iid_count)
        }
    }
    pub fn sid_count(&mut self) -> Result<usize, BedErrorPlus> {
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
        fam_or_bim: FamOrBim,
        field_index: usize,
    ) -> Result<nd::Array1<String>, BedErrorPlus> {
        let path_buf = match fam_or_bim {
            FamOrBim::Fam => {
                if let Some(fam_path) = &self.fam_path {
                    fam_path.clone()
                } else {
                    self.path.with_extension("fam")
                }
            }
            FamOrBim::Bim => {
                if let Some(bim_path) = &self.bim_path {
                    bim_path.clone()
                } else {
                    self.path.with_extension("bim")
                }
            }
        };

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

    pub fn iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.iid {
            OptionOrSkip::Some(ref iid) => Ok(iid),
            OptionOrSkip::None => {
                let iid = self.read_fam_or_bim(FamOrBim::Fam, 1)?;
                self.iid = OptionOrSkip::Some(iid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.iid.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("iid".to_string()).into()),
        }
    }

    pub fn sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.sid {
            OptionOrSkip::Some(ref sid) => Ok(sid),
            OptionOrSkip::None => {
                let sid = self.read_fam_or_bim(FamOrBim::Bim, 2)?;
                self.sid = OptionOrSkip::Some(sid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.sid.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("sid".to_string()).into()),
        }
    }

    pub fn chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.chromosome {
            OptionOrSkip::Some(ref chromosome) => Ok(chromosome),
            OptionOrSkip::None => {
                let chromosome = self.read_fam_or_bim(FamOrBim::Bim, 0)?;
                self.chromosome = OptionOrSkip::Some(chromosome);
                // This unwrap is safe because we just created 'Some'
                Ok(self.chromosome.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("chromosome".to_string()).into()),
        }
    }

    pub fn allele_1(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_1 {
            OptionOrSkip::Some(ref allele_1) => Ok(allele_1),
            OptionOrSkip::None => {
                let allele_1 = self.read_fam_or_bim(FamOrBim::Bim, 4)?;
                self.allele_1 = OptionOrSkip::Some(allele_1);
                // This unwrap is safe because we just created 'Some'
                Ok(self.allele_1.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("allele_1".to_string()).into()),
        }
    }

    pub fn allele_2(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_2 {
            OptionOrSkip::Some(ref allele_2) => Ok(allele_2),
            OptionOrSkip::None => {
                let allele_2 = self.read_fam_or_bim(FamOrBim::Bim, 5)?;
                self.allele_2 = OptionOrSkip::Some(allele_2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.allele_2.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("allele_2".to_string()).into()),
        }
    }

    pub fn fid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.fid {
            OptionOrSkip::Some(ref fid) => Ok(fid),
            OptionOrSkip::None => {
                let fid = self.read_fam_or_bim(FamOrBim::Fam, 0)?;
                self.fid = OptionOrSkip::Some(fid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.fid.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("fid".to_string()).into()),
        }
    }

    pub fn father(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.father {
            OptionOrSkip::Some(ref father) => Ok(father),
            OptionOrSkip::None => {
                let father = self.read_fam_or_bim(FamOrBim::Fam, 2)?;
                self.father = OptionOrSkip::Some(father);
                // This unwrap is safe because we just created 'Some'
                Ok(self.father.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("father".to_string()).into()),
        }
    }

    pub fn mother(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.mother {
            OptionOrSkip::Some(ref mother) => Ok(mother),
            OptionOrSkip::None => {
                let mother = self.read_fam_or_bim(FamOrBim::Fam, 3)?;
                self.mother = OptionOrSkip::Some(mother);
                // This unwrap is safe because we just created 'Some'
                Ok(self.mother.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("mother".to_string()).into()),
        }
    }

    pub fn sex(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.sex {
            OptionOrSkip::Some(ref sex) => Ok(sex),
            OptionOrSkip::None => {
                let sex: Result<nd::Array1<i32>, _> = self
                    .read_fam_or_bim(FamOrBim::Fam, 4)?
                    .into_iter()
                    .map(|s| match s.parse::<i32>() {
                        Err(_) => Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                            "!!!cmk0".to_string(),
                        ))),
                        Ok(i) => Ok(i),
                    })
                    .collect();
                let sex2 = sex?;
                self.sex = OptionOrSkip::Some(sex2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.sex.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("sex".to_string()).into()),
        }
    }

    pub fn cm_position(&mut self) -> Result<&nd::Array1<f32>, BedErrorPlus> {
        match self.cm_position {
            OptionOrSkip::Some(ref cm_position) => Ok(cm_position),
            OptionOrSkip::None => {
                let cm_position: Result<nd::Array1<f32>, _> = self
                    .read_fam_or_bim(FamOrBim::Bim, 2)?
                    .into_iter()
                    .map(|s| match s.parse::<f32>() {
                        Err(err) => {
                            println!("!!!cmk{}", s);
                            Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                                err.to_string(),
                            )))
                        }
                        Ok(i) => Ok(i),
                    })
                    .collect();
                let cm_position2 = cm_position?;
                self.cm_position = OptionOrSkip::Some(cm_position2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.cm_position.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("cm_position".to_string()).into()),
        }
    }

    pub fn bp_position(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.bp_position {
            OptionOrSkip::Some(ref bp_position) => Ok(bp_position),
            OptionOrSkip::None => {
                let bp_position: Result<nd::Array1<i32>, _> = self
                    .read_fam_or_bim(FamOrBim::Bim, 3)?
                    .into_iter()
                    .map(|s| match s.parse::<i32>() {
                        Err(_) => Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                            "!!!cmk0".to_string(),
                        ))),
                        Ok(i) => Ok(i),
                    })
                    .collect();
                let bp_position2 = bp_position?;
                self.bp_position = OptionOrSkip::Some(bp_position2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.bp_position.as_ref().unwrap())
            }
            OptionOrSkip::Skip => Err(BedError::PropertySkipped("bp_position".to_string()).into()),
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
        let iid_count = self.iid_count()?;
        let sid_count = self.sid_count()?;

        let num_threads = if let Some(num_threads) = read_options.num_threads {
            num_threads
        } else {
            if let Ok(num_threads) = env::var("BED_READER_NUM_THREADS") {
                num_threads.parse::<usize>()?
            } else if let Ok(num_threads) = env::var("NUM_THREADS") {
                num_threads.parse::<usize>()?
            } else {
                0
            }
        };

        let iid_index = read_options.iid_index.to_vec(iid_count);
        let sid_index = read_options.sid_index.to_vec(sid_count);

        let shape = ShapeBuilder::set_f((iid_index.len(), sid_index.len()), read_options.is_f);
        let mut val = nd::Array2::<TOut>::default(shape);

        read_no_alloc(
            &self.path,
            iid_count,
            sid_count,
            self.allele_to_count == Allele::A1,
            &iid_index,
            &sid_index,
            read_options.missing_value,
            num_threads,
            &mut val.view_mut(),
        )?;

        Ok(val)
    }
}

impl Index {
    // !!!cmk later test every case
    // We can't define a 'From' because we want to add count at the last moment.
    // Would be nice to not always allocate a new vec, maybe with Rc<[T]>?
    // Even better would be to support an iterator from Index (an enum with fields).
    pub fn to_vec(&self, count: usize) -> Vec<usize> {
        match self {
            Index::None => (0..count).collect(),
            Index::Vec(vec) => vec.to_vec(),
            Index::NDArrayBool(nd_array_bool) => {
                // !!!cmk later check that bool_index.len() == iid_count
                nd_array_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i)
                    .collect()
            }
            // !!!cmk later can we implement this without two allocations?
            Index::NDSliceInfo(nd_slice_info) => {
                // https://docs.rs/ndarray/0.15.4/ndarray/struct.ArrayBase.html#slicing
                let mut array: nd::Array1<usize> = (0..count).collect();
                array.slice_collapse(nd_slice_info);
                array.to_vec()
            }
            Index::Range(range) => range.clone().collect(),
            Index::NDArray(nd_array) => nd_array.to_vec(),
            Index::One(one) => vec![*one],
            Index::VecBool(vec_bool) => {
                // !!!cmk later similar code elsewhere
                // !!!cmk later check that vec_bool.len() == iid_count
                vec_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i)
                    .collect()
            }
        }
    }
}

pub(crate) type SliceInfo1 =
    nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>;

// Could implement an enumerator, but it is complex and requires a 'match' on each next()
//     https://stackoverflow.com/questions/65272613/how-to-implement-intoiterator-for-an-enum-of-iterable-variants
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
#[derive(Debug, Clone, Builder)]
pub struct ReadOptions<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing> {
    #[builder(default = "TOut::missing()")]
    missing_value: TOut,

    #[builder(default = "Index::None")]
    #[builder(setter(into))]
    iid_index: Index,

    #[builder(default = "Index::None")]
    #[builder(setter(into))]
    sid_index: Index,

    // !!!cmk 0 play with having an enum and extra set methods
    #[builder(default = "true")]
    is_f: bool, // !!!cmk later use enum or .f()

    #[builder(default, setter(strip_option))]
    pub num_threads: Option<usize>,
}

impl<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing + Clone> ReadOptions<TOut> {
    pub fn builder() -> ReadOptionsBuilder<TOut> {
        ReadOptionsBuilder::default()
    }
}

// !!!cmk 0 alias "Copy + Default + From<i8> + Debug + Sync + Send + Missing + Clone"?
impl<TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing + Clone>
    ReadOptionsBuilder<TOut>
{
    pub fn read(&self, bed: &mut Bed) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_option = self.build()?;
        bed.read_with_options(read_option)
    }

    pub fn c(&mut self) -> &Self {
        self.is_f(false);
        self
    }
}

impl ReadOptionsBuilder<i8> {
    pub fn i8(&self) -> &Self {
        self
    }
}

impl ReadOptionsBuilder<f32> {
    pub fn f32(&self) -> &Self {
        self
    }
}

impl ReadOptionsBuilder<f64> {
    pub fn f64(&self) -> &Self {
        self
    }
}

// !!!cmk later could a macro likes be nice?
// #[macro_export]
// macro_rules! read {
//     ($bed:expr) => {
//         $bed.read()
//     };
//     ($bed:expr, $option: expr) => {
//         ReadOptions::builder().$option.read($bed)
//     };
// }
