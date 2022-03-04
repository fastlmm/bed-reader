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
use std::io::Write;
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    ops::Range,
    path::{Path, PathBuf},
};

use derive_builder::Builder;
// !!! might want to use this instead use typed_builder::TypedBuilder;

use crate::{
    count_lines, read_no_alloc, write_val, BedError, BedErrorPlus, BedVal, BED_FILE_MAGIC1,
    BED_FILE_MAGIC2, CB_HEADER_USIZE,
};

#[derive(Debug, Clone)]
pub enum LazyOrSkip<T> {
    Lazy,
    Skip,
    Some(T),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Skippable<T> {
    Some(T),
    Skip,
}

impl<T> Skippable<T> {
    pub fn unwrap(self) -> T {
        match self {
            Skippable::Some(some) => some,
            Skippable::Skip => {
                todo!("Skippable::Skip")
            }
        }
    }
}

impl<T> LazyOrSkip<T> {
    pub const fn as_ref(&self) -> LazyOrSkip<&T> {
        match *self {
            LazyOrSkip::Some(ref x) => LazyOrSkip::Some(x),
            LazyOrSkip::Lazy => LazyOrSkip::Lazy,
            LazyOrSkip::Skip => LazyOrSkip::Skip,
        }
    }
    pub fn unwrap(self) -> T {
        match self {
            LazyOrSkip::Some(val) => val,
            LazyOrSkip::Lazy => panic!("called `LazyOrSkip::::unwrap()` on a `Lazy` value"),
            LazyOrSkip::Skip => panic!("called `LazyOrSkip::::unwrap()` on a `Skip` value"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Metadata {
    pub fid: Skippable<nd::Array1<String>>,
    pub iid: Skippable<nd::Array1<String>>,
    pub father: Skippable<nd::Array1<String>>,
    pub mother: Skippable<nd::Array1<String>>,
    pub sex: Skippable<nd::Array1<i32>>,
    pub pheno: Skippable<nd::Array1<String>>,

    pub chromosome: Skippable<nd::Array1<String>>,
    pub sid: Skippable<nd::Array1<String>>,
    pub cm_position: Skippable<nd::Array1<f32>>,
    pub bp_position: Skippable<nd::Array1<i32>>,
    pub allele_1: Skippable<nd::Array1<String>>,
    pub allele_2: Skippable<nd::Array1<String>>,
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

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,

    #[builder(default = "true")]
    is_a1_counted: bool,

    #[builder(default = "true")]
    is_checked_early: bool,

    #[builder(default, setter(strip_option))]
    iid_count: Option<usize>,

    #[builder(default, setter(strip_option))]
    sid_count: Option<usize>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    fid: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    iid: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    father: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    mother: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    sex: LazyOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    pheno: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    chromosome: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    sid: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    cm_position: LazyOrSkip<nd::Array1<f32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    bp_position: LazyOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_1: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_2: LazyOrSkip<nd::Array1<String>>,

    #[builder(private)]
    fam_file: PathBuf,

    #[builder(private)]
    bim_file: PathBuf,
}

impl BedBuilder {
    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        let path: PathBuf = path.as_ref().into();

        let fam_file = to_metadata_path(&path, &None, "fam");
        let bim_file = to_metadata_path(&path, &None, "bim");

        Self {
            path: Some(path),
            fam_path: None,
            bim_path: None,

            is_a1_counted: None,
            is_checked_early: None,
            iid_count: None,
            sid_count: None,

            fid: None,
            iid: None,
            father: None,
            mother: None,
            sex: None,
            pheno: None,

            chromosome: None,
            sid: None,
            cm_position: None,
            bp_position: None,
            allele_1: None,
            allele_2: None,

            // !!!cmk later give better name
            fam_file: Some(fam_file),
            bim_file: Some(bim_file),
        }
    }

    // !!!cmk later play with aliasing this as bed_reader::Result<T>
    pub fn build(&self) -> Result<Bed, BedErrorPlus> {
        let mut bed = self.build_no_file_check()?;

        bed.fam_file = to_metadata_path(&bed.path, &bed.fam_path, "fam");
        bed.bim_file = to_metadata_path(&bed.path, &bed.bim_path, "bim");

        if bed.is_checked_early {
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

    pub fn fam_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.fam_path = Some(Some(path.as_ref().into()));
        self
    }

    pub fn bim_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.bim_path = Some(Some(path.as_ref().into()));
        self
    }

    pub fn count_a1(&mut self) -> &mut Self {
        self.is_a1_counted = Some(true);
        self
    }

    pub fn count_a2(&mut self) -> &mut Self {
        self.is_a1_counted = Some(false);
        self
    }
    pub fn skip_early_check(mut self) -> Self {
        self.is_checked_early = Some(false);
        self
    }

    pub fn skip_iid(&mut self) -> &Self {
        self.iid = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_fid(mut self) -> Self {
        self.fid = Some(LazyOrSkip::Skip);
        self
    }
    pub fn skip_father(mut self) -> Self {
        self.father = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_mother(mut self) -> Self {
        self.mother = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_sex(mut self) -> Self {
        self.sex = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_pheno(mut self) -> Self {
        self.pheno = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_chromosome(mut self) -> Self {
        self.chromosome = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_sid(mut self) -> Self {
        self.sid = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_cm_position(mut self) -> Self {
        self.cm_position = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_bp_position(mut self) -> Self {
        self.bp_position = Some(LazyOrSkip::Skip);
        self
    }
    pub fn skip_allele_1(mut self) -> Self {
        self.allele_1 = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_allele_2(mut self) -> Self {
        self.allele_2 = Some(LazyOrSkip::Skip);
        self
    }

    // https://stackoverflow.com/questions/38183551/concisely-initializing-a-vector-of-strings
    // https://stackoverflow.com/questions/65250496/how-to-convert-intoiteratoritem-asrefstr-to-iteratoritem-str-in-rust
    pub fn fid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, fid: I) -> Self {
        self.fid = Some(LazyOrSkip::Some(
            fid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    pub fn iid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, iid: I) -> Self {
        self.iid = Some(LazyOrSkip::Some(
            iid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }
    pub fn father<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, father: I) -> Self {
        self.father = Some(LazyOrSkip::Some(
            father.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }
    pub fn mother<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, mother: I) -> Self {
        self.mother = Some(LazyOrSkip::Some(
            mother.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn sex<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, sex: I) -> Self {
        self.sex = Some(LazyOrSkip::Some(
            sex.into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    pub fn pheno<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, pheno: I) -> Self {
        self.pheno = Some(LazyOrSkip::Some(
            pheno.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    pub fn chromosome<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, chromosome: I) -> Self {
        self.chromosome = Some(LazyOrSkip::Some(
            chromosome
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }

    pub fn sid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, sid: I) -> Self {
        self.sid = Some(LazyOrSkip::Some(
            sid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn cm_position<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, cm_position: I) -> Self {
        self.cm_position = Some(LazyOrSkip::Some(
            cm_position
                .into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn bp_position<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, bp_position: I) -> Self {
        self.bp_position = Some(LazyOrSkip::Some(
            bp_position
                .into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    pub fn allele_1<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, allele_1: I) -> Self {
        self.allele_1 = Some(LazyOrSkip::Some(
            allele_1
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }

    pub fn allele_2<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, allele_2: I) -> Self {
        self.allele_2 = Some(LazyOrSkip::Some(
            allele_2
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }
}

fn to_metadata_path(
    bed_path: &PathBuf,
    metadata_path: &Option<PathBuf>,
    extension: &str,
) -> PathBuf {
    if let Some(metadata_path) = metadata_path {
        metadata_path.clone()
    } else {
        bed_path.with_extension(extension)
    }
}

fn to_skippable<T: Clone>(lazy_or_skip: &LazyOrSkip<T>) -> Skippable<T> {
    match lazy_or_skip {
        LazyOrSkip::Lazy => {
            todo!("!!!cmk later")
        }
        LazyOrSkip::Skip => Skippable::Skip,
        LazyOrSkip::Some(some) => Skippable::Some(some.clone()),
    }
}

impl Bed {
    pub fn builder<P: AsRef<Path>>(path: P) -> BedBuilder {
        BedBuilder::new(path)
    }

    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedErrorPlus> {
        Bed::builder(path).build()
    }

    // !!!cmk later is this how you do lazy accessors?
    // !!!cmk later should this be "try_get_..." or just "iid_count" or as is
    pub fn iid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let iid_count = count_lines(&self.fam_file)?;
            self.iid_count = Some(iid_count);
            Ok(iid_count)
        }
    }
    pub fn sid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let sid_count = count_lines(&self.bim_file)?;
            self.sid_count = Some(sid_count);
            Ok(sid_count)
        }
    }

    /// !!!cmk later don't re-read for every column

    pub fn metadata(&mut self) -> Result<Metadata, BedErrorPlus> {
        let _ = self.fid();
        let _ = self.iid();
        let _ = self.father();
        let _ = self.mother();
        let _ = self.sex();
        let _ = self.pheno();

        let _ = self.chromosome();
        let _ = self.sid();
        let _ = self.cm_position();
        let _ = self.bp_position();
        let _ = self.allele_1();
        let _ = self.allele_2();

        let metadata = Metadata {
            fid: to_skippable(&self.fid),
            iid: to_skippable(&self.iid),
            father: to_skippable(&self.father),
            mother: to_skippable(&self.mother),
            sex: to_skippable(&self.sex),
            pheno: to_skippable(&self.pheno),

            chromosome: to_skippable(&self.chromosome),
            sid: to_skippable(&self.sid),
            cm_position: to_skippable(&self.cm_position),
            bp_position: to_skippable(&self.bp_position),
            allele_1: to_skippable(&self.allele_1),
            allele_2: to_skippable(&self.allele_2),
        };
        Ok(metadata)
    }

    pub fn fid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.fid, &self.fam_file, 0, "fid")
    }
    pub fn iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.iid, &self.fam_file, 1, "iid")
    }
    pub fn father(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.father, &self.fam_file, 2, "father")
    }
    pub fn mother(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.mother, &self.fam_file, 3, "mother")
    }
    pub fn sex(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        unlazy(&mut self.sex, &self.fam_file, 4, "sex")
    }
    pub fn pheno(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.pheno, &self.fam_file, 5, "pheno")
    }

    pub fn chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.chromosome, &self.bim_file, 0, "chromosome")
    }
    pub fn sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.sid, &self.bim_file, 1, "sid")
    }
    pub fn cm_position(&mut self) -> Result<&nd::Array1<f32>, BedErrorPlus> {
        unlazy(&mut self.cm_position, &self.bim_file, 2, "cm_position")
    }
    pub fn bp_position(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        unlazy(&mut self.bp_position, &self.bim_file, 3, "bp_position")
    }
    pub fn allele_1(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.allele_1, &self.bim_file, 4, "allele_1")
    }
    pub fn allele_2(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        unlazy(&mut self.allele_2, &self.bim_file, 5, "allele_2")
    }

    // !!!cmk later rename TOut to TVal
    pub fn read<TOut: BedVal>(&mut self) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_options = ReadOptions::<TOut>::builder().build()?;
        self.read_with_options(read_options)
    }

    pub(crate) fn read_with_options<TOut: BedVal>(
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
            self.is_a1_counted,
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
pub struct ReadOptions<TOut: BedVal> {
    #[builder(default = "TOut::missing()")]
    missing_value: TOut,

    #[builder(default = "Index::None")]
    #[builder(setter(into))]
    iid_index: Index,

    #[builder(default = "Index::None")]
    #[builder(setter(into))]
    sid_index: Index,

    #[builder(default = "true")]
    is_f: bool,

    #[builder(default, setter(strip_option))]
    pub num_threads: Option<usize>,
}

impl<TOut: BedVal> ReadOptions<TOut> {
    pub fn builder() -> ReadOptionsBuilder<TOut> {
        ReadOptionsBuilder::default()
    }
}

impl<TOut: BedVal> ReadOptionsBuilder<TOut> {
    pub fn read(&self, bed: &mut Bed) -> Result<nd::Array2<TOut>, BedErrorPlus> {
        let read_option = self.build()?;
        bed.read_with_options(read_option)
    }

    pub fn f(&mut self) -> &Self {
        self.is_f(true);
        self
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

#[derive(Clone, Debug, Builder)]
#[builder(build_fn(private, name = "write_no_file_check", error = "BedErrorPlus"))]
pub struct WriteOptions {
    #[builder(setter(custom))]
    pub path: PathBuf, // !!!cmk later always clone?

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    fid: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    iid: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    father: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    mother: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    sex: Skippable<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    pheno: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    chromosome: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    sid: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    cm_position: Skippable<nd::Array1<f32>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    bp_position: Skippable<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    allele_1: Skippable<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "Skippable::Skip")]
    allele_2: Skippable<nd::Array1<String>>,
}

impl WriteOptions {
    pub fn builder<P: AsRef<Path>>(path: P) -> WriteOptionsBuilder {
        WriteOptionsBuilder::new(path)
    }

    pub fn fid(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.fid {
            Skippable::Some(ref fid) => Ok(fid),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("fid".to_string()).into()),
        }
    }
    pub fn iid(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.iid {
            Skippable::Some(ref iid) => Ok(iid),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("iid".to_string()).into()),
        }
    }
    pub fn father(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.father {
            Skippable::Some(ref father) => Ok(father),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("father".to_string()).into()),
        }
    }
    pub fn mother(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.mother {
            Skippable::Some(ref mother) => Ok(mother),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("mother".to_string()).into()),
        }
    }
    pub fn sex(&self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.sex {
            Skippable::Some(ref sex) => Ok(sex),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("sex".to_string()).into()),
        }
    }
    pub fn pheno(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.pheno {
            Skippable::Some(ref pheno) => Ok(pheno),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("pheno".to_string()).into()),
        }
    }

    pub fn chromosome(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.chromosome {
            Skippable::Some(ref chromosome) => Ok(chromosome),
            Skippable::Skip => {
                Err(BedError::CannotUseSkippedMetadata("chromosome".to_string()).into())
            }
        }
    }
    pub fn sid(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.sid {
            Skippable::Some(ref sid) => Ok(sid),
            Skippable::Skip => Err(BedError::CannotUseSkippedMetadata("sid".to_string()).into()),
        }
    }
    pub fn cm_position(&self) -> Result<&nd::Array1<f32>, BedErrorPlus> {
        match self.cm_position {
            Skippable::Some(ref cm_position) => Ok(cm_position),
            Skippable::Skip => {
                Err(BedError::CannotUseSkippedMetadata("cm_position".to_string()).into())
            }
        }
    }

    pub fn bp_position(&self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.bp_position {
            Skippable::Some(ref bp_position) => Ok(bp_position),
            Skippable::Skip => {
                Err(BedError::CannotUseSkippedMetadata("bp_position".to_string()).into())
            }
        }
    }
    pub fn allele_1(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_1 {
            Skippable::Some(ref allele_1) => Ok(allele_1),
            Skippable::Skip => {
                Err(BedError::CannotUseSkippedMetadata("allele_1".to_string()).into())
            }
        }
    }
    pub fn allele_2(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_2 {
            Skippable::Some(ref allele_2) => Ok(allele_2),
            Skippable::Skip => {
                Err(BedError::CannotUseSkippedMetadata("allele_2".to_string()).into())
            }
        }
    }
}

impl WriteOptionsBuilder {
    // !!! cmk later just use the default builder?
    pub fn build(&self) -> Result<WriteOptions, BedErrorPlus> {
        let write_options = self.write_no_file_check()?;

        Ok(write_options)
    }

    // !!!cmk later should check that metadata agrees with val size
    // !!!cmk later maybe use the default builder?
    pub fn write<TVal: BedVal>(&self, val: &nd::Array2<TVal>) -> Result<(), BedErrorPlus> {
        let write_options = self.build()?;
        write_with_options(val, &write_options)?;

        Ok(())
    }

    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        Self {
            path: Some(path.as_ref().into()),
            fam_path: None,
            bim_path: None,

            fid: None,
            iid: None,
            father: None,
            mother: None,
            sex: None,
            pheno: None,

            chromosome: None,
            sid: None,
            cm_position: None,
            bp_position: None,
            allele_1: None,
            allele_2: None,
        }
    }

    pub fn fam_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.fam_path = Some(Some(path.as_ref().into()));
        self
    }

    pub fn bim_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.bim_path = Some(Some(path.as_ref().into()));
        self
    }

    pub fn metadata(mut self, metadata: &Metadata) -> Self {
        if let Skippable::Some(iid) = &metadata.iid {
            self.iid = Some(Skippable::Some(iid.clone()));
        }
        if let Skippable::Some(sid) = &metadata.sid {
            self.sid = Some(Skippable::Some(sid.clone()));
        }
        self
    }

    pub fn fid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, fid: I) -> Self {
        self.fid = Some(Skippable::Some(
            fid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    pub fn iid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, iid: I) -> Self {
        self.iid = Some(Skippable::Some(
            iid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }
    pub fn father<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, father: I) -> Self {
        self.father = Some(Skippable::Some(
            father.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }
    pub fn mother<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, mother: I) -> Self {
        self.mother = Some(Skippable::Some(
            mother.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn sex<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, sex: I) -> Self {
        self.sex = Some(Skippable::Some(
            sex.into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    pub fn pheno<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, pheno: I) -> Self {
        self.pheno = Some(Skippable::Some(
            pheno.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    pub fn chromosome<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, chromosome: I) -> Self {
        self.chromosome = Some(Skippable::Some(
            chromosome
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }

    pub fn sid<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, sid: I) -> Self {
        self.sid = Some(Skippable::Some(
            sid.into_iter().map(|s| s.as_ref().to_string()).collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn cm_position<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, cm_position: I) -> Self {
        self.cm_position = Some(Skippable::Some(
            cm_position
                .into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    // !!!cmk later fix the unwraps
    pub fn bp_position<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, bp_position: I) -> Self {
        self.bp_position = Some(Skippable::Some(
            bp_position
                .into_iter()
                .map(|s| s.as_ref().parse().unwrap())
                .collect(),
        ));
        self
    }

    pub fn allele_1<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, allele_1: I) -> Self {
        self.allele_1 = Some(Skippable::Some(
            allele_1
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }

    pub fn allele_2<I: IntoIterator<Item = T>, T: AsRef<str>>(mut self, allele_2: I) -> Self {
        self.allele_2 = Some(Skippable::Some(
            allele_2
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        ));
        self
    }
}

pub fn write<TVal: BedVal>(val: &nd::Array2<TVal>, path: &Path) -> Result<(), BedErrorPlus> {
    WriteOptions::builder(path).write(val)
}

pub fn write_with_options<TVal: BedVal>(
    val: &nd::Array2<TVal>,
    write_options: &WriteOptions,
) -> Result<(), BedErrorPlus> {
    // !!!cmk later can this be done in one step??
    let shape = val.shape();
    let iid_count = shape[0];
    let sid_count = shape[1];
    // !!!cmk later if something goes wrong, clean up the files?
    let path = &write_options.path;

    // !!!cmk 0 set is_a1_count
    // !!!cmk 0 set missing
    // !!!cmk set num_threads
    write_val(path, &val.view(), true, TVal::missing(), 0)?;

    let fam_path = to_metadata_path(path, &write_options.fam_path, "fam");
    let bim_path = to_metadata_path(path, &write_options.bim_path, "bim");

    // !!!cmk later can we avoid always allocating the default, maybe with cow?
    let iid_default: nd::Array1<String> = (0..iid_count).map(|i| format!("iid{}", i + 1)).collect();
    let iid = match write_options.iid {
        // !!!cmk later check that length agrees with shape
        Skippable::Some(ref iid) => iid,
        Skippable::Skip => &iid_default,
    };

    // !!!cmk later can we avoid always allocating the default, maybe with cow?
    let sid_default: nd::Array1<String> = (0..sid_count).map(|i| format!("sid{}", i + 1)).collect();
    let sid = match write_options.sid {
        Skippable::Some(ref sid) => sid,
        Skippable::Skip => &sid_default,
    };

    println!("{:?}", iid);
    for ix in iid {
        println!("{:?}", ix);
    }

    // !!!cmk 0 finish this code
    let iid_zeros = nd::Array1::from_elem(iid_count, "0");
    let sid_zero = nd::Array1::from_elem(sid_count, "0");

    let file = File::create(fam_path)?;
    let mut writer = BufWriter::new(file);
    nd::azip!(iid).for_each(|iid_zero_x| {
        println!("{:?}", *iid_zero_x);
    });
    nd::azip!(&iid_zeros)
        .and(iid)
        .and(&iid_zeros)
        .and(&iid_zeros)
        .and(&iid_zeros)
        .and(&iid_zeros)
        .for_each(|fid, iid, father, mother, sex, pheno| {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                *fid, *iid, *father, *mother, *sex, *pheno
            )
            .unwrap(); // !!!cmk later throw error
        });

    let sid_zeros = nd::Array1::from_elem(sid_count, "0");
    let allele_1_default = nd::Array1::from_elem(sid_count, "A1");
    let allele_2_default = nd::Array1::from_elem(sid_count, "A2");

    let file = File::create(bim_path)?;
    let mut writer = BufWriter::new(file);
    nd::azip!(&sid_zeros)
        .and(sid)
        .and(&sid_zeros)
        .and(&sid_zeros)
        .and(&allele_1_default)
        .and(&allele_2_default)
        .for_each(
            |chromosome, sid, cm_position, bp_position, allele_1, allele_2| {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    *chromosome, *sid, *cm_position, *bp_position, *allele_1, *allele_2
                )
                .unwrap(); // !!!cmk later throw error
            },
        );
    Ok(())
}

trait FromStringArray<T> {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<Self>, BedErrorPlus>
    where
        Self: Sized;
}

impl FromStringArray<String> for String {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<String>, BedErrorPlus> {
        Ok(string_array)
    }
}

// !!!cmk 0 fix these so they return any errors
impl FromStringArray<f64> for f64 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<f64>, BedErrorPlus> {
        Ok(string_array.map(|s| s.parse::<f64>().unwrap()))
    }
}
impl FromStringArray<f32> for f32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<f32>, BedErrorPlus> {
        Ok(string_array.map(|s| s.parse::<f32>().unwrap()))
    }
}
impl FromStringArray<i32> for i32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<i32>, BedErrorPlus> {
        Ok(string_array.map(|s| s.parse::<i32>().unwrap()))
    }
}
fn unlazy<'a, T: FromStringArray<T>>(
    field: &'a mut LazyOrSkip<nd::Array1<T>>,
    path_buf: &PathBuf,
    field_index: usize,
    name: &str,
) -> Result<&'a nd::Array1<T>, BedErrorPlus> {
    match field {
        LazyOrSkip::Some(ref array) => Ok(array),
        LazyOrSkip::Lazy => {
            let array: nd::Array1<T> =
                T::from_string_array(read_fam_or_bim(path_buf, field_index)?)?;
            *field = LazyOrSkip::Some(array);
            match field {
                LazyOrSkip::Some(ref array) => Ok(array),
                _ => panic!("impossible"),
            }
        }
        LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata(name.to_string()).into()),
    }
}

fn read_fam_or_bim(
    path_buf: &PathBuf,
    field_index: usize,
) -> Result<nd::Array1<String>, BedErrorPlus> {
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
