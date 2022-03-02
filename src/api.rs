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
    counts, read_no_alloc, write, BedError, BedErrorPlus, BedVal, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
    CB_HEADER_USIZE,
};

// !!!cmk 0 replace this with Option<Skippable>
#[derive(Debug, Clone)]
pub enum LazyOrSkip<T> {
    Lazy,
    Skip,
    Some(T),
}

#[derive(Debug, Clone)]
pub enum Skippable<T> {
    Some(T),
    Skip,
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

// impl<T> Skippable<T> {
//     pub const fn as_ref(&self) -> Skippable<&T> {
//         match *self {
//             Skippable::Some(ref x) => Skippable::Some(x),
//             Skippable::Skip => Skippable::Skip,
//         }
//     }
//     pub fn unwrap(self) -> T {
//         match self {
//             Skippable::Some(val) => val,
//             Skippable::Skip => panic!("called `Skippable::::unwrap()` on a `Skip` value"),
//         }
//     }
// }

#[derive(Clone, Debug)]
pub struct Metadata {
    pub iid: LazyOrSkip<nd::Array1<String>>,

    pub sid: LazyOrSkip<nd::Array1<String>>,
}

impl Default for Metadata {
    fn default() -> Self {
        Self {
            iid: LazyOrSkip::Lazy,
            sid: LazyOrSkip::Lazy,
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

    #[builder(default = "true")]
    is_a1_counted: bool,

    #[builder(default = "true")]
    is_checked_early: bool,

    #[builder(default, setter(strip_option))]
    iid_count: Option<usize>,

    #[builder(default, setter(strip_option))]
    sid_count: Option<usize>,

    metadata: Metadata,

    // #[builder(setter(custom))]
    // #[builder(default = "LazyOrSkip::Lazy")]
    // iid: LazyOrSkip<nd::Array1<String>>,

    // #[builder(setter(custom))]
    // #[builder(default = "LazyOrSkip::Lazy")]
    // sid: LazyOrSkip<nd::Array1<String>>,
    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    chromosome: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_1: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_2: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    sex: LazyOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    fid: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    father: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    mother: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    cm_position: LazyOrSkip<nd::Array1<f32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    bp_position: LazyOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,
}

impl BedBuilder {
    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        Self {
            path: Some(path.as_ref().into()),
            is_a1_counted: None,
            is_checked_early: None,
            iid_count: None,
            sid_count: None,
            fid: None,
            metadata: Some(Metadata::default()),
            father: None,
            mother: None,
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
        // !!!cmk 0 is that unwrap safe?
        self.metadata.as_mut().unwrap().iid = LazyOrSkip::Skip;
        self
    }

    pub fn skip_sid(mut self) -> Self {
        self.metadata.as_mut().unwrap().sid = LazyOrSkip::Skip;
        self
    }

    pub fn skip_chromosome(mut self) -> Self {
        self.chromosome = Some(LazyOrSkip::Skip);
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

    pub fn skip_cm_position(mut self) -> Self {
        self.cm_position = Some(LazyOrSkip::Skip);
        self
    }

    pub fn skip_bp_position(mut self) -> Self {
        self.bp_position = Some(LazyOrSkip::Skip);
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
        self.metadata.as_mut().unwrap().iid = LazyOrSkip::Some(new_string_array);
        self
    }

    pub fn sid<I, T>(mut self, sid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            sid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.metadata.as_mut().unwrap().sid = LazyOrSkip::Some(new_string_array);
        self
    }

    pub fn chromosome<I, T>(mut self, chromosome: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> = chromosome
            .into_iter()
            .map(|s| s.as_ref().to_string())
            .collect();
        self.chromosome = Some(LazyOrSkip::Some(new_string_array));
        self
    }

    pub fn fid<I, T>(mut self, fid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            fid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.fid = Some(LazyOrSkip::Some(new_string_array));
        self
    }
    pub fn father<I, T>(mut self, father: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            father.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.father = Some(LazyOrSkip::Some(new_string_array));
        self
    }
    pub fn mother<I, T>(mut self, mother: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            mother.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.mother = Some(LazyOrSkip::Some(new_string_array));
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

// fn skip_to_none(
//     result: Result<&nd::Array1<String>, BedErrorPlus>,
// ) -> Result<Option<nd::Array1<String>>, BedErrorPlus> {
//     match result {
//         // !!!cmk 0 understand to_owned. Is there a copy?
//         Ok(array) => Ok(Some(array.to_owned())),
//         Err(BedErrorPlus::BedError(BedError::CannotUseSkippedMetadata(_))) => Ok(None),
//         Err(e) => Err(e),
//     }
// }

impl Bed {
    pub fn builder<P: AsRef<Path>>(path: P) -> BedBuilder {
        BedBuilder::new(path)
    }

    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedErrorPlus> {
        Bed::builder(path).build()
    }

    // pub fn metadata(&mut self) -> Result<Metadata, BedErrorPlus> {
    //     Ok(Metadata {
    //         iid: skip_to_none(self.iid())?,
    //         sid: skip_to_none(self.sid())?,
    //         // chromosome: &self.chromosome,
    //         // allele_1: &self.allele_1,
    //         // allele_2: &self.allele_2,
    //         // fid: &self.fid,
    //         // father: &self.father,
    //         // mother: &self.mother,
    //     })
    // }

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
            FamOrBim::Fam => to_metadata_path(&self.path, &self.fam_path, "fam"),
            FamOrBim::Bim => to_metadata_path(&self.path, &self.bim_path, "bim"),
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

    pub fn metadata(&mut self) -> Result<&Metadata, BedErrorPlus> {
        self.iid()?;
        self.sid()?;
        Ok(&self.metadata)
    }

    pub fn iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.metadata.iid {
            LazyOrSkip::Some(ref iid) => Ok(iid),
            LazyOrSkip::Lazy => {
                let iid = self.read_fam_or_bim(FamOrBim::Fam, 1)?;
                self.metadata.iid = LazyOrSkip::Some(iid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.metadata.iid.as_ref().unwrap())
            }
            LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata("iid".to_string()).into()),
        }
    }

    pub fn sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.metadata.sid {
            LazyOrSkip::Some(ref sid) => Ok(sid),
            LazyOrSkip::Lazy => {
                let sid = self.read_fam_or_bim(FamOrBim::Bim, 2)?;
                self.metadata.sid = LazyOrSkip::Some(sid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.metadata.sid.as_ref().unwrap())
            }
            LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata("sid".to_string()).into()),
        }
    }

    // pub fn metadata(&mut self) -> Result<&Metadata, BedErrorPlus> {
    //     Ok(&self.metadata)
    // }

    pub fn chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.chromosome {
            LazyOrSkip::Some(ref chromosome) => Ok(chromosome),
            LazyOrSkip::Lazy => {
                let chromosome = self.read_fam_or_bim(FamOrBim::Bim, 0)?;
                self.chromosome = LazyOrSkip::Some(chromosome);
                // This unwrap is safe because we just created 'Some'
                Ok(self.chromosome.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("chromosome".to_string()).into())
            }
        }
    }

    pub fn allele_1(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_1 {
            LazyOrSkip::Some(ref allele_1) => Ok(allele_1),
            LazyOrSkip::Lazy => {
                let allele_1 = self.read_fam_or_bim(FamOrBim::Bim, 4)?;
                self.allele_1 = LazyOrSkip::Some(allele_1);
                // This unwrap is safe because we just created 'Some'
                Ok(self.allele_1.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("allele_1".to_string()).into())
            }
        }
    }

    pub fn allele_2(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.allele_2 {
            LazyOrSkip::Some(ref allele_2) => Ok(allele_2),
            LazyOrSkip::Lazy => {
                let allele_2 = self.read_fam_or_bim(FamOrBim::Bim, 5)?;
                self.allele_2 = LazyOrSkip::Some(allele_2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.allele_2.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("allele_2".to_string()).into())
            }
        }
    }

    pub fn fid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.fid {
            LazyOrSkip::Some(ref fid) => Ok(fid),
            LazyOrSkip::Lazy => {
                let fid = self.read_fam_or_bim(FamOrBim::Fam, 0)?;
                self.fid = LazyOrSkip::Some(fid);
                // This unwrap is safe because we just created 'Some'
                Ok(self.fid.as_ref().unwrap())
            }
            LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata("fid".to_string()).into()),
        }
    }

    pub fn father(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.father {
            LazyOrSkip::Some(ref father) => Ok(father),
            LazyOrSkip::Lazy => {
                let father = self.read_fam_or_bim(FamOrBim::Fam, 2)?;
                self.father = LazyOrSkip::Some(father);
                // This unwrap is safe because we just created 'Some'
                Ok(self.father.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("father".to_string()).into())
            }
        }
    }

    pub fn mother(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.mother {
            LazyOrSkip::Some(ref mother) => Ok(mother),
            LazyOrSkip::Lazy => {
                let mother = self.read_fam_or_bim(FamOrBim::Fam, 3)?;
                self.mother = LazyOrSkip::Some(mother);
                // This unwrap is safe because we just created 'Some'
                Ok(self.mother.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("mother".to_string()).into())
            }
        }
    }

    pub fn sex(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.sex {
            LazyOrSkip::Some(ref sex) => Ok(sex),
            LazyOrSkip::Lazy => {
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
                self.sex = LazyOrSkip::Some(sex2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.sex.as_ref().unwrap())
            }
            LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata("sex".to_string()).into()),
        }
    }

    pub fn cm_position(&mut self) -> Result<&nd::Array1<f32>, BedErrorPlus> {
        match self.cm_position {
            LazyOrSkip::Some(ref cm_position) => Ok(cm_position),
            LazyOrSkip::Lazy => {
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
                self.cm_position = LazyOrSkip::Some(cm_position2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.cm_position.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("cm_position".to_string()).into())
            }
        }
    }

    pub fn bp_position(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        match self.bp_position {
            LazyOrSkip::Some(ref bp_position) => Ok(bp_position),
            LazyOrSkip::Lazy => {
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
                self.bp_position = LazyOrSkip::Some(bp_position2);
                // This unwrap is safe because we just created 'Some'
                Ok(self.bp_position.as_ref().unwrap())
            }
            LazyOrSkip::Skip => {
                Err(BedError::CannotUseSkippedMetadata("bp_position".to_string()).into())
            }
        }
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
    metadata: Metadata,
}

impl WriteOptions {
    pub fn builder<P: AsRef<Path>>(path: P) -> WriteOptionsBuilder {
        WriteOptionsBuilder::new(path)
    }

    pub fn iid(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.metadata.iid {
            LazyOrSkip::Some(ref iid) => Ok(iid),
            LazyOrSkip::Lazy => {
                panic!("!!!cmk 0 raise error here");
            }
            LazyOrSkip::Skip => {
                panic!("!!!cmk 0 maybe return default value");
            }
        }
    }

    pub fn sid(&self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        match self.metadata.sid {
            LazyOrSkip::Some(ref iid) => Ok(iid),
            LazyOrSkip::Lazy => {
                panic!("!!!cmk 0 raise error here");
            }
            LazyOrSkip::Skip => {
                panic!("!!!cmk 0 maybe return default value");
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
            metadata: Some(Metadata::default()),
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

    pub fn metadata(mut self, metadata: Metadata) -> Self {
        // !!!cmk later check that no metadata is still 'Lazy'
        if let LazyOrSkip::Some(iid) = &metadata.iid {
            self = self.iid(iid);
        }
        if let LazyOrSkip::Some(sid) = &metadata.sid {
            self = self.sid(sid);
        }
        self
    }

    pub fn iid<I, T>(mut self, iid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            iid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.metadata.as_mut().unwrap().iid = LazyOrSkip::Some(new_string_array);
        self
    }

    pub fn sid<I, T>(mut self, sid: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<str>,
    {
        let new_string_array: nd::Array1<String> =
            sid.into_iter().map(|s| s.as_ref().to_string()).collect();
        self.metadata.as_mut().unwrap().sid = LazyOrSkip::Some(new_string_array);
        self
    }
}

pub fn write_with_options<TVal: BedVal>(
    val: &nd::Array2<TVal>,
    write_options: &WriteOptions,
) -> Result<(), BedErrorPlus> {
    // !!!cmk later if something goes wrong, clean up the files?
    let path = &write_options.path;

    // !!!cmk 0 use fam/bim
    // !!!cmk 0 set is_a1_count
    // !!!cmk 0 set missing
    // !!!cmk set num_threads
    write(path, &val.view(), true, TVal::missing(), 0)?;

    let _fam_path = to_metadata_path(path, &write_options.fam_path, "fam");
    let _bim_path = to_metadata_path(path, &write_options.bim_path, "bim");

    let _iid = write_options.iid()?;
    let _sid = write_options.sid()?;

    todo!("!!!cmk 0 write metadata");

    // Ok(())
}
