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
use std::ops::{Bound, Range, RangeBounds, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive};
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter},
    ops::RangeFull,
    path::{Path, PathBuf},
};

use derive_builder::Builder;
// !!! might want to use this instead use typed_builder::TypedBuilder;

use crate::{
    count_lines, open_and_check, read_no_alloc, write_val, BedError, BedErrorPlus, BedVal,
};

#[derive(Debug, Clone)]
pub enum LazyOrSkip<T> {
    Lazy,
    Skip,
    Some(T),
}

impl<T> LazyOrSkip<T> {
    pub fn is_lazy(&self) -> bool {
        match self {
            LazyOrSkip::Lazy => true,
            _ => false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Skippable<T> {
    Some(T),
    Skip,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Metadata<'a> {
    pub fid: Skippable<&'a nd::Array1<String>>,
    pub iid: Skippable<&'a nd::Array1<String>>,
    pub father: Skippable<&'a nd::Array1<String>>,
    pub mother: Skippable<&'a nd::Array1<String>>,
    pub sex: Skippable<&'a nd::Array1<i32>>,
    pub pheno: Skippable<&'a nd::Array1<String>>,

    pub chromosome: Skippable<&'a nd::Array1<String>>,
    pub sid: Skippable<&'a nd::Array1<String>>,
    pub cm_position: Skippable<&'a nd::Array1<f32>>,
    pub bp_position: Skippable<&'a nd::Array1<i32>>,
    pub allele_1: Skippable<&'a nd::Array1<String>>,
    pub allele_2: Skippable<&'a nd::Array1<String>>,
}

fn option_count<T>(array: &LazyOrSkip<nd::Array1<T>>) -> Option<usize> {
    match array {
        LazyOrSkip::Some(array) => Some(array.len()),
        LazyOrSkip::Skip => None,
        LazyOrSkip::Lazy => None,
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

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,

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

    // i32 based on https://www.cog-genomics.org/plink2/formats#bim
    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    bp_position: LazyOrSkip<nd::Array1<i32>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_1: LazyOrSkip<nd::Array1<String>>,

    #[builder(setter(custom))]
    #[builder(default = "LazyOrSkip::Lazy")]
    allele_2: LazyOrSkip<nd::Array1<String>>,
}

impl BedBuilder {
    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        let path: PathBuf = path.as_ref().into();

        Self {
            path: Some(path),
            fam_path: None,
            bim_path: None,

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
        }
    }

    fn check_counts(
        count_vec: Vec<Option<usize>>,
        option_xid_count: &mut Option<usize>,
        prefix: &str,
    ) -> Result<(), BedErrorPlus> {
        for option_count in count_vec {
            if let Some(count) = option_count {
                match option_xid_count {
                    Some(xid_count) => {
                        if *xid_count != count {
                            return Err(BedError::InconsistentCount(
                                prefix.to_string(),
                                *xid_count,
                                count,
                            )
                            .into());
                        }
                    }
                    None => {
                        *option_xid_count = Some(count);
                    }
                }
            }
        }
        Ok(())
    }

    pub fn build(&self) -> Result<Bed, BedErrorPlus> {
        let mut bed = self.build_no_file_check()?;

        if bed.is_checked_early {
            open_and_check(&bed.path)?;
        }

        BedBuilder::check_counts(
            vec![
                option_count(&bed.fid),
                option_count(&bed.iid),
                option_count(&bed.father),
                option_count(&bed.mother),
                option_count(&bed.sex),
                option_count(&bed.pheno),
            ],
            &mut bed.iid_count,
            "iid",
        )?;

        BedBuilder::check_counts(
            vec![
                option_count(&bed.chromosome),
                option_count(&bed.sid),
                option_count(&bed.cm_position),
                option_count(&bed.bp_position),
                option_count(&bed.allele_1),
                option_count(&bed.allele_2),
            ],
            &mut bed.sid_count,
            "sid",
        )?;

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

    pub fn skip_early_check(mut self) -> Self {
        self.is_checked_early = Some(false);
        self
    }

    pub fn skip_iid(mut self) -> Self {
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

    pub fn sex<I: IntoIterator<Item = i32>>(mut self, sex: I) -> Self {
        self.sex = Some(LazyOrSkip::Some(sex.into_iter().map(|s| s).collect()));
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

    pub fn cm_position<I: IntoIterator<Item = f32>>(mut self, cm_position: I) -> Self {
        self.cm_position = Some(LazyOrSkip::Some(
            cm_position.into_iter().map(|s| s).collect(),
        ));
        self
    }

    pub fn bp_position<I: IntoIterator<Item = i32>>(mut self, bp_position: I) -> Self {
        self.bp_position = Some(LazyOrSkip::Some(
            bp_position.into_iter().map(|s| s).collect(),
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

fn to_skippable<'a, T>(lazy_or_skip: &'a LazyOrSkip<T>) -> Skippable<&'a T> {
    match lazy_or_skip {
        LazyOrSkip::Lazy => panic!("assert: impossible"),
        LazyOrSkip::Skip => Skippable::Skip,
        LazyOrSkip::Some(some) => Skippable::Some(some),
    }
}

impl Bed {
    pub fn builder<P: AsRef<Path>>(path: P) -> BedBuilder {
        BedBuilder::new(path)
    }

    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedErrorPlus> {
        Bed::builder(path).build()
    }

    pub fn fam_path(&mut self) -> Result<PathBuf, BedErrorPlus> {
        if let Some(path) = &self.fam_path {
            Ok(path.clone())
        } else {
            let path = to_metadata_path(&self.path, &self.fam_path, "fam");
            self.fam_path = Some(path.clone());
            Ok(path)
        }
    }

    pub fn bim_path(&mut self) -> Result<PathBuf, BedErrorPlus> {
        if let Some(path) = &self.bim_path {
            Ok(path.clone())
        } else {
            let path = to_metadata_path(&self.path, &self.bim_path, "bim");
            self.bim_path = Some(path.clone());
            Ok(path)
        }
    }

    pub fn iid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let fam_path = self.fam_path()?;
            let iid_count = count_lines(fam_path)?;
            self.iid_count = Some(iid_count);
            Ok(iid_count)
        }
    }
    pub fn sid_count(&mut self) -> Result<usize, BedErrorPlus> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let bim_path = self.bim_path()?;
            let sid_count = count_lines(bim_path)?;
            self.sid_count = Some(sid_count);
            Ok(sid_count)
        }
    }

    pub fn metadata(&mut self) -> Result<Metadata, BedErrorPlus> {
        self.fam()?;
        self.bim()?;

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

    fn fam(&mut self) -> Result<(), BedErrorPlus> {
        let mut field_vec: Vec<usize> = Vec::new();
        if self.fid.is_lazy() {
            field_vec.push(0);
        }
        if self.iid.is_lazy() {
            field_vec.push(1);
        }
        if self.father.is_lazy() {
            field_vec.push(2);
        }
        if self.mother.is_lazy() {
            field_vec.push(3);
        }
        if self.sex.is_lazy() {
            field_vec.push(4);
        }
        if self.pheno.is_lazy() {
            field_vec.push(5);
        }

        let fam_path = self.fam_path()?;
        let (mut vec_of_vec, count) = self.read_fam_or_bim(&field_vec, &fam_path)?;
        match self.iid_count {
            Some(iid_count) => {
                if iid_count != count {
                    return Err(
                        BedError::InconsistentCount("iid".to_string(), iid_count, count).into(),
                    );
                }
            }
            None => {
                self.iid_count = Some(count);
            }
        }

        // unwraps are safe because we pop once for every push
        if self.pheno.is_lazy() {
            self.pheno = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.sex.is_lazy() {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?; // !!!cmk later test this error
            self.sex = LazyOrSkip::Some(array);
        }
        if self.mother.is_lazy() {
            self.mother = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.father.is_lazy() {
            self.father = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.iid.is_lazy() {
            self.iid = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.fid.is_lazy() {
            self.fid = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }

        Ok(())
    }
    fn bim(&mut self) -> Result<(), BedErrorPlus> {
        let mut field_vec: Vec<usize> = Vec::new();
        if self.chromosome.is_lazy() {
            field_vec.push(0);
        }
        if self.sid.is_lazy() {
            field_vec.push(1);
        }
        if self.cm_position.is_lazy() {
            field_vec.push(2);
        }
        if self.bp_position.is_lazy() {
            field_vec.push(3);
        }
        if self.allele_1.is_lazy() {
            field_vec.push(4);
        }
        if self.allele_2.is_lazy() {
            field_vec.push(5);
        }

        let bim_path = self.bim_path()?;
        let (mut vec_of_vec, count) = self.read_fam_or_bim(&field_vec, &bim_path)?;
        match self.sid_count {
            Some(sid_count) => {
                if sid_count != count {
                    return Err(
                        BedError::InconsistentCount("sid".to_string(), sid_count, count).into(),
                    );
                }
            }
            None => {
                self.sid_count = Some(count);
            }
        }

        // unwraps are safe because we pop once for every push
        if self.allele_2.is_lazy() {
            self.allele_2 = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.allele_1.is_lazy() {
            self.allele_1 = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.bp_position.is_lazy() {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?; // !!!cmk later test this error
            self.bp_position = LazyOrSkip::Some(array);
        }
        if self.cm_position.is_lazy() {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<f32>())
                .collect::<Result<nd::Array1<f32>, _>>()?; // !!!cmk later test this error
            self.cm_position = LazyOrSkip::Some(array);
        }

        if self.sid.is_lazy() {
            self.sid = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }
        if self.chromosome.is_lazy() {
            self.chromosome = LazyOrSkip::Some(nd::Array::from_vec(vec_of_vec.pop().unwrap()));
        }

        Ok(())
    }

    fn read_fam_or_bim(
        &self,
        field_vec: &Vec<usize>,
        path_buf: &PathBuf,
    ) -> Result<(Vec<Vec<String>>, usize), BedErrorPlus> {
        let mut vec_of_vec = vec![vec![]; field_vec.len()];

        let file = File::open(&path_buf)?;

        let reader = BufReader::new(file);
        let mut count = 0;
        for line in reader.lines() {
            let line = line?;
            count += 1;
            let field = line.split_whitespace();

            let mut field_count = 0;
            let mut of_interest_count = 0;
            for field in field {
                if field_vec.contains(&field_count) {
                    vec_of_vec[of_interest_count].push(field.to_string());
                    of_interest_count += 1;
                }
                field_count += 1;
            }
            if field_count != 6 {
                return Err(BedError::MetadataFieldCount(
                    6,
                    field_count,
                    path_buf.to_str().unwrap().to_string(),
                )
                .into());
            }
        }

        Ok((vec_of_vec, count))
    }

    fn unlazy_fam<T: FromStringArray<T>>(&mut self, is_lazy: bool) -> Result<(), BedErrorPlus> {
        if is_lazy {
            self.fam()?
        }
        Ok(())
    }
    fn unlazy_bim<T: FromStringArray<T>>(&mut self, is_lazy: bool) -> Result<(), BedErrorPlus> {
        if is_lazy {
            self.bim()?
        }
        Ok(())
    }

    fn get_some_field<'a, T: FromStringArray<T>>(
        &'a self,
        field: &'a LazyOrSkip<nd::Array1<T>>,
        name: &str,
    ) -> Result<&'a nd::Array1<T>, BedErrorPlus> {
        match field {
            LazyOrSkip::Some(array) => Ok(array),
            LazyOrSkip::Skip => Err(BedError::CannotUseSkippedMetadata(name.to_string()).into()),
            LazyOrSkip::Lazy => panic!("impossible"),
        }
    }

    pub fn fid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.fid.is_lazy())?;
        self.get_some_field(&self.fid, "fid")
    }

    pub fn iid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.iid.is_lazy())?;
        self.get_some_field(&self.iid, "iid")
    }
    pub fn father(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.father.is_lazy())?;
        self.get_some_field(&self.father, "father")
    }
    pub fn mother(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.mother.is_lazy())?;
        self.get_some_field(&self.mother, "mother")
    }
    pub fn sex(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.sex.is_lazy())?;
        self.get_some_field(&self.sex, "sex")
    }
    pub fn pheno(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_fam::<String>(self.pheno.is_lazy())?;
        self.get_some_field(&self.pheno, "pheno")
    }

    pub fn chromosome(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.chromosome.is_lazy())?;
        self.get_some_field(&self.chromosome, "chromosome")
    }
    pub fn sid(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.sid.is_lazy())?;
        self.get_some_field(&self.sid, "sid")
    }
    pub fn cm_position(&mut self) -> Result<&nd::Array1<f32>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.cm_position.is_lazy())?;
        self.get_some_field(&self.cm_position, "cm_position")
    }
    pub fn bp_position(&mut self) -> Result<&nd::Array1<i32>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.bp_position.is_lazy())?;
        self.get_some_field(&self.bp_position, "bp_position")
    }
    pub fn allele_1(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.allele_1.is_lazy())?;
        self.get_some_field(&self.allele_1, "allele_1")
    }
    pub fn allele_2(&mut self) -> Result<&nd::Array1<String>, BedErrorPlus> {
        self.unlazy_bim::<String>(self.allele_2.is_lazy())?;
        self.get_some_field(&self.allele_2, "allele_2")
    }

    pub fn read<TVal: BedVal>(&mut self) -> Result<nd::Array2<TVal>, BedErrorPlus> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_with_options(read_options)
    }

    pub fn read_and_fill_with_options<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
        read_options: ReadOptions<TVal>,
    ) -> Result<(), BedErrorPlus> {
        let iid_count = self.iid_count()?;
        let sid_count = self.sid_count()?;

        let num_threads = compute_num_threads(read_options.num_threads)?;

        let iid_index = read_options.iid_index.to_vec(iid_count)?;
        let sid_index = read_options.sid_index.to_vec(sid_count)?;

        // !!!cmk 0 check that val has right (iid_index.len(), sid_index.len()), read_options.is_f);

        read_no_alloc(
            &self.path,
            iid_count,
            sid_count,
            read_options.is_a1_counted,
            &iid_index,
            &sid_index,
            read_options.missing_value,
            num_threads,
            &mut val.view_mut(),
        )?;

        Ok(())
    }

    pub fn read_and_fill<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
    ) -> Result<(), BedErrorPlus> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        let num_threads = compute_num_threads(read_options.num_threads)?;

        let iid_count = self.iid_count()?;
        let sid_count = self.sid_count()?;
        let iid_index = read_options.iid_index.to_vec(iid_count)?;
        let sid_index = read_options.sid_index.to_vec(sid_count)?;

        read_no_alloc(
            &self.path,
            iid_count,
            sid_count,
            read_options.is_a1_counted,
            &iid_index,
            &sid_index,
            read_options.missing_value,
            num_threads,
            &mut val.view_mut(),
        )?;

        Ok(())
    }

    pub fn read_with_options<TVal: BedVal>(
        &mut self,
        read_options: ReadOptions<TVal>,
    ) -> Result<nd::Array2<TVal>, BedErrorPlus> {
        let iid_count_in = self.iid_count()?;
        let sid_count_in = self.sid_count()?;
        let iid_count_out = read_options.iid_index.len(iid_count_in)?;
        let sid_count_out = read_options.sid_index.len(sid_count_in)?;
        let shape = ShapeBuilder::set_f((iid_count_out, sid_count_out), read_options.is_f);
        let mut val = nd::Array2::<TVal>::default(shape);

        self.read_and_fill_with_options(&mut val.view_mut(), read_options)?;

        Ok(val)
    }
}

fn compute_num_threads(option_num_threads: Option<usize>) -> Result<usize, BedErrorPlus> {
    let num_threads = if let Some(num_threads) = option_num_threads {
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
    Ok(num_threads)
}

impl Index {
    // !!!cmk later test every case
    // We can't define a 'From' because we want to add count at the last moment.
    // Would be nice to not always allocate a new vec, maybe with Rc<[T]>?
    // Even better would be to support an iterator from Index (an enum with fields).
    pub fn to_vec(&self, count: usize) -> Result<Vec<usize>, BedErrorPlus> {
        match self {
            Index::All => Ok((0..count).collect()),
            Index::Vec(vec) => Ok(vec.to_vec()),
            Index::NDArrayBool(nd_array_bool) => {
                if nd_array_bool.len() != count {
                    return Err(
                        BedError::BoolArrayVectorWrongLength(count, nd_array_bool.len()).into(),
                    );
                }
                Ok(nd_array_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i)
                    .collect())
            }
            // !!!cmk later can we implement this without two allocations?
            Index::NDSliceInfo(nd_slice_info) => {
                Ok(RangeNdSlice::new(nd_slice_info, count)?.to_vec())
            }
            Index::RangeAny(range_any) => {
                let range = range_any.to_range(count)?;
                Ok(range.collect::<Vec<usize>>())
            }
            Index::NDArray(nd_array) => Ok(nd_array.to_vec()),
            Index::One(one) => Ok(vec![*one]),
            Index::VecBool(vec_bool) => {
                if vec_bool.len() != count {
                    return Err(BedError::BoolArrayVectorWrongLength(count, vec_bool.len()).into());
                }
                Ok(vec_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i)
                    .collect())
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
    All,
    One(usize),
    Vec(Vec<usize>),
    NDArray(nd::Array1<usize>),
    VecBool(Vec<bool>),
    NDArrayBool(nd::Array1<bool>),
    NDSliceInfo(SliceInfo1),
    RangeAny(RangeAny),
}

#[derive(Debug, Clone)]
pub struct RangeAny {
    start: Option<usize>,
    end: Option<usize>,
}

impl RangeAny {
    // https://stackoverflow.com/questions/55925523/array-cannot-be-indexed-by-rangefull
    fn to_range(&self, count: usize) -> Result<Range<usize>, BedErrorPlus> {
        let start = if let Some(start) = self.start {
            start
        } else {
            0
        };
        let end = if let Some(end) = self.end { end } else { count };
        if start > end {
            Err(BedError::StartGreaterThanEnd(start, end).into())
        } else {
            Ok(Range {
                start: start,
                end: end,
            })
        }
    }

    fn len(&self, count: usize) -> Result<usize, BedErrorPlus> {
        let range = self.to_range(count)?;
        Ok(range.end - range.start)
    }
}

#[derive(Debug, Clone)]
pub struct RangeNdSlice {
    start: usize,
    end: usize,
    step: usize,
    is_reversed: bool,
}

// https://www.geeksforgeeks.org/find-ceil-ab-without-using-ceil-function/
fn div_ceil(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

impl RangeNdSlice {
    fn len(&self) -> usize {
        if self.start > self.end {
            0
        } else {
            div_ceil(self.end - self.start, self.step)
        }
    }

    // https://docs.rs/ndarray/0.15.4/ndarray/struct.ArrayBase.html#slicing
    fn to_vec(&self) -> Vec<usize> {
        if self.start > self.end {
            Vec::new()
        } else {
            let mut vec: Vec<usize> = (self.start..self.end).step_by(self.step).collect();
            if self.is_reversed {
                vec.reverse();
            }
            vec
        }
    }

    fn new(nd_slice_info: &SliceInfo1, count: usize) -> Result<Self, BedErrorPlus> {
        //  self.to_vec(count).len(),
        // https://docs.rs/ndarray/0.15.4/ndarray/struct.ArrayBase.html#method.slice_collapse
        // Error in the following cases
        // * SliceInfo is not a 1-dimensional or is a NewAxis
        // * Step is 0
        // * Start is greater than count
        // * End is greater than count
        // As with ndarray, Start can be greater than End is allowed
        // and means the slice is empty.
        if nd_slice_info.in_ndim() != 1 || nd_slice_info.out_ndim() != 1 {
            return Err(BedError::NdSliceInfoNot1D.into());
        }

        let slice_info_elem = nd_slice_info[0];
        match slice_info_elem {
            nd::SliceInfoElem::Slice { start, end, step } => {
                // https://docs.rs/ndarray/0.15.4/ndarray/enum.SliceInfoElem.html
                // s![..], 0,None,1
                // s![a..b;2] a,b,2
                // s![a..;-1], from a to end in reverse order
                // start index; negative are counted from the back of the axis
                // end index; negative are counted from the back of the axis; when not present the default is the full length of the axis.
                // step size in elements; the default is 1, for every element.
                // A range with step size. end is an exclusive index. Negative start or end indexes are counted from the back of the axis. If end is None, the slice extends to the end of the axis.
                let step2: usize;
                let is_reverse2: bool;
                if step > 0 {
                    step2 = step as usize;
                    is_reverse2 = false;
                } else if step < 0 {
                    step2 = (-step) as usize;
                    is_reverse2 = true;
                } else {
                    return Err(BedError::StepZero.into());
                }

                let start2 = if start >= 0 {
                    let start3 = start as usize;
                    if start3 > count {
                        return Err(BedError::StartGreaterThanCount(start3, count).into());
                    } else {
                        start3
                    }
                } else {
                    let start3 = (-start) as usize;
                    if start3 > count {
                        return Err(BedError::StartGreaterThanCount(start3, count).into());
                    }
                    count - start3
                };

                let end2 = if let Some(end) = end {
                    if end >= 0 {
                        let end3 = end as usize;
                        if end3 > count {
                            return Err(BedError::EndGreaterThanCount(end3, count).into());
                        } else {
                            end3
                        }
                    } else {
                        let end3 = (-end) as usize;
                        if end3 > count {
                            return Err(BedError::EndGreaterThanCount(end3, count).into());
                        }
                        count - end3
                    }
                } else {
                    count
                };

                Ok(RangeNdSlice {
                    start: start2,
                    end: end2,
                    step: step2,
                    is_reversed: is_reverse2,
                })
            }
            nd::SliceInfoElem::Index(index) => Ok(RangeNdSlice {
                start: index as usize,
                end: index as usize + 1,
                step: 1,
                is_reversed: false,
            }),
            nd::SliceInfoElem::NewAxis => {
                return Err(BedError::NewAxis.into());
            }
        }
    }
}

impl Index {
    pub fn len(&self, count: usize) -> Result<usize, BedErrorPlus> {
        match self {
            Index::All => Ok(count),
            Index::One(_) => Ok(1),
            Index::Vec(vec) => Ok(vec.len()),
            Index::NDArray(nd_array) => Ok(nd_array.len()),
            Index::VecBool(vec_bool) => Ok(vec_bool.iter().filter(|&b| *b).count()),
            Index::NDArrayBool(nd_array_bool) => Ok(nd_array_bool.iter().filter(|&b| *b).count()),
            Index::NDSliceInfo(nd_slice_info) => Ok(RangeNdSlice::new(nd_slice_info, count)?.len()),
            Index::RangeAny(range_any) => range_any.len(count),
        }
    }
}

// !!!cmk later see if what ref conversions. See https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust
impl From<SliceInfo1> for Index {
    fn from(slice_info: SliceInfo1) -> Index {
        Index::NDSliceInfo(slice_info)
    }
}

fn to_range_any<T: RangeBounds<usize>>(range_thing: T) -> RangeAny {
    let start_bound = range_thing.start_bound();
    let start = match start_bound {
        Bound::Included(&start) => Some(start),
        Bound::Excluded(&start) => Some(start + 1),
        Bound::Unbounded => None,
    };

    let end_bound = range_thing.end_bound();
    let end = match end_bound {
        Bound::Included(&end) => Some(end + 1),
        Bound::Excluded(&end) => Some(end),
        Bound::Unbounded => None,
    };
    RangeAny { start, end }
}

impl From<RangeFull> for RangeAny {
    fn from(range_thing: RangeFull) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<RangeFull> for Index {
    fn from(range_thing: RangeFull) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<Range<usize>> for RangeAny {
    fn from(range_thing: Range<usize>) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<Range<usize>> for Index {
    fn from(range_thing: Range<usize>) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<RangeFrom<usize>> for RangeAny {
    fn from(range_thing: RangeFrom<usize>) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<RangeFrom<usize>> for Index {
    fn from(range_thing: RangeFrom<usize>) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<RangeInclusive<usize>> for RangeAny {
    fn from(range_thing: RangeInclusive<usize>) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<RangeInclusive<usize>> for Index {
    fn from(range_thing: RangeInclusive<usize>) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<RangeTo<usize>> for RangeAny {
    fn from(range_thing: RangeTo<usize>) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<RangeTo<usize>> for Index {
    fn from(range_thing: RangeTo<usize>) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<RangeToInclusive<usize>> for RangeAny {
    fn from(range_thing: RangeToInclusive<usize>) -> RangeAny {
        to_range_any(range_thing)
    }
}

impl From<RangeToInclusive<usize>> for Index {
    fn from(range_thing: RangeToInclusive<usize>) -> Index {
        Index::RangeAny(range_thing.into())
    }
}

impl From<&[usize]> for Index {
    fn from(array: &[usize]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl<const N: usize> From<[usize; N]> for Index {
    fn from(array: [usize; N]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl<const N: usize> From<&[usize; N]> for Index {
    fn from(array: &[usize; N]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl From<&nd::ArrayView1<'_, usize>> for Index {
    fn from(view: &nd::ArrayView1<usize>) -> Index {
        Index::NDArray(view.to_owned())
    }
}

impl From<&Vec<usize>> for Index {
    fn from(vec_ref: &Vec<usize>) -> Index {
        Index::Vec(vec_ref.clone())
    }
}

// !!!cmk later is ref &ndarray const array and bool OK
impl From<&nd::ArrayView1<'_, bool>> for Index {
    fn from(view: &nd::ArrayView1<bool>) -> Index {
        Index::NDArrayBool(view.to_owned())
    }
}

impl From<&Vec<bool>> for Index {
    fn from(vec_ref: &Vec<bool>) -> Index {
        Index::VecBool(vec_ref.clone())
    }
}

impl From<&[bool]> for Index {
    fn from(array: &[bool]) -> Index {
        Index::VecBool(array.to_vec())
    }
}

impl<const N: usize> From<[bool; N]> for Index {
    fn from(array: [bool; N]) -> Index {
        Index::VecBool(array.to_vec())
    }
}

impl<const N: usize> From<&[bool; N]> for Index {
    fn from(array: &[bool; N]) -> Index {
        Index::VecBool(array.to_vec())
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
        Index::All
    }
}

// See https://nullderef.com/blog/rust-parameters/
#[derive(Debug, Clone, Builder)]
pub struct ReadOptions<TVal: BedVal> {
    #[builder(default = "TVal::missing()")]
    missing_value: TVal,

    #[builder(default = "Index::All")]
    #[builder(setter(into))]
    iid_index: Index,

    #[builder(default = "Index::All")]
    #[builder(setter(into))]
    sid_index: Index,

    #[builder(default = "true")]
    is_f: bool,

    #[builder(default = "true")]
    is_a1_counted: bool,

    #[builder(default, setter(strip_option))]
    pub num_threads: Option<usize>,
}

impl<TVal: BedVal> ReadOptions<TVal> {
    pub fn builder() -> ReadOptionsBuilder<TVal> {
        ReadOptionsBuilder::default()
    }
}

impl<TVal: BedVal> ReadOptionsBuilder<TVal> {
    pub fn read(&self, bed: &mut Bed) -> Result<nd::Array2<TVal>, BedErrorPlus> {
        let read_option = self.build()?;
        bed.read_with_options(read_option)
    }
    pub fn read_and_fill(
        &self,
        bed: &mut Bed,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
    ) -> Result<(), BedErrorPlus> {
        let read_option = self.build()?;
        bed.read_and_fill_with_options(val, read_option)
    }

    pub fn f(&mut self) -> &Self {
        self.is_f(true);
        self
    }

    pub fn c(&mut self) -> &Self {
        self.is_f(false);
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
pub struct WriteOptions<TVal>
where
    TVal: BedVal,
{
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

    #[builder(default = "true")]
    is_a1_counted: bool,

    #[builder(default, setter(strip_option))]
    pub num_threads: Option<usize>,

    #[builder(default = "TVal::missing()")]
    missing_value: TVal,
}

impl<TVal> WriteOptions<TVal>
where
    TVal: BedVal,
{
    pub fn builder<P: AsRef<Path>>(path: P) -> WriteOptionsBuilder<TVal> {
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

    pub fn is_a1_counted(&self) -> bool {
        self.is_a1_counted
    }
}

impl<TVal> WriteOptionsBuilder<TVal>
where
    TVal: BedVal,
{
    // !!! cmk later just use the default builder?
    pub fn build(&self) -> Result<WriteOptions<TVal>, BedErrorPlus> {
        let write_options = self.write_no_file_check()?;

        Ok(write_options)
    }

    // !!!cmk later should check that metadata agrees with val size
    // !!!cmk later maybe use the default builder?
    pub fn write<S: nd::Data<Elem = TVal>>(
        &self,
        val: &nd::ArrayBase<S, nd::Ix2>,
    ) -> Result<(), BedErrorPlus> {
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

            is_a1_counted: None,
            num_threads: None,
            missing_value: None,
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

    // !!!cmk later can we also extract a metadata property from write options?
    pub fn metadata(mut self, metadata: &Metadata) -> Self {
        if let Skippable::Some(fid) = &metadata.fid {
            self.fid = Some(Skippable::Some((*fid).clone()));
        }
        if let Skippable::Some(iid) = &metadata.iid {
            self.iid = Some(Skippable::Some((*iid).clone()));
        }
        if let Skippable::Some(father) = &metadata.father {
            self.father = Some(Skippable::Some((*father).clone()));
        }
        if let Skippable::Some(mother) = &metadata.mother {
            self.mother = Some(Skippable::Some((*mother).clone()));
        }
        if let Skippable::Some(sex) = &metadata.sex {
            self.sex = Some(Skippable::Some((*sex).clone()));
        }
        if let Skippable::Some(pheno) = &metadata.pheno {
            self.pheno = Some(Skippable::Some((*pheno).clone()));
        }

        if let Skippable::Some(chromosome) = &metadata.chromosome {
            self.chromosome = Some(Skippable::Some((*chromosome).clone()));
        }
        if let Skippable::Some(sid) = &metadata.sid {
            self.sid = Some(Skippable::Some((*sid).clone()));
        }
        if let Skippable::Some(cm_position) = &metadata.cm_position {
            self.cm_position = Some(Skippable::Some((*cm_position).clone()));
        }
        if let Skippable::Some(bp_position) = &metadata.bp_position {
            self.bp_position = Some(Skippable::Some((*bp_position).clone()));
        }
        if let Skippable::Some(allele_1) = &metadata.allele_1 {
            self.allele_1 = Some(Skippable::Some((*allele_1).clone()));
        }
        if let Skippable::Some(allele_2) = &metadata.allele_2 {
            self.allele_2 = Some(Skippable::Some((*allele_2).clone()));
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

    pub fn sex<I: IntoIterator<Item = i32>>(mut self, sex: I) -> Self {
        self.sex = Some(Skippable::Some(sex.into_iter().map(|i| i).collect()));
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

    pub fn cm_position<I: IntoIterator<Item = f32>>(mut self, cm_position: I) -> Self {
        self.cm_position = Some(Skippable::Some(
            cm_position.into_iter().map(|s| s).collect(),
        ));
        self
    }

    pub fn bp_position<I: IntoIterator<Item = i32>>(mut self, bp_position: I) -> Self {
        self.bp_position = Some(Skippable::Some(
            bp_position.into_iter().map(|s| s).collect(),
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

    pub fn count_a1(&mut self) -> &mut Self {
        self.is_a1_counted = Some(true);
        self
    }

    pub fn count_a2(&mut self) -> &mut Self {
        self.is_a1_counted = Some(false);
        self
    }
}

pub fn write<S: nd::Data<Elem = TVal>, TVal: BedVal>(
    val: &nd::ArrayBase<S, nd::Ix2>,
    path: &Path,
) -> Result<(), BedErrorPlus> {
    WriteOptions::builder(path).write(val)
}

// !!!cmk 0 rename
// !!!cmk later do this without a "clone"
fn compute_field<T, F>(field: &Skippable<nd::Array1<T>>, count: usize, lambda: F) -> nd::Array1<T>
where
    T: Clone + Default + Debug,
    F: Fn(usize) -> T,
{
    match field {
        Skippable::Some(array) => array.clone(),
        Skippable::Skip => (0..count).map(|_| lambda(0)).collect::<nd::Array1<T>>(),
    }
}

pub fn write_with_options<S, TVal>(
    val: &nd::ArrayBase<S, nd::Ix2>,
    write_options: &WriteOptions<TVal>,
) -> Result<(), BedErrorPlus>
where
    S: nd::Data<Elem = TVal>,
    TVal: BedVal,
{
    // !!!cmk later can this be done in one step??
    let shape = val.shape();
    let iid_count = shape[0];
    let sid_count = shape[1];
    // !!!cmk 0 if something goes wrong, clean up the files?
    let path = &write_options.path;

    let num_threads = compute_num_threads(write_options.num_threads)?;
    write_val(
        path,
        val,
        write_options.is_a1_counted,
        write_options.missing_value,
        num_threads,
    )?;

    let fam_path = to_metadata_path(path, &write_options.fam_path, "fam");
    let bim_path = to_metadata_path(path, &write_options.bim_path, "bim");

    let fid = compute_field(&write_options.fid, iid_count, |_| "0".to_string());
    let iid = compute_field(&write_options.iid, iid_count, |i| format!("iid{}", i + 1));
    let father = compute_field(&write_options.father, iid_count, |_| "0".to_string());
    let mother = compute_field(&write_options.mother, iid_count, |_| "0".to_string());
    let sex = compute_field(&write_options.sex, iid_count, |_| 0);
    let pheno = compute_field(&write_options.pheno, iid_count, |_| "0".to_string());

    let chromosome = compute_field(&write_options.chromosome, sid_count, |_| "0".to_string());
    let sid = compute_field(&write_options.sid, sid_count, |i| format!("sid{}", i + 1));
    let cm_position = compute_field(&write_options.cm_position, sid_count, |_| 0.0);
    let bp_position = compute_field(&write_options.bp_position, sid_count, |_| 0);
    let allele_1 = compute_field(&write_options.allele_1, sid_count, |_| "A1".to_string());
    let allele_2 = compute_field(&write_options.allele_2, sid_count, |_| "A2".to_string());

    let file = File::create(fam_path)?;
    let mut writer = BufWriter::new(file);
    let mut result: Result<(), BedErrorPlus> = Ok(());
    nd::azip!((fid in &fid, iid in &iid, father in &father, mother in &mother, sex in &sex, pheno in &pheno)
    {
        if result.is_ok() {
            if let Err(e) = writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            *fid, *iid, *father, *mother, *sex, *pheno
        )
        {
         result = Err(BedErrorPlus::IOError(e)); // !!!cmk later test this
        }
     }});
    result?;

    let file = File::create(bim_path)?;
    let mut writer = BufWriter::new(file);
    let mut result: Result<(), BedErrorPlus> = Ok(());
    nd::azip!((chromosome in &chromosome, sid in &sid, cm_position in &cm_position, bp_position in &bp_position, allele_1 in &allele_1, allele_2 in &allele_2)
    {
        if result.is_ok() {
            if let Err(e) = writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            *chromosome, *sid, *cm_position, *bp_position, *allele_1, *allele_2
        )
        {
         result = Err(BedErrorPlus::IOError(e)); // !!!cmk later test this
        }
     }});
    result?;

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

// !!!cmk later test these
impl FromStringArray<f32> for f32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<f32>, BedErrorPlus> {
        let result = string_array
            .iter()
            .map(|s| s.parse::<f32>())
            .collect::<Result<nd::Array1<f32>, _>>();
        match result {
            Ok(array) => Ok(array),
            Err(e) => Err(BedErrorPlus::ParseFloatError(e)),
        }
    }
}
impl FromStringArray<i32> for i32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<i32>, BedErrorPlus> {
        let result = string_array
            .iter()
            .map(|s| s.parse::<i32>())
            .collect::<Result<nd::Array1<i32>, _>>();
        match result {
            Ok(array) => Ok(array),
            Err(e) => Err(BedErrorPlus::ParseIntError(e)),
        }
    }
}
