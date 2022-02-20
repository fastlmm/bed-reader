pub mod api {
    use core::fmt::Debug;
    use nd::ShapeBuilder;
    use ndarray as nd;
    use std::{
        fs::File,
        io::{BufRead, BufReader, Read},
        path::Path,
    };
    use derive_builder::Builder;
    use crate::{
        counts, read_no_alloc, BedError, BedErrorPlus, Missing, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
        CB_HEADER_USIZE,
    };
    #[builder(build_fn(skip))]
    pub struct Bed {
        pub filename: String,
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
    #[allow(clippy::all)]
    ///Builder for [`Bed`](struct.Bed.html).
    pub struct BedBuilder {
        filename: ::derive_builder::export::core::option::Option<String>,
        count_a1: ::derive_builder::export::core::option::Option<bool>,
        iid_count: ::derive_builder::export::core::option::Option<Option<usize>>,
        sid_count: ::derive_builder::export::core::option::Option<Option<usize>>,
        iid: ::derive_builder::export::core::option::Option<Option<nd::Array1<String>>>,
        sid: ::derive_builder::export::core::option::Option<Option<nd::Array1<String>>>,
        chromosome: ::derive_builder::export::core::option::Option<Option<nd::Array1<String>>>,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    #[allow(clippy::all)]
    impl ::core::clone::Clone for BedBuilder {
        #[inline]
        fn clone(&self) -> BedBuilder {
            match *self {
                BedBuilder {
                    filename: ref __self_0_0,
                    count_a1: ref __self_0_1,
                    iid_count: ref __self_0_2,
                    sid_count: ref __self_0_3,
                    iid: ref __self_0_4,
                    sid: ref __self_0_5,
                    chromosome: ref __self_0_6,
                } => BedBuilder {
                    filename: ::core::clone::Clone::clone(&(*__self_0_0)),
                    count_a1: ::core::clone::Clone::clone(&(*__self_0_1)),
                    iid_count: ::core::clone::Clone::clone(&(*__self_0_2)),
                    sid_count: ::core::clone::Clone::clone(&(*__self_0_3)),
                    iid: ::core::clone::Clone::clone(&(*__self_0_4)),
                    sid: ::core::clone::Clone::clone(&(*__self_0_5)),
                    chromosome: ::core::clone::Clone::clone(&(*__self_0_6)),
                },
            }
        }
    }
    #[allow(clippy::all)]
    #[allow(dead_code)]
    impl BedBuilder {
        #[allow(unused_mut)]
        pub fn filename(&mut self, value: String) -> &mut Self {
            let mut new = self;
            new.filename = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn count_a1(&mut self, value: bool) -> &mut Self {
            let mut new = self;
            new.count_a1 = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn iid_count(&mut self, value: usize) -> &mut Self {
            let mut new = self;
            new.iid_count = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
        #[allow(unused_mut)]
        pub fn sid_count(&mut self, value: usize) -> &mut Self {
            let mut new = self;
            new.sid_count = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
        #[allow(unused_mut)]
        pub fn iid(&mut self, value: nd::Array1<String>) -> &mut Self {
            let mut new = self;
            new.iid = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
        #[allow(unused_mut)]
        pub fn sid(&mut self, value: nd::Array1<String>) -> &mut Self {
            let mut new = self;
            new.sid = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
        #[allow(unused_mut)]
        pub fn chromosome(&mut self, value: nd::Array1<String>) -> &mut Self {
            let mut new = self;
            new.chromosome = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
    }
    impl ::derive_builder::export::core::default::Default for BedBuilder {
        fn default() -> Self {
            Self {
                filename: ::derive_builder::export::core::default::Default::default(),
                count_a1: ::derive_builder::export::core::default::Default::default(),
                iid_count: ::derive_builder::export::core::default::Default::default(),
                sid_count: ::derive_builder::export::core::default::Default::default(),
                iid: ::derive_builder::export::core::default::Default::default(),
                sid: ::derive_builder::export::core::default::Default::default(),
                chromosome: ::derive_builder::export::core::default::Default::default(),
            }
        }
    }
    ///Error type for BedBuilder
    #[non_exhaustive]
    pub enum BedBuilderError {
        /// Uninitialized field
        UninitializedField(&'static str),
        /// Custom validation error
        ValidationError(::derive_builder::export::core::string::String),
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for BedBuilderError {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match (&*self,) {
                (&BedBuilderError::UninitializedField(ref __self_0),) => {
                    let debug_trait_builder =
                        &mut ::core::fmt::Formatter::debug_tuple(f, "UninitializedField");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
                (&BedBuilderError::ValidationError(ref __self_0),) => {
                    let debug_trait_builder =
                        &mut ::core::fmt::Formatter::debug_tuple(f, "ValidationError");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
            }
        }
    }
    impl ::derive_builder::export::core::convert::From<::derive_builder::UninitializedFieldError>
        for BedBuilderError
    {
        fn from(s: ::derive_builder::UninitializedFieldError) -> Self {
            Self::UninitializedField(s.field_name())
        }
    }
    impl
        ::derive_builder::export::core::convert::From<
            ::derive_builder::export::core::string::String,
        > for BedBuilderError
    {
        fn from(s: ::derive_builder::export::core::string::String) -> Self {
            Self::ValidationError(s)
        }
    }
    impl ::derive_builder::export::core::fmt::Display for BedBuilderError {
        fn fmt(
            &self,
            f: &mut ::derive_builder::export::core::fmt::Formatter,
        ) -> ::derive_builder::export::core::fmt::Result {
            match self {
                Self::UninitializedField(ref field) => f.write_fmt(::core::fmt::Arguments::new_v1(
                    &["`", "` must be initialized"],
                    &[::core::fmt::ArgumentV1::new_display(&field)],
                )),
                Self::ValidationError(ref error) => f.write_fmt(::core::fmt::Arguments::new_v1(
                    &[""],
                    &[::core::fmt::ArgumentV1::new_display(&error)],
                )),
            }
        }
    }
    impl std::error::Error for BedBuilderError {}
    impl Bed {
        pub fn builder(filename: String) -> BedBuilder {
            BedBuilder::new(filename)
        }
        pub fn read_options<TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing>(
            &self,
        ) -> ReadOptionBuilder<TOut> {
            ReadOptionBuilder::default()
        }
    }
    impl BedBuilder {
        pub fn new(filename: String) -> Self {
            Self {
                filename: Some(filename),
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
                filename: match self.filename {
                    Some(ref value) => Clone::clone(value),
                    None => {
                        ::core::panicking::panic_fmt(::core::fmt::Arguments::new_v1(
                            &["filename is required"],
                            &[],
                        ));
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
            let filename = bed.filename.to_string();
            let mut buf_reader = BufReader::new(File::open(&filename).unwrap());
            let mut bytes_vector: Vec<u8> = ::alloc::vec::from_elem(0, CB_HEADER_USIZE);
            buf_reader.read_exact(&mut bytes_vector).unwrap();
            if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
                return Err(BedError::IllFormed(filename.to_string()).into());
            }
            Result::Ok(bed)
        }
    }
    impl Bed {
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
            let path_buf = Path::new(&self.filename).with_extension(suffix);
            let file = match File::open(&path_buf) {
                Err(_) => {
                    let string_path = path_buf.to_string_lossy().to_string();
                    return Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                        string_path,
                    )));
                }
                Ok(file) => file,
            };
            let reader = BufReader::new(file);
            let iter = reader.lines().map(move |line| {
                let line = line.unwrap();
                let field = line.split_whitespace().nth(field_index).unwrap();
                field.to_string()
            });
            Ok(iter)
        }
        pub fn get_iid(&mut self) -> &nd::Array1<String> {
            if self.iid.is_some() {
                self.iid.as_ref().unwrap()
            } else {
                let iid = self.read_fam_or_bim("fam", 1).unwrap().collect();
                self.iid = Some(iid);
                self.iid.as_ref().unwrap()
            }
        }
        pub fn get_sid(&mut self) -> &nd::Array1<String> {
            if self.sid.is_some() {
                self.sid.as_ref().unwrap()
            } else {
                let sid = self.read_fam_or_bim("bim", 1).unwrap().collect();
                self.sid = Some(sid);
                self.sid.as_ref().unwrap()
            }
        }
        pub fn get_chromosome(&mut self) -> &nd::Array1<String> {
            if self.chromosome.is_some() {
                self.chromosome.as_ref().unwrap()
            } else {
                let chromosome = self.read_fam_or_bim("bim", 0).unwrap().collect();
                self.chromosome = Some(chromosome);
                self.chromosome.as_ref().unwrap()
            }
        }
        pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send + Missing>(
            &self,
            read_arg: ReadOption<TOut>,
        ) -> Result<nd::Array2<TOut>, BedErrorPlus> {
            let (iid_count, sid_count) = counts(&self.filename)?;
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
    pub fn to_vec(index: Index, count: usize) -> Vec<usize> {
        match index {
            Index::None => (0..count).collect(),
            Index::Vec(iid_index) => iid_index,
            Index::Bool(bool_index) => (0..count)
                .zip(bool_index)
                .filter(|(_, b)| *b)
                .map(|(i, _)| i)
                .collect(),
            Index::NDSlice(slice_index) => {
                let full_array: nd::Array1<usize> = (0..count).collect();
                let array = full_array.slice(slice_index);
                array.to_vec()
            }
        }
    }
    pub(crate) type SliceInfo1 =
        nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>;
    pub enum Index {
        None,
        Vec(Vec<usize>),
        Bool(nd::Array1<bool>),
        NDSlice(SliceInfo1),
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for Index {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match (&*self,) {
                (&Index::None,) => ::core::fmt::Formatter::write_str(f, "None"),
                (&Index::Vec(ref __self_0),) => {
                    let debug_trait_builder = &mut ::core::fmt::Formatter::debug_tuple(f, "Vec");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
                (&Index::Bool(ref __self_0),) => {
                    let debug_trait_builder = &mut ::core::fmt::Formatter::debug_tuple(f, "Bool");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
                (&Index::NDSlice(ref __self_0),) => {
                    let debug_trait_builder =
                        &mut ::core::fmt::Formatter::debug_tuple(f, "NDSlice");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for Index {
        #[inline]
        fn clone(&self) -> Index {
            match (&*self,) {
                (&Index::None,) => Index::None,
                (&Index::Vec(ref __self_0),) => {
                    Index::Vec(::core::clone::Clone::clone(&(*__self_0)))
                }
                (&Index::Bool(ref __self_0),) => {
                    Index::Bool(::core::clone::Clone::clone(&(*__self_0)))
                }
                (&Index::NDSlice(ref __self_0),) => {
                    Index::NDSlice(::core::clone::Clone::clone(&(*__self_0)))
                }
            }
        }
    }
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
    pub struct ReadOption<'a, TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing> {
        bed: &'a Bed,
        #[builder(default = "TOut::missing()")]
        missing_value: TOut,
        #[builder(default = "Index::None")]
        iid_index: Index,
        #[builder(default = "Index::None")]
        sid_index: Index,
        #[builder(default = "true")]
        output_is_orderf: bool,
        #[builder(default, setter(strip_option))]
        pub num_threads: Option<usize>,
    }
    #[allow(clippy::all)]
    ///Builder for [`ReadOption`](struct.ReadOption.html).
    pub struct ReadOptionBuilder<
        'a,
        TOut: Copy + Default + From<i8> + Debug + Sync + Send + Missing,
    > {
        bed: ::derive_builder::export::core::option::Option<&'a Bed>,
        missing_value: ::derive_builder::export::core::option::Option<TOut>,
        iid_index: ::derive_builder::export::core::option::Option<Index>,
        sid_index: ::derive_builder::export::core::option::Option<Index>,
        output_is_orderf: ::derive_builder::export::core::option::Option<bool>,
        num_threads: ::derive_builder::export::core::option::Option<Option<usize>>,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    #[allow(clippy::all)]
    impl<
            'a,
            TOut: ::core::clone::Clone + Copy + Default + From<i8> + Debug + Sync + Send + Missing,
        > ::core::clone::Clone for ReadOptionBuilder<'a, TOut>
    {
        #[inline]
        fn clone(&self) -> ReadOptionBuilder<'a, TOut> {
            match *self {
                ReadOptionBuilder {
                    bed: ref __self_0_0,
                    missing_value: ref __self_0_1,
                    iid_index: ref __self_0_2,
                    sid_index: ref __self_0_3,
                    output_is_orderf: ref __self_0_4,
                    num_threads: ref __self_0_5,
                } => ReadOptionBuilder {
                    bed: ::core::clone::Clone::clone(&(*__self_0_0)),
                    missing_value: ::core::clone::Clone::clone(&(*__self_0_1)),
                    iid_index: ::core::clone::Clone::clone(&(*__self_0_2)),
                    sid_index: ::core::clone::Clone::clone(&(*__self_0_3)),
                    output_is_orderf: ::core::clone::Clone::clone(&(*__self_0_4)),
                    num_threads: ::core::clone::Clone::clone(&(*__self_0_5)),
                },
            }
        }
    }
    #[allow(clippy::all)]
    #[allow(dead_code)]
    impl<
            'a,
            TOut: Copy
                + Default
                + From<i8>
                + Debug
                + Sync
                + Send
                + Missing
                + ::derive_builder::export::core::clone::Clone,
        > ReadOptionBuilder<'a, TOut>
    {
        #[allow(unused_mut)]
        pub fn bed(&mut self, value: &'a Bed) -> &mut Self {
            let mut new = self;
            new.bed = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn missing_value(&mut self, value: TOut) -> &mut Self {
            let mut new = self;
            new.missing_value = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn iid_index(&mut self, value: Index) -> &mut Self {
            let mut new = self;
            new.iid_index = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn sid_index(&mut self, value: Index) -> &mut Self {
            let mut new = self;
            new.sid_index = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn output_is_orderf(&mut self, value: bool) -> &mut Self {
            let mut new = self;
            new.output_is_orderf = ::derive_builder::export::core::option::Option::Some(value);
            new
        }
        #[allow(unused_mut)]
        pub fn num_threads(&mut self, value: usize) -> &mut Self {
            let mut new = self;
            new.num_threads = ::derive_builder::export::core::option::Option::Some(
                ::derive_builder::export::core::option::Option::Some(value),
            );
            new
        }
        ///Builds a new `ReadOption`.
        ///
        ///# Errors
        ///
        ///If a required field has not been initialized.
        pub fn build(
            &self,
        ) -> ::derive_builder::export::core::result::Result<
            ReadOption<'a, TOut>,
            ReadOptionBuilderError,
        > {
            Ok(ReadOption {
                bed: match self.bed {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => {
                        return ::derive_builder::export::core::result::Result::Err(
                            ::derive_builder::export::core::convert::Into::into(
                                ::derive_builder::UninitializedFieldError::from("bed"),
                            ),
                        )
                    }
                },
                missing_value: match self.missing_value {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => TOut::missing(),
                },
                iid_index: match self.iid_index {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => Index::None,
                },
                sid_index: match self.sid_index {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => Index::None,
                },
                output_is_orderf: match self.output_is_orderf {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => true,
                },
                num_threads: match self.num_threads {
                    Some(ref value) => ::derive_builder::export::core::clone::Clone::clone(value),
                    None => ::std::default::Default::default(),
                },
            })
        }
    }
    impl<
            'a,
            TOut: Copy
                + Default
                + From<i8>
                + Debug
                + Sync
                + Send
                + Missing
                + ::derive_builder::export::core::clone::Clone,
        > ::derive_builder::export::core::default::Default for ReadOptionBuilder<'a, TOut>
    {
        fn default() -> Self {
            Self {
                bed: ::derive_builder::export::core::default::Default::default(),
                missing_value: ::derive_builder::export::core::default::Default::default(),
                iid_index: ::derive_builder::export::core::default::Default::default(),
                sid_index: ::derive_builder::export::core::default::Default::default(),
                output_is_orderf: ::derive_builder::export::core::default::Default::default(),
                num_threads: ::derive_builder::export::core::default::Default::default(),
            }
        }
    }
    ///Error type for ReadOptionBuilder
    #[non_exhaustive]
    pub enum ReadOptionBuilderError {
        /// Uninitialized field
        UninitializedField(&'static str),
        /// Custom validation error
        ValidationError(::derive_builder::export::core::string::String),
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ReadOptionBuilderError {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match (&*self,) {
                (&ReadOptionBuilderError::UninitializedField(ref __self_0),) => {
                    let debug_trait_builder =
                        &mut ::core::fmt::Formatter::debug_tuple(f, "UninitializedField");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
                (&ReadOptionBuilderError::ValidationError(ref __self_0),) => {
                    let debug_trait_builder =
                        &mut ::core::fmt::Formatter::debug_tuple(f, "ValidationError");
                    let _ = ::core::fmt::DebugTuple::field(debug_trait_builder, &&(*__self_0));
                    ::core::fmt::DebugTuple::finish(debug_trait_builder)
                }
            }
        }
    }
    impl ::derive_builder::export::core::convert::From<::derive_builder::UninitializedFieldError>
        for ReadOptionBuilderError
    {
        fn from(s: ::derive_builder::UninitializedFieldError) -> Self {
            Self::UninitializedField(s.field_name())
        }
    }
    impl
        ::derive_builder::export::core::convert::From<
            ::derive_builder::export::core::string::String,
        > for ReadOptionBuilderError
    {
        fn from(s: ::derive_builder::export::core::string::String) -> Self {
            Self::ValidationError(s)
        }
    }
    impl ::derive_builder::export::core::fmt::Display for ReadOptionBuilderError {
        fn fmt(
            &self,
            f: &mut ::derive_builder::export::core::fmt::Formatter,
        ) -> ::derive_builder::export::core::fmt::Result {
            match self {
                Self::UninitializedField(ref field) => f.write_fmt(::core::fmt::Arguments::new_v1(
                    &["`", "` must be initialized"],
                    &[::core::fmt::ArgumentV1::new_display(&field)],
                )),
                Self::ValidationError(ref error) => f.write_fmt(::core::fmt::Arguments::new_v1(
                    &[""],
                    &[::core::fmt::ArgumentV1::new_display(&error)],
                )),
            }
        }
    }
    impl std::error::Error for ReadOptionBuilderError {}
}
