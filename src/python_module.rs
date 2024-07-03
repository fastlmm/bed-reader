#![cfg(feature = "extension-module")]

use crate::{BedCloud, CloudFile};
use crate::{
    BedError, BedErrorPlus, Dist, _file_ata_piece_internal, create_pool, encode1, file_aat_piece,
    file_ata_piece, file_b_less_aatbx, impute_and_zero_mean_snps, matrix_subset_no_alloc,
    read_into_f32, read_into_f64, Bed, ReadOptions, WriteOptions,
};
use numpy::PyArrayMethods;
use numpy::{PyArray1, PyArray2, PyArray3};
use pyo3::{
    exceptions::PyIOError,
    exceptions::PyIndexError,
    exceptions::PyValueError,
    prelude::*,
    types::{PyDict, PyModule},
    Bound, PyErr,
};
use std::collections::HashMap;
use tokio::runtime;

#[pymodule]
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn bed_reader(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // See User's guide: https://pyo3.rs/v0.15.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.13.1/exception.html

    mod io {
        pyo3::import_exception!(io, UnsupportedOperation);
    }

    impl std::convert::From<Box<BedErrorPlus>> for PyErr {
        fn from(err: Box<BedErrorPlus>) -> PyErr {
            match *err {
                BedErrorPlus::BedError(
                    BedError::IidIndexTooBig(_)
                    | BedError::SidIndexTooBig(_)
                    | BedError::IndexMismatch(_, _, _, _)
                    | BedError::IndexesTooBigForFiles(_, _)
                    | BedError::SubsetMismatch(_, _, _, _),
                ) => PyIndexError::new_err(err.to_string()),

                BedErrorPlus::IOError(_) => PyIOError::new_err(err.to_string()),

                _ => PyValueError::new_err(err.to_string()),
            }
        }
    }

    #[pyfn(m)]
    fn url_to_bytes(location: &str, options: &Bound<'_, PyDict>) -> Result<Vec<u8>, PyErr> {
        let options: HashMap<String, String> = options.extract()?;
        let cloud_file = CloudFile::new_with_options(location, options)
            .map_err(|e| Box::new(BedErrorPlus::CloudFileError(e)))?;
        let rt = runtime::Runtime::new()?;
        rt.block_on(async {
            let all = cloud_file.read_all().await.map_err(|e| {
                PyErr::new::<PyValueError, _>(format!(
                    "Error retrieving bytes for '{location}: {e}"
                ))
            })?;
            let vec: Vec<u8> = all.to_vec();
            Ok(vec)
        })
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_pass_by_value)]
    fn read_f64(
        filename: &str,
        _ignore: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let mut bed = Bed::builder(filename)
            .iid_count(iid_count)
            .sid_count(sid_count)
            .build()?;

        ReadOptions::builder()
            .iid_index(*ii)
            .sid_index(*si)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .read_and_fill(&mut bed, &mut val.view_mut())?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_pass_by_value)]
    fn read_f32(
        filename: &str,
        _ignore: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let mut bed = Bed::builder(filename)
            .iid_count(iid_count)
            .sid_count(sid_count)
            .build()?;

        ReadOptions::builder()
            .iid_index(*ii)
            .sid_index(*si)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .read_and_fill(&mut bed, &mut val.view_mut())?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_pass_by_value)]
    fn read_i8(
        filename: &str,
        _ignore: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<i8>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let mut bed = Bed::builder(filename)
            .iid_count(iid_count)
            .sid_count(sid_count)
            .build()?;

        ReadOptions::builder()
            .iid_index(*ii)
            .sid_index(*si)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .read_and_fill(&mut bed, &mut val.view_mut())?;

        Ok(())
    }

    #[pyfn(m)]
    fn check_file_cloud(location: &str, options: &Bound<'_, PyDict>) -> Result<(), PyErr> {
        let options: HashMap<String, String> = options.extract()?;
        runtime::Runtime::new()?.block_on(async {
            BedCloud::new_with_options(location, options).await?;
            Ok(())
        })
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn read_cloud_i8(
        location: &str,
        options: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<i8>>,
        num_threads: usize,
        max_concurrent_requests: usize,
        max_chunk_bytes: usize,
    ) -> Result<(), PyErr> {
        let options: HashMap<String, String> = options.extract()?;

        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let rt = runtime::Runtime::new()?;
        rt.block_on(async {
            let mut bed_cloud = BedCloud::builder_with_options(location, options)?
                .iid_count(iid_count)
                .sid_count(sid_count)
                .build()
                .await?;

            ReadOptions::builder()
                .iid_index(*ii)
                .sid_index(*si)
                .is_a1_counted(is_a1_counted)
                .num_threads(num_threads)
                .max_concurrent_requests(max_concurrent_requests)
                .max_chunk_bytes(max_chunk_bytes)
                .read_and_fill_cloud(&mut bed_cloud, &mut val.view_mut())
                .await?;

            Ok(())
        })
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn read_cloud_f32(
        location: &str,
        options: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
        max_concurrent_requests: usize,
        max_chunk_bytes: usize,
    ) -> Result<(), PyErr> {
        let options: HashMap<String, String> = options.extract()?;

        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let rt = runtime::Runtime::new()?;
        rt.block_on(async {
            let mut bed_cloud = BedCloud::builder_with_options(location, options)?
                .iid_count(iid_count)
                .sid_count(sid_count)
                .build()
                .await?;

            ReadOptions::builder()
                .iid_index(*ii)
                .sid_index(*si)
                .is_a1_counted(is_a1_counted)
                .num_threads(num_threads)
                .max_concurrent_requests(max_concurrent_requests)
                .max_chunk_bytes(max_chunk_bytes)
                .read_and_fill_cloud(&mut bed_cloud, &mut val.view_mut())
                .await?;

            Ok(())
        })
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn read_cloud_f64(
        location: &str,
        options: &Bound<'_, PyDict>,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        val: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
        max_concurrent_requests: usize,
        max_chunk_bytes: usize,
    ) -> Result<(), PyErr> {
        let options: HashMap<String, String> = options.extract()?;

        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = val.readwrite();
        let mut val = val.as_array_mut();

        let rt = runtime::Runtime::new()?;
        rt.block_on(async {
            let mut bed_cloud = BedCloud::builder_with_options(location, options)?
                .iid_count(iid_count)
                .sid_count(sid_count)
                .build()
                .await?;

            ReadOptions::builder()
                .iid_index(*ii)
                .sid_index(*si)
                .is_a1_counted(is_a1_counted)
                .num_threads(num_threads)
                .max_concurrent_requests(max_concurrent_requests)
                .max_chunk_bytes(max_chunk_bytes)
                .read_and_fill_cloud(&mut bed_cloud, &mut val.view_mut())
                .await?;

            Ok(())
        })
    }

    #[pyfn(m)]
    fn write_f64(
        filename: &str,
        is_a1_counted: bool,
        val: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        // LATER: If these methods are later changed
        // to support major="individual", be sure to
        // change python function 'to_bed' which
        // currently uses a work-around.
        let mut val = val.readwrite();
        let val = val.as_array_mut();

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    fn write_f32(
        filename: &str,
        is_a1_counted: bool,
        val: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = val.readwrite();
        let val = val.as_array_mut();

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    fn write_i8(
        filename: &str,
        is_a1_counted: bool,
        val: &Bound<'_, PyArray2<i8>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = val.readwrite();
        let val = val.as_array_mut();

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(unused_variables)]
    fn encode1_i8(
        is_a1_counted: bool,
        val: &Bound<'_, PyArray1<i8>>,
        bytes_vector: &Bound<'_, PyArray1<u8>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = val.readonly();
        let val = val.as_array();
        let mut bytes_vector = bytes_vector.readwrite();
        let mut bytes_vector = bytes_vector.as_array_mut();
        let bytes_vector = bytes_vector
            .as_slice_mut()
            .ok_or_else(|| Box::new(BedError::EncodingContiguous().into()))?;

        Ok(encode1(&val, bytes_vector, is_a1_counted, -127)?)
    }

    #[pyfn(m)]
    #[allow(unused_variables)]
    fn encode1_f32(
        is_a1_counted: bool,
        val: &Bound<'_, PyArray1<f32>>,
        bytes_vector: &Bound<'_, PyArray1<u8>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = val.readonly();
        let val = val.as_array();
        let mut bytes_vector = bytes_vector.readwrite();
        let mut bytes_vector = bytes_vector.as_array_mut();
        let bytes_vector = bytes_vector
            .as_slice_mut()
            .ok_or_else(|| Box::new(BedError::EncodingContiguous().into()))?;

        Ok(encode1(&val, bytes_vector, is_a1_counted, f32::NAN)?)
    }

    #[pyfn(m)]
    #[allow(unused_variables)]
    fn encode1_f64(
        is_a1_counted: bool,
        val: &Bound<'_, PyArray1<f64>>,
        bytes_vector: &Bound<'_, PyArray1<u8>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = val.readonly();
        let val = val.as_array();
        let mut bytes_vector = bytes_vector.readwrite();
        let mut bytes_vector = bytes_vector.as_array_mut();
        let bytes_vector = bytes_vector
            .as_slice_mut()
            .ok_or_else(|| Box::new(BedError::EncodingContiguous().into()))?;

        Ok(encode1(&val, bytes_vector, is_a1_counted, f64::NAN)?)
    }

    #[pyfn(m)]
    fn subset_f64_f64(
        val_in: &Bound<'_, PyArray3<f64>>,
        iid_index: &Bound<'_, PyArray1<usize>>,
        sid_index: &Bound<'_, PyArray1<usize>>,
        val_out: &Bound<'_, PyArray3<f64>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();

        let val_in = val_in.readonly();
        let val_in = val_in.as_array();
        let mut val_out = val_out.readwrite();
        let mut val_out = val_out.as_array_mut();

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    fn subset_f32_f64(
        val_in: &Bound<'_, PyArray3<f32>>,
        iid_index: &Bound<'_, PyArray1<usize>>,
        sid_index: &Bound<'_, PyArray1<usize>>,
        val_out: &Bound<'_, PyArray3<f64>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();

        let val_in = val_in.readonly();
        let val_in = val_in.as_array();
        let mut val_out = val_out.readwrite();
        let mut val_out = val_out.as_array_mut();

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    fn subset_f32_f32(
        val_in: &Bound<'_, PyArray3<f32>>,
        iid_index: &Bound<'_, PyArray1<usize>>,
        sid_index: &Bound<'_, PyArray1<usize>>,
        val_out: &Bound<'_, PyArray3<f32>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();

        let val_in = val_in.readonly();
        let val_in = val_in.as_array();
        let mut val_out = val_out.readwrite();
        let mut val_out = val_out.as_array_mut();

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn standardize_f32(
        val: &Bound<'_, PyArray2<f32>>,
        beta_not_unit_variance: bool,
        beta_a: f64,
        beta_b: f64,
        apply_in_place: bool,
        use_stats: bool,
        stats: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = val.readwrite();
        let mut val = val.as_array_mut();
        let mut stats = stats.readwrite();
        let mut stats = stats.as_array_mut();
        let dist = create_dist(beta_not_unit_variance, beta_a, beta_b);
        create_pool(num_threads)?.install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                &dist,
                apply_in_place,
                use_stats,
                &mut stats.view_mut(),
            )
        })?;
        Ok(())
    }

    fn create_dist(beta_not_unit_variance: bool, a: f64, b: f64) -> Dist {
        if beta_not_unit_variance {
            Dist::Beta { a, b }
        } else {
            Dist::Unit
        }
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn standardize_f64(
        val: &Bound<'_, PyArray2<f64>>,
        beta_not_unit_variance: bool,
        beta_a: f64,
        beta_b: f64,
        apply_in_place: bool,
        use_stats: bool,
        stats: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = val.readwrite();
        let mut val = val.as_array_mut();
        let mut stats = stats.readwrite();
        let mut stats = stats.as_array_mut();
        let dist = create_dist(beta_not_unit_variance, beta_a, beta_b);

        create_pool(num_threads)?.install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                &dist,
                apply_in_place,
                use_stats,
                &mut stats.view_mut(),
            )
        })?;
        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn file_ata_piece_f32_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        col_start: usize,
        ata_piece: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = ata_piece.readwrite();
        let mut ata_piece = ata_piece.as_array_mut();

        create_pool(num_threads)?.install(|| {
            file_ata_piece(
                filename,
                offset,
                row_count,
                col_count,
                col_start,
                &mut ata_piece,
                log_frequency,
                read_into_f32,
            )
        })?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn file_ata_piece_f64_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        col_start: usize,
        ata_piece: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = ata_piece.readwrite();
        let mut ata_piece = ata_piece.as_array_mut();

        create_pool(num_threads)?.install(|| {
            file_ata_piece(
                filename,
                offset,
                row_count,
                col_count,
                col_start,
                &mut ata_piece,
                log_frequency,
                read_into_f64,
            )
        })?;

        Ok(())
    }

    // Old version of function for backwards compatibility
    #[pyfn(m)]
    fn file_dot_piece(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_start: usize,
        ata_piece: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = ata_piece.readwrite();
        let mut ata_piece = ata_piece.as_array_mut();

        create_pool(num_threads)?.install(|| {
            _file_ata_piece_internal(
                filename,
                offset,
                row_count,
                col_start,
                &mut ata_piece,
                log_frequency,
                read_into_f64,
            )
        })?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn file_aat_piece_f32_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        row_start: usize,
        aat_piece: &Bound<'_, PyArray2<f32>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut aat_piece = aat_piece.readwrite();
        let mut aat_piece = aat_piece.as_array_mut();

        create_pool(num_threads)?.install(|| {
            file_aat_piece(
                filename,
                offset,
                row_count,
                col_count,
                row_start,
                &mut aat_piece,
                log_frequency,
                read_into_f32,
            )
        })?;

        Ok(())
    }

    #[pyfn(m)]
    #[allow(clippy::too_many_arguments)]
    fn file_aat_piece_f64_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        row_start: usize,
        aat_piece: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut aat_piece = aat_piece.readwrite();
        let mut aat_piece = aat_piece.as_array_mut();

        create_pool(num_threads)?.install(|| {
            file_aat_piece(
                filename,
                offset,
                row_count,
                col_count,
                row_start,
                &mut aat_piece,
                log_frequency,
                read_into_f64,
            )
        })?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "file_b_less_aatbx")]
    #[allow(clippy::too_many_arguments)]
    fn file_b_less_aatbx_translator(
        a_filename: &str,
        offset: u64,
        iid_count: usize,
        b1: &Bound<'_, PyArray2<f64>>,
        aatb: &Bound<'_, PyArray2<f64>>,
        atb: &Bound<'_, PyArray2<f64>>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut b1 = b1.readwrite();
        let mut b1 = b1.as_array_mut();
        let mut aatb = aatb.readwrite();
        let mut aatb = aatb.as_array_mut();
        let mut atb = atb.readwrite();
        let mut atb = atb.as_array_mut();

        create_pool(num_threads)?.install(|| {
            file_b_less_aatbx(
                a_filename,
                offset,
                iid_count,
                &mut b1,
                &mut aatb,
                &mut atb,
                log_frequency,
            )
        })?;

        Ok(())
    }
    Ok(())
}

// LATER on both rust and python side, when counting bim and fam files, also parse them -- don't read them twice.
