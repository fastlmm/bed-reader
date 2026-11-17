#![cfg(feature = "extension-module")]

use crate::{
    create_pool, encode1, file_aat_piece, file_ata_piece, file_ata_piece_internal,
    file_b_less_aatbx, impute_and_zero_mean_snps, matrix_subset_no_alloc, read_into_f32,
    read_into_f64, Bed, BedError, BedErrorPlus, Dist, ReadOptions, WriteOptions,
};
use crate::{BedCloud, CloudFile};
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

#[pymodule]
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn bed_reader(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // See User's guide: https://pyo3.rs/v0.15.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.13.1/exception.html

    mod io {
        pyo3::import_exception!(io, UnsupportedOperation);
    }

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
    fn check_file_cloud(location: &str, options: &Bound<'_, PyDict>) -> Result<(), PyErr> {
        let options: HashMap<String, String> = options.extract()?;
        runtime::Runtime::new()?.block_on(async {
            BedCloud::new_with_options(location, options).await?;
            Ok(())
        })
    }

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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
    #[pyfunction]
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
            file_ata_piece_internal(
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

    #[pyfunction]
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

    #[pyfunction]
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

    #[pyfunction]
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

    m.add_function(wrap_pyfunction!(url_to_bytes, m)?)?;
    m.add_function(wrap_pyfunction!(read_f64, m)?)?;
    m.add_function(wrap_pyfunction!(read_f32, m)?)?;
    m.add_function(wrap_pyfunction!(read_i8, m)?)?;
    m.add_function(wrap_pyfunction!(check_file_cloud, m)?)?;
    m.add_function(wrap_pyfunction!(read_cloud_i8, m)?)?;
    m.add_function(wrap_pyfunction!(read_cloud_f32, m)?)?;
    m.add_function(wrap_pyfunction!(read_cloud_f64, m)?)?;
    m.add_function(wrap_pyfunction!(write_f64, m)?)?;
    m.add_function(wrap_pyfunction!(write_f32, m)?)?;
    m.add_function(wrap_pyfunction!(write_i8, m)?)?;
    m.add_function(wrap_pyfunction!(encode1_i8, m)?)?;
    m.add_function(wrap_pyfunction!(encode1_f32, m)?)?;
    m.add_function(wrap_pyfunction!(encode1_f64, m)?)?;
    m.add_function(wrap_pyfunction!(subset_f64_f64, m)?)?;
    m.add_function(wrap_pyfunction!(subset_f32_f64, m)?)?;
    m.add_function(wrap_pyfunction!(subset_f32_f32, m)?)?;
    m.add_function(wrap_pyfunction!(standardize_f32, m)?)?;
    m.add_function(wrap_pyfunction!(standardize_f64, m)?)?;
    m.add_function(wrap_pyfunction!(file_ata_piece_f32_orderf, m)?)?;
    m.add_function(wrap_pyfunction!(file_ata_piece_f64_orderf, m)?)?;
    m.add_function(wrap_pyfunction!(file_dot_piece, m)?)?;
    m.add_function(wrap_pyfunction!(file_aat_piece_f32_orderf, m)?)?;
    m.add_function(wrap_pyfunction!(file_aat_piece_f64_orderf, m)?)?;
    m.add_function(wrap_pyfunction!(file_b_less_aatbx_translator, m)?)?;
    #[cfg(feature = "cubecl-cuda")]
    gpu_python::register(m)?;

    Ok(())
}

// LATER on both rust and python side, when counting bim and fam files, also parse them -- don't read them twice.

/// Python bindings for GPU reading (Cargo feature `cubecl-cuda`): `read_gpu` returns a `GpuArray` that `CuPy`
/// (or any `__cuda_array_interface__` consumer) wraps without copying; `standardize_gpu`
/// standardizes it in place first. Low-level plumbing, like `read_f32`: `open_bed.read(xp=...)` and
/// `PySnpTools` call these.
#[cfg(feature = "cubecl-cuda")]
mod gpu_python {
    use crate::backend::{
        standardize::{standardize_beta, standardize_unit},
        CubeCl, CubeClMatrix, CudaPointer,
    };
    use crate::{Bed, ReadOptions};
    use cudarc::driver::CudaContext;
    use std::sync::{Arc, Mutex};
    use numpy::{PyArray1, PyArrayMethods, ToPyArray};
    use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

    enum AnyMatrix {
        I8(CubeClMatrix<i8>),
        F32(CubeClMatrix<f32>),
        F64(CubeClMatrix<f64>),
    }

    /// A matrix in GPU memory. `CuPy` arrays made from it hold a reference to it, which keeps
    /// the GPU buffer alive.
    #[pyclass(weakref)]
    struct GpuArray {
        matrix: AnyMatrix,
        /// Keeps the allocation in place while consumers use the pointer it reports.
        // CudaPointer is Send, but its erased resource is not Sync. A mutex
        // makes the owner transferable between Python threads without unsafe.
        // Access is always through &mut self (also guarded by PyO3's mutable
        // borrow), using get_mut; this mutex is never locked or poisoned.
        pointer: Mutex<Option<CudaPointer>>,
        // CubeCl selects CUDA device 0 and its primary context. Keep the same
        // context alive so destruction can wait for consumers on any stream.
        context: Arc<CudaContext>,
    }

    impl Drop for GpuArray {
        fn drop(&mut self) {
            if self.pointer.get_mut().expect("pointer mutex is never locked").is_some() {
                // CuPy enqueues work asynchronously and can drop its last Python
                // reference before that work completes. CubeCL cannot see those
                // streams: wait before its pool can recycle the shared allocation.
                self.context.synchronize().expect(
                    "CUDA failed while waiting for array consumers before releasing GPU memory",
                );
            }
        }
    }

    #[pymethods]
    impl GpuArray {
        #[getter]
        fn shape(&self) -> (usize, usize) {
            match &self.matrix {
                AnyMatrix::I8(matrix) => matrix.dim(),
                AnyMatrix::F32(matrix) => matrix.dim(),
                AnyMatrix::F64(matrix) => matrix.dim(),
            }
        }

        #[getter]
        fn __cuda_array_interface__<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
            // `cuda_pointer` finishes the matrix's pending GPU work before another library reads it.
            let (dim, is_f, itemsize, typestr, pointer) = match &self.matrix {
                AnyMatrix::I8(matrix) => (matrix.dim(), matrix.is_f(), 1usize, "|i1", matrix.cuda_pointer()?),
                AnyMatrix::F32(matrix) => (matrix.dim(), matrix.is_f(), 4usize, "<f4", matrix.cuda_pointer()?),
                AnyMatrix::F64(matrix) => (matrix.dim(), matrix.is_f(), 8usize, "<f8", matrix.cuda_pointer()?),
            };
            let (rows, cols) = dim;
            let strides = if is_f { (itemsize, itemsize * rows) } else { (itemsize * cols, itemsize) };
            // Even an empty matrix owns a small native allocation. CuPy needs
            // its non-null pointer to determine the device for an empty view.
            let ptr = pointer.ptr();
            *self.pointer.get_mut().expect("pointer mutex is never locked") = Some(pointer);
            let interface = PyDict::new(py);
            interface.set_item("shape", (rows, cols))?;
            interface.set_item("typestr", typestr)?;
            interface.set_item("data", (ptr, false))?;
            interface.set_item("strides", strides)?;
            interface.set_item("version", 3)?;
            Ok(interface)
        }
    }

    #[pyfunction]
    #[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
    fn read_gpu(
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &Bound<'_, PyArray1<isize>>,
        sid_index: &Bound<'_, PyArray1<isize>>,
        dtype: &str,
        is_f: bool,
    ) -> PyResult<GpuArray> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let (ii, si) = (iid_index.as_slice()?, sid_index.as_slice()?);
        let mut bed = Bed::builder(filename).iid_count(iid_count).sid_count(sid_count).build()?;
        let matrix = match dtype {
            "int8" => AnyMatrix::I8(
                ReadOptions::<i8>::builder()
                    .iid_index(ii)
                    .sid_index(si)
                    .is_a1_counted(is_a1_counted)
                    .is_f(is_f)
                    .cubecl()
                    .read(&mut bed)?,
            ),
            "float32" => AnyMatrix::F32(
                ReadOptions::<f32>::builder()
                    .iid_index(ii)
                    .sid_index(si)
                    .is_a1_counted(is_a1_counted)
                    .is_f(is_f)
                    .cubecl()
                    .read(&mut bed)?,
            ),
            "float64" => AnyMatrix::F64(
                ReadOptions::<f64>::builder()
                    .iid_index(ii)
                    .sid_index(si)
                    .is_a1_counted(is_a1_counted)
                    .is_f(is_f)
                    .cubecl()
                    .read(&mut bed)?,
            ),
            _ => return Err(PyValueError::new_err(format!("dtype must be 'int8', 'float32', or 'float64', not '{dtype}'"))),
        };
        let context = CudaContext::new(0).map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(GpuArray { matrix, pointer: Mutex::new(None), context })
    }

    #[pyfunction]
    fn allocate_gpu(rows: usize, columns: usize, dtype: &str, is_f: bool) -> PyResult<GpuArray> {
        let matrix = match dtype {
            "int8" => AnyMatrix::I8(CubeClMatrix::empty(rows, columns, is_f)?),
            "float32" => AnyMatrix::F32(CubeClMatrix::empty(rows, columns, is_f)?),
            "float64" => AnyMatrix::F64(CubeClMatrix::empty(rows, columns, is_f)?),
            _ => return Err(PyValueError::new_err("GPU dtype must be int8, float32, or float64")),
        };
        let context = CudaContext::new(0).map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(GpuArray { matrix, pointer: Mutex::new(None), context })
    }

    /// Standardizes `array` in place on the GPU and returns each SNP's (mean, std), like
    /// `standardize_f32/f64` with `apply_in_place=True, use_stats=False`.
    #[pyfunction]
    fn standardize_gpu<'py>(
        py: Python<'py>,
        array: &mut GpuArray,
        beta_not_unit_variance: bool,
        beta_a: f64,
        beta_b: f64,
    ) -> PyResult<Bound<'py, PyAny>> {
        // A caller may have changed the exported CuPy view on any CUDA stream.
        // Finish those writes before CubeCL reads the shared allocation.
        if array.pointer.get_mut().expect("pointer mutex is never locked").is_some() {
            array.context.synchronize().map_err(|error| PyValueError::new_err(error.to_string()))?;
        }
        let stats = match &mut array.matrix {
            AnyMatrix::I8(_) => return Err(PyValueError::new_err("Standardizing needs float32 or float64, not int8")),
            AnyMatrix::F32(matrix) => {
                let stats = if beta_not_unit_variance { standardize_beta(matrix, beta_a, beta_b)? } else { standardize_unit(matrix)? };
                stats.to_pyarray(py).into_any()
            }
            AnyMatrix::F64(matrix) => {
                let stats = if beta_not_unit_variance { standardize_beta(matrix, beta_a, beta_b)? } else { standardize_unit(matrix)? };
                stats.to_pyarray(py).into_any()
            }
        };
        if array.pointer.get_mut().expect("pointer mutex is never locked").is_some() {
            // Existing CuPy views do not ask for the array interface again.
            // Make queued writes visible before returning to those consumers.
            let client = CubeCl::client()?;
            pollster::block_on(client.sync()).map_err(|error| PyValueError::new_err(error.to_string()))?;
        }
        Ok(stats)
    }

    /// Asks `CubeCL` to return the GPU memory of dropped arrays, so that other GPU libraries in the
    /// process (for example, `CuPy`) can use it.
    #[pyfunction]
    fn gpu_memory_cleanup() -> PyResult<()> {
        CubeCl::memory_cleanup()?;
        Ok(())
    }

    pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<GpuArray>()?;
        m.add_function(wrap_pyfunction!(read_gpu, m)?)?;
        m.add_function(wrap_pyfunction!(allocate_gpu, m)?)?;
        m.add_function(wrap_pyfunction!(standardize_gpu, m)?)?;
        m.add_function(wrap_pyfunction!(gpu_memory_cleanup, m)?)?;
        Ok(())
    }
}
