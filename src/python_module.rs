#![cfg(feature = "extension-module")]
use numpy::{PyArray1, PyArray2, PyArray3};

use pyo3::{
    exceptions::PyIOError,
    exceptions::PyIndexError,
    exceptions::PyValueError,
    prelude::{pymodule, PyModule, PyResult, Python},
    PyErr,
};

use crate::{
    BedError, BedErrorPlus, Dist, _file_ata_piece_internal, create_pool, file_aat_piece,
    file_ata_piece, file_b_less_aatbx, impute_and_zero_mean_snps, matrix_subset_no_alloc,
    read_into_f32, read_into_f64, Bed, ReadOptions, WriteOptions,
};

#[pymodule]
fn bed_reader(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // See User's guide: https://pyo3.rs/v0.15.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.13.1/exception.html

    mod io {
        pyo3::import_exception!(io, UnsupportedOperation);
    }

    impl std::convert::From<BedErrorPlus> for PyErr {
        fn from(err: BedErrorPlus) -> PyErr {
            match err {
                BedErrorPlus::BedError(BedError::IidIndexTooBig(_))
                | BedErrorPlus::BedError(BedError::SidIndexTooBig(_))
                | BedErrorPlus::BedError(BedError::IndexMismatch(_, _, _, _))
                | BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))
                | BedErrorPlus::BedError(BedError::SubsetMismatch(_, _, _, _)) => {
                    PyIndexError::new_err(err.to_string())
                }

                BedErrorPlus::IOError(_) => PyIOError::new_err(err.to_string()),

                _ => PyValueError::new_err(err.to_string()),
            }
        }
    }

    #[pyfn(m)]
    #[pyo3(name = "read_f64")]
    #[allow(clippy::too_many_arguments)]
    fn read_f64(
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &PyArray1<isize>,
        sid_index: &PyArray1<isize>,
        val: &PyArray2<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = unsafe { val.as_array_mut() };

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
    #[pyo3(name = "read_f32")]
    #[allow(clippy::too_many_arguments)]
    fn read_f32(
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &PyArray1<isize>,
        sid_index: &PyArray1<isize>,
        val: &PyArray2<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = unsafe { val.as_array_mut() };

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
    #[pyo3(name = "read_i8")]
    #[allow(clippy::too_many_arguments)]
    fn read_i8(
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        is_a1_counted: bool,
        iid_index: &PyArray1<isize>,
        sid_index: &PyArray1<isize>,
        val: &PyArray2<i8>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let mut val = unsafe { val.as_array_mut() };

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
    #[pyo3(name = "write_f64")]
    fn write_f64(
        filename: &str,
        is_a1_counted: bool,
        val: &PyArray2<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "write_f32")]
    fn write_f32(
        filename: &str,
        is_a1_counted: bool,
        val: &PyArray2<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "write_i8")]
    fn write_i8(
        filename: &str,
        is_a1_counted: bool,
        val: &PyArray2<i8>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        WriteOptions::builder(filename)
            .is_a1_counted(is_a1_counted)
            .num_threads(num_threads)
            .skip_fam()
            .skip_bim()
            .write(&val)?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "subset_f64_f64")]
    fn subset_f64_f64(
        val_in: &PyArray3<f64>,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val_out: &PyArray3<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "subset_f32_f64")]
    fn subset_f32_f64(
        val_in: &PyArray3<f32>,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val_out: &PyArray3<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "subset_f32_f32")]
    fn subset_f32_f32(
        val_in: &PyArray3<f32>,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val_out: &PyArray3<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out))?;

        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "standardize_f32")]
    #[allow(clippy::too_many_arguments)]
    fn standardize_f32(
        val: &PyArray2<f32>,
        beta_not_unit_variance: bool,
        beta_a: f64,
        beta_b: f64,
        apply_in_place: bool,
        use_stats: bool,
        stats: &PyArray2<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = unsafe { val.as_array_mut() };
        let mut stats = unsafe { stats.as_array_mut() };
        let dist = create_dist(beta_not_unit_variance, beta_a, beta_b);
        create_pool(num_threads)?.install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                dist,
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
    #[pyo3(name = "standardize_f64")]
    #[allow(clippy::too_many_arguments)]
    fn standardize_f64(
        val: &PyArray2<f64>,
        beta_not_unit_variance: bool,
        beta_a: f64,
        beta_b: f64,
        apply_in_place: bool,
        use_stats: bool,
        stats: &PyArray2<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let mut val = unsafe { val.as_array_mut() };
        let mut stats = unsafe { stats.as_array_mut() };
        let dist = create_dist(beta_not_unit_variance, beta_a, beta_b);

        create_pool(num_threads)?.install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                dist,
                apply_in_place,
                use_stats,
                &mut stats.view_mut(),
            )
        })?;
        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "file_ata_piece_f32_orderf")]
    #[allow(clippy::too_many_arguments)]
    fn file_ata_piece_f32_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        col_start: usize,
        ata_piece: &PyArray2<f32>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = unsafe { ata_piece.as_array_mut() };

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
    #[pyo3(name = "file_ata_piece_f64_orderf")]
    #[allow(clippy::too_many_arguments)]
    fn file_ata_piece_f64_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        col_start: usize,
        ata_piece: &PyArray2<f64>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = unsafe { ata_piece.as_array_mut() };

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
    #[pyo3(name = "file_dot_piece")]
    fn file_dot_piece(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_start: usize,
        ata_piece: &PyArray2<f64>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut ata_piece = unsafe { ata_piece.as_array_mut() };

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
    #[pyo3(name = "file_aat_piece_f32_orderf")]
    #[allow(clippy::too_many_arguments)]
    fn file_aat_piece_f32_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        row_start: usize,
        aat_piece: &PyArray2<f32>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut aat_piece = unsafe { aat_piece.as_array_mut() };

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
    #[pyo3(name = "file_aat_piece_f64_orderf")]
    #[allow(clippy::too_many_arguments)]
    fn file_aat_piece_f64_orderf(
        filename: &str,
        offset: u64,
        row_count: usize,
        col_count: usize,
        row_start: usize,
        aat_piece: &PyArray2<f64>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut aat_piece = unsafe { aat_piece.as_array_mut() };

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
        b1: &PyArray2<f64>,
        aatb: &PyArray2<f64>,
        atb: &PyArray2<f64>,
        num_threads: usize,
        log_frequency: usize,
    ) -> Result<(), PyErr> {
        let mut b1 = unsafe { b1.as_array_mut() };
        let mut aatb = unsafe { aatb.as_array_mut() };
        let mut atb = unsafe { atb.as_array_mut() };

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
