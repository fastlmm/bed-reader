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
    read_into_f32, read_into_f64, read_no_alloc, write,
};

#[pymodule]
fn bed_reader(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // See User's guide: https://pyo3.rs/v0.13.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.13.1/exception.html

    mod io {
        pyo3::import_exception!(io, UnsupportedOperation);
    }

    impl std::convert::From<BedErrorPlus> for PyErr {
        fn from(err: BedErrorPlus) -> PyErr {
            match err {
                BedErrorPlus::IOError(_) => PyIOError::new_err(err.to_string()),
                BedErrorPlus::ThreadPoolError(_) => PyValueError::new_err(err.to_string()),
                BedErrorPlus::BedError(BedError::IidIndexTooBig(_))
                | BedErrorPlus::BedError(BedError::SidIndexTooBig(_))
                | BedErrorPlus::BedError(BedError::IndexMismatch(_, _, _, _))
                | BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))
                | BedErrorPlus::BedError(BedError::SubsetMismatch(_, _, _, _)) => {
                    PyIndexError::new_err(err.to_string())
                }
                _ => PyValueError::new_err(err.to_string()),
            }
        }
    }

    #[pyfn(m, "read_f64")]
    fn read_f64_py(
        _py: Python<'_>,
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        count_a1: bool,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val: &PyArray2<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let mut val = unsafe { val.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?.install(|| {
            read_no_alloc(
                filename,
                iid_count,
                sid_count,
                count_a1,
                ii,
                si,
                f64::NAN,
                &mut val,
            )
        })?;

        Ok(())
    }

    #[pyfn(m, "read_f32")]
    fn read_f32_py(
        _py: Python<'_>,
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        count_a1: bool,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val: &PyArray2<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let mut val = unsafe { val.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?.install(|| {
            read_no_alloc(
                filename,
                iid_count,
                sid_count,
                count_a1,
                ii,
                si,
                f32::NAN,
                &mut val,
            )
        })?;

        Ok(())
    }

    #[pyfn(m, "read_i8")]
    fn read_i8_py(
        _py: Python<'_>,
        filename: &str,
        iid_count: usize,
        sid_count: usize,
        count_a1: bool,
        iid_index: &PyArray1<usize>,
        sid_index: &PyArray1<usize>,
        val: &PyArray2<i8>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let mut val = unsafe { val.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        create_pool(num_threads)?.install(|| {
            read_no_alloc(
                filename, iid_count, sid_count, count_a1, ii, si, -127i8, &mut val,
            )
        })?;

        Ok(())
    }

    #[pyfn(m, "write_f64")]
    fn write_f64_py(
        _py: Python<'_>,
        filename: &str,
        count_a1: bool,
        val: &PyArray2<f64>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        create_pool(num_threads)?.install(|| write(filename, &val, count_a1, f64::NAN))?;

        Ok(())
    }

    #[pyfn(m, "write_f32")]
    fn write_f32_py(
        _py: Python<'_>,
        filename: &str,
        count_a1: bool,
        val: &PyArray2<f32>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        create_pool(num_threads)?.install(|| write(filename, &val, count_a1, f32::NAN))?;

        Ok(())
    }

    #[pyfn(m, "write_i8")]
    fn write_i8_py(
        _py: Python<'_>,
        filename: &str,
        count_a1: bool,
        val: &PyArray2<i8>,
        num_threads: usize,
    ) -> Result<(), PyErr> {
        let val = unsafe { val.as_array() };

        create_pool(num_threads)?.install(|| write(filename, &val, count_a1, -127))?;

        Ok(())
    }

    #[pyfn(m, "subset_f64_f64")]
    fn subset_f64_f64(
        _py: Python<'_>,
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

    #[pyfn(m, "subset_f32_f64")]
    fn subset_f32_f64(
        _py: Python<'_>,
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

    #[pyfn(m, "subset_f32_f32")]
    fn subset_f32_f32(
        _py: Python<'_>,
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

    #[pyfn(m, "standardize_f32")]
    fn standardize_f32(
        _py: Python<'_>,
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
            return Dist::Beta { a: a, b: b };
        } else {
            return Dist::Unit;
        };
    }

    #[pyfn(m, "standardize_f64")]
    fn standardize_f64(
        _py: Python<'_>,
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

    #[pyfn(m, "file_ata_piece_f32_orderf")]
    fn file_ata_piece_f32_py(
        _py: Python<'_>,
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

    #[pyfn(m, "file_ata_piece_f64_orderf")]
    fn file_ata_piece_f64_py(
        _py: Python<'_>,
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
    #[pyfn(m, "file_dot_piece")]
    fn file_dot_piece_py(
        _py: Python<'_>,
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

    #[pyfn(m, "file_aat_piece_f32_orderf")]
    fn file_aat_piece_f32_py(
        _py: Python<'_>,
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

    #[pyfn(m, "file_aat_piece_f64_orderf")]
    fn file_aat_piece_f64_py(
        _py: Python<'_>,
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

    #[pyfn(m, "file_b_less_aatbx")]
    fn file_b_less_aatbx_py(
        _py: Python<'_>,
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
