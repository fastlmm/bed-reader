use numpy::{PyArray1, PyArray2, PyArray3};
use pyo3::{
    exceptions::{PyOSError, PyValueError}, // !!!cmk make sure give right Python error for each bed error
    prelude::{pymodule, PyModule, PyResult, Python},
    PyErr,
};

use crate::{
    create_pool, impute_and_zero_mean_snps, matrix_subset_no_alloc, read_no_alloc, write,
    BedErrorPlus,
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
            PyOSError::new_err(err.to_string())
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

        let result = create_pool(num_threads).install(|| {
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
        });

        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
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

        let result = create_pool(num_threads).install(|| {
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
        });

        match result {
            Err(result) => Err(PyValueError::new_err(result.to_string())), // !!!cmk make all ValueError not Os error?
            _ => Ok(()),
        }
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

        let result = create_pool(num_threads).install(|| {
            read_no_alloc(
                filename, iid_count, sid_count, count_a1, ii, si, -127i8, &mut val,
            )
        });

        match result {
            Err(result) => Err(PyValueError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        let result =
            create_pool(num_threads).install(|| write(filename, &val, count_a1, (true, f64::NAN)));

        match result {
            Err(result) => Err(PyValueError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        let result =
            create_pool(num_threads).install(|| write(filename, &val, count_a1, (true, f32::NAN)));

        match result {
            Err(result) => Err(PyValueError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        let result =
            create_pool(num_threads).install(|| write(filename, &val, count_a1, (false, -127)));

        match result {
            Err(result) => Err(PyValueError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        println!("cmk subset6464 {:?}", num_threads);

        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let result = create_pool(num_threads)
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out));

        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        println!("cmk subset3264 {:?}", num_threads);
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let result = create_pool(num_threads)
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out));

        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        println!("cmk subset3232 {:?}", num_threads);
        let iid_index = iid_index.readonly();
        let sid_index = sid_index.readonly();
        let val_in = unsafe { val_in.as_array() };
        let mut val_out = unsafe { val_out.as_array_mut() };

        let ii = &iid_index.as_slice()?;
        let si = &sid_index.as_slice()?;

        let result = create_pool(num_threads)
            .install(|| matrix_subset_no_alloc(&val_in, ii, si, &mut val_out));

        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        let result = create_pool(num_threads).install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                beta_not_unit_variance,
                beta_a,
                beta_b,
                apply_in_place,
                use_stats,
                &mut stats.view_mut(),
            )
        });
        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
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
        println!("cmk std64 {:?}", num_threads);

        let mut val = unsafe { val.as_array_mut() };
        let mut stats = unsafe { stats.as_array_mut() };
        // println!("cmk a");
        let result = create_pool(num_threads).install(|| {
            impute_and_zero_mean_snps(
                &mut val.view_mut(),
                beta_not_unit_variance,
                beta_a,
                beta_b,
                apply_in_place,
                use_stats,
                &mut stats.view_mut(),
            )
        });
        match result {
            Err(result) => Err(PyOSError::new_err(result.to_string())),
            _ => Ok(()),
        }
    }
    Ok(())
}
