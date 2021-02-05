// !!!cmk https://stackoverflow.com/questions/32900809/how-to-suppress-function-is-never-used-warning-for-a-function-used-by-tests

#[cfg(test)]
use crate::try_div_4;
#[cfg(test)]
use crate::{
    counts, impute_and_zero_mean_snps, matrix_subset_no_alloc, read, read_with_indexes, write,
};
#[cfg(test)]
use crate::{internal_read_no_alloc, read_no_alloc, BedError, BedErrorPlus};
#[cfg(test)]
use ndarray as nd;
#[cfg(test)]
use ndarray::ShapeBuilder;
#[cfg(test)]
use ndarray_npy::read_npy;
#[cfg(test)]
use num_traits::{abs, Signed};
#[cfg(test)]
use std::path::Path;

// fn big1() {
//     let bigfile = r"M:\deldir\genbgen\2\merged_487400x220000.1.bed"; // !!!cmk not in datafolder
//                                                                      //slicer = np.s_[4000:6000,:20000]
//                                                                      //with open_bed(bigfile,num_threads=None) as bed:
//     for fortran_order in [false, true].iter() {
//         for dtype in ["f64"].iter() {
//             let start = Instant::now();
//             let iid_index = (4000..6000).collect::<Vec<usize>>();
//             let sid_index = (0..20_000).collect::<Vec<usize>>();
//             let _val1 = read_with_indexes(
//                 bigfile,
//                 &iid_index,
//                 &sid_index,
//                 *fortran_order,
//                 true,
//                 f64::NAN,
//             );
//             println!(
//                 "{},{},{}",
//                 fortran_order,
//                 dtype,
//                 start.elapsed().as_secs_f32()
//             );
//         }
//     }
// }

#[test]
fn best_int8() {
    let path = std::env::current_dir().unwrap();
    println!("cmk{}", path.display());
    let filename = "bed_reader/tests/data/some_missing.bed";

    for output_order_is_f in [true, false].iter() {
        let val = read(filename, *output_order_is_f, true, -127).unwrap();
        let ref_val_i8 = reference_val_i8(true);
        assert_eq!(val, ref_val_i8);
    }
}

#[cfg(test)]
fn reference_val_i8(count_a1: bool) -> nd::Array2<i8> {
    let ref_val = reference_val(count_a1);

    let (row_count, col_count) = ref_val.dim();
    let mut ref_val_i8 = nd::Array2::<i8>::zeros((row_count, col_count));
    for i in 0..row_count {
        for j in 0..col_count {
            // !!!cmk use map?
            if ref_val[[i, j]].is_nan() {
                ref_val_i8[[i, j]] = -127i8;
            } else {
                ref_val_i8[[i, j]] = ref_val[[i, j]] as i8;
            }
        }
    }
    return ref_val_i8;
}

#[test]
fn read1() {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let (iid_count, sid_count) = counts(file).unwrap();
    assert!(iid_count == 10);
    assert!(sid_count == 100);
    let val = read(file, true, true, -127).unwrap();
    let val_f64 = val.mapv(|elem| elem as f64);
    let mean_ = val_f64.mean().unwrap();
    assert!(mean_ == -13.142); // really shouldn't do mean on data where -127 represents missing

    let result = read(
        "bed_reader/tests/data/small_too_short.bed",
        true,
        true,
        -127,
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed)) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn rest_reader_bed() {
    let file = "bed_reader/tests/data/some_missing.bed";
    let count_a1 = false;

    let ref_val = reference_val(count_a1);
    let ref_val_i8 = reference_val_i8(count_a1);

    let val = read(file, true, count_a1, f32::NAN).unwrap();
    assert!(allclose(&ref_val.view(), &val.view(), 1e-08, true));

    let val_f64 = read(file, true, count_a1, f64::NAN).unwrap();
    assert!(allclose(&ref_val.view(), &val_f64.view(), 1e-08, true));

    let val2 = read(file, true, count_a1, -127).unwrap();
    assert_eq!(val2, ref_val_i8);
}

#[cfg(test)]
fn reference_val(count_a1: bool) -> nd::Array2<f64> {
    let file = "bed_reader/tests/data/some_missing.val.npy";

    let mut val: nd::Array2<f64> = read_npy(file).unwrap();
    if !count_a1 {
        val = val * -1.0 + 2.0;
    }

    return val;
}

#[cfg(test)]
fn allclose<
    T1: 'static + Copy + PartialEq + PartialOrd + Signed,
    T2: 'static + Copy + PartialEq + PartialOrd + Signed + Into<T1>,
>(
    val1: &nd::ArrayView2<'_, T1>,
    val2: &nd::ArrayView2<'_, T2>,
    atol: T1,
    equal_nan: bool,
) -> bool {
    assert!(val1.dim() == val2.dim());
    // !!!cmk could be run in parallel
    let result = nd::Zip::from(val1)
        .and(val2)
        .fold(true, |acc, ptr_a, ptr_b| -> bool {
            if !acc {
                return false;
            }
            let a_nan = *ptr_a != *ptr_a;
            let b_nan = *ptr_b != *ptr_b;

            if a_nan || b_nan {
                if equal_nan {
                    return a_nan == b_nan;
                } else {
                    return false;
                }
            } else {
                let c: T1 = abs(*ptr_a - T2::into(*ptr_b));
                return c <= atol;
            }
        });
    return result;
}

#[test]
fn index() {
    let filename = "bed_reader/tests/data/some_missing.bed";
    let (iid_count, sid_count) = counts(filename).unwrap();
    let iid_index_full = (0..iid_count).collect::<Vec<usize>>();
    let ref_val_float = reference_val(true);

    let val = read(filename, true, true, f32::NAN).unwrap();
    assert!(allclose(&ref_val_float.view(), &val.view(), 1e-08, true));

    let val = read_with_indexes(filename, &iid_index_full, &[2], true, true, f32::NAN).unwrap();
    assert!(allclose(
        &(ref_val_float.slice(nd::s![.., 2..3])),
        &val.view(),
        1e-08,
        true
    ));

    let val = read_with_indexes(filename, &[1usize], &[2usize], true, true, f32::NAN).unwrap();
    assert!(allclose(
        &ref_val_float.slice(nd::s![1..2, 2..3]),
        &val.view(),
        1e-08,
        true
    ));

    // val = bed.read([2, -2])
    let val = read_with_indexes(
        filename,
        &iid_index_full,
        &[2, (sid_count - 2)],
        true,
        true,
        f32::NAN,
    )
    .unwrap();

    // !!!cmk look for nd::ndarray macro/function to stack columns
    let col1 = ref_val_float.slice(nd::s![.., 2..3]);
    let col2 = ref_val_float.slice(nd::s![.., -2..-1]);
    assert!(allclose(&col1, &val.slice(nd::s![.., 0..1]), 1e-08, true));
    assert!(allclose(&col2, &val.slice(nd::s![.., 1..2]), 1e-08, true));

    let result = read_with_indexes(filename, &[usize::MAX], &[2], true, true, f32::NAN);
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::IidIndexTooBig)) => (),
        _ => panic!("test failure"),
    };

    let result = read_with_indexes(
        "bed_reader/tests/data/small_no_fam.bed",
        &[0],
        &[0],
        true,
        true,
        f32::NAN,
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::IOError(_)) => (),
        _ => panic!("test failure"),
    };

    let result = read_with_indexes(
        "bed_reader/tests/data/small_no_bim.bed",
        &[0],
        &[0],
        true,
        true,
        f32::NAN,
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::IOError(_)) => (),
        _ => panic!("test failure"),
    };

    let result = read_with_indexes(filename, &[2], &[usize::MAX], true, true, f32::NAN);
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::SidIndexTooBig)) => (),
        _ => panic!("test failure"),
    };

    let mut ignore_val = nd::Array2::zeros((1, 1));
    let result = internal_read_no_alloc(
        &"ignore",
        usize::MAX,
        usize::MAX,
        true,
        &[usize::MAX - 1],
        &[usize::MAX - 1],
        f64::NAN,
        &mut ignore_val.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };

    let result = internal_read_no_alloc(
        &"bed_reader/tests/data/no_such_file.nsf",
        1,
        1,
        true,
        &[0],
        &[0],
        f64::NAN,
        &mut ignore_val.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::IOError(_)) => (),
        _ => panic!("test failure"),
    };
}

// Python has more tests, but they aren't needed in Rust until it gets fancy indexing

#[test]
fn writer() {
    let filename = "bed_reader/tests/data/some_missing.bed";
    let val = read(filename, false, true, -127).unwrap();

    let filename2 = r"m:\deldir\rusttest\rust_bed_reader_writer_test.bed";
    write(filename2, &val.view(), true, (false, -127)).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }
    println!("cmk {}", filename2);

    let val2 = read(filename2, false, true, -127).unwrap();

    println!("cmk {}", val);
    println!("cmk {}", val2);
    assert!(allclose(&val.view(), &val2.view(), 0, true));

    let val = read(filename, false, true, f64::NAN).unwrap();

    let filename2 = r"m:\deldir\rusttest\rust_bed_reader_writer_testf64.bed";
    write(filename2, &val.view(), true, (true, f64::NAN)).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }
    println!("cmk {}", filename2);

    let val2 = read(filename2, false, true, f64::NAN).unwrap();

    println!("cmk {}", val);
    println!("cmk {}", val2);
    assert!(allclose(&val.view(), &val2.view(), 1e-8, true));

    let mut val = read(filename, false, true, f64::NAN).unwrap();
    val[(0, 0)] = 5.0;
    let filename = r"m:\deldir\rusttest\rust_bed_reader_writer_testf64_5.bed";
    let result = write(filename, &val.view(), true, (true, f64::NAN));
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::BadValue)) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn subset1() {
    let in_val1 = nd::arr3(&[
        [[0.0], [1.0], [2.0]],
        [[3.0], [4.0], [5.0]],
        [[6.0], [7.0], [8.0]],
    ]);
    let iid_index = [0usize, 2, 1];
    let sid_index = [2usize, 2, 1, 0];
    let mut out_val1 = nd::Array3::<f32>::zeros((iid_index.len(), sid_index.len(), 1));

    let _ = matrix_subset_no_alloc(
        &in_val1.view(),
        &iid_index,
        &sid_index,
        &mut out_val1.view_mut(),
    )
    .unwrap();

    let answer64 = nd::array![
        [[2.0], [2.0], [1.0], [0.0],],
        [[8.0], [8.0], [7.0], [6.0],],
        [[5.0], [5.0], [4.0], [3.0],]
    ];

    assert_eq!(out_val1, answer64);

    let shape_in = ShapeBuilder::set_f((3, 3, 1), true);
    let mut in_val2 = nd::Array3::<f32>::default(shape_in);
    in_val2.assign(&in_val1);
    let shape_out = ShapeBuilder::set_f((3, 4, 1), true);
    let mut out_val2 = nd::Array3::<f64>::zeros(shape_out);

    let _ = matrix_subset_no_alloc(
        &in_val2.view(),
        &iid_index,
        &sid_index,
        &mut out_val2.view_mut(),
    )
    .unwrap();

    let answer32 = nd::array![
        [[2.0], [2.0], [1.0], [0.0],],
        [[8.0], [8.0], [7.0], [6.0],],
        [[5.0], [5.0], [4.0], [3.0],]
    ];

    assert_eq!(out_val2, answer32);
}

#[test]
fn fill_in() {
    let filename = "bed_reader/tests/data/some_missing.bed";

    for output_is_order_f_ptr in [false, true].iter() {
        let mut val = read(filename, *output_is_order_f_ptr, true, f64::NAN).unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));

        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            false,
            0.0,
            0.0,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();

        assert!((val[(0, 0)] - 0.16783627165933704).abs() < 1e-8);
    }
}

#[test]
fn standardize_unit() {
    for output_is_order_f_ptr in [true, false].iter() {
        let mut val = read(
            r"D:\OneDrive\programs\fastlmm\fastlmm\feature_selection\examples\toydata.5chrom.bed",
            *output_is_order_f_ptr,
            false,
            f64::NAN,
        )
        .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            false,
            0.0,
            0.0,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();

        assert!((val[(0, 0)] - -0.3050261832617668).abs() < 1e-8);
    }
}

#[test]
fn div_4() {
    match try_div_4(0, 16, 0u8) {
        Ok(tup) => assert_eq!(tup, (0usize, 0u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(1, 16, 0u8) {
        Ok(tup) => assert_eq!(tup, (1usize, 1u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(4, 16, 0u8) {
        Ok(tup) => assert_eq!(tup, (1usize, 1u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(5, 16, 0u8) {
        Ok(tup) => assert_eq!(tup, (2usize, 2u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(2000, 0, 0u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };
    match try_div_4(0, 256, 0u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 10, 5u8) {
        Ok(tup) => assert_eq!(tup, (25usize, 25u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(25 * 4 + 1, 10, 5u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 11, 5u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 10, 6u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles)) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn standardize_beta() {
    for output_is_order_f_ptr in [true, false].iter() {
        let mut val = read(
            r"D:\OneDrive\programs\fastlmm\fastlmm\feature_selection\examples\toydata.5chrom.bed",
            *output_is_order_f_ptr,
            false,
            f64::NAN,
        )
        .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            true,
            1.0,
            25.0,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();

        assert!((val[(0, 0)] - -0.000031887380905091765).abs() < 1e-8);
    }
}

#[test]
fn read_errors() {
    let iid_count = 100;
    let sid_count = 200;
    let iid_index = (0..iid_count).collect::<Vec<usize>>();
    let sid_index = (0..iid_count).collect::<Vec<usize>>();
    let output_is_order_f = true;
    let shape = ShapeBuilder::set_f((iid_index.len(), sid_index.len()), output_is_order_f);
    let mut val = nd::Array2::<f64>::default(shape);

    match read_no_alloc(
        &"bed_reader/tests/data/no_such_file.nsf",
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        &mut val.view_mut(),
    ) {
        Err(BedErrorPlus::IOError(_)) => (),
        _ => panic!("test failure"),
    };

    let result = read_no_alloc(
        &"bed_reader/tests/data/some_missing.fam",
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        &mut val.view_mut(),
    );
    // !!!cmk println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed)) => (),
        _ => panic!("test failure"),
    };

    let result = read_no_alloc(
        &"bed_reader/tests/data/empty.bed",
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        &mut val.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::IOError(_)) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn read_modes() {
    let filename = "bed_reader/tests/data/small.bed";
    let (iid_count_s1, sid_count_s1) = counts(filename).unwrap();
    let mut val_small_mode_1 = nd::Array2::<i8>::default((iid_count_s1, sid_count_s1));
    let iid_index_s1 = (0..iid_count_s1).collect::<Vec<usize>>();
    let sid_index_s1 = (0..sid_count_s1).collect::<Vec<usize>>();

    let result = read_no_alloc(
        &filename,
        iid_count_s1,
        sid_count_s1,
        true,
        &iid_index_s1,
        &sid_index_s1,
        -127i8,
        &mut val_small_mode_1.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Ok(_) => (),
        _ => panic!("test failure"),
    };

    let result = read_no_alloc(
        &"bed_reader/tests/data/small_too_short.bed",
        iid_count_s1,
        sid_count_s1,
        true,
        &iid_index_s1,
        &sid_index_s1,
        -127i8,
        &mut val_small_mode_1.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed)) => (),
        _ => panic!("test failure"),
    };

    let mut val_small_mode_0 = nd::Array2::<i8>::default((sid_count_s1, iid_count_s1));
    let result = read_no_alloc(
        &"bed_reader/tests/data/smallmode0.bed",
        sid_count_s1,
        iid_count_s1,
        true,
        &sid_index_s1,
        &iid_index_s1,
        -127i8,
        &mut val_small_mode_0.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Ok(_) => (),
        _ => panic!("test failure"),
    };

    assert_eq!(val_small_mode_0.t(), val_small_mode_1);

    let result = read_no_alloc(
        &"bed_reader/tests/data/smallmodebad.bed",
        iid_count_s1,
        sid_count_s1,
        true,
        &iid_index_s1,
        &sid_index_s1,
        -127i8,
        &mut val_small_mode_1.view_mut(),
    );
    println!("cmk {:?}", result);
    match result {
        Err(BedErrorPlus::BedError(BedError::BadMode)) => (),
        _ => panic!("test failure"),
    };
}

// cmk What does  pyo3::Python::with_gil mean?
// // !!!cmk add more tests from the python/c++ project
// // 2nd half of test_bed_int8
// // test_respect_read_inputs
// // What happens if index is out of range (too big? too small?)
// // Coverage
