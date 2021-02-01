// !!!cmk https://stackoverflow.com/questions/32900809/how-to-suppress-function-is-never-used-warning-for-a-function-used-by-tests

#[cfg(test)]
use crate::{
    counts, impute_and_zero_mean_snps, matrix_subset_no_alloc, read, read_with_indexes, write,
};
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
        let val = read(filename, *output_order_is_f, true, -127);
        let ref_val_i8 = reference_val_i8(true);
        assert_eq!(val, ref_val_i8);
    }
}

#[cfg(test)]
fn reference_val_i8(count_a1: bool) -> nd::Array2<i8> {
    let ref_val = reference_val(count_a1);

    let row_count = ref_val.shape()[0]; // !!!cmk must be a more concise way to do this
    let col_count = ref_val.shape()[1];
    let mut ref_val_i8 = nd::Array2::<i8>::zeros((row_count, col_count));
    for i in 0..ref_val.shape()[0] {
        for j in 0..ref_val.shape()[1] {
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
    let (iid_count, sid_count) = counts(file);
    assert!(iid_count == 10);
    assert!(sid_count == 100);
    let val = read(file, true, true, -127);
    let val_f64 = val.mapv(|elem| elem as f64);
    let mean_ = val_f64.mean().unwrap();
    assert!(mean_ == -13.142); // really shouldn't do mean on data where -127 represents missing
}

#[test]
fn rest_reader_bed() {
    let file = "bed_reader/tests/data/some_missing.bed";
    let count_a1 = false;

    let ref_val = reference_val(count_a1);
    let ref_val_i8 = reference_val_i8(count_a1);

    let val = read(file, true, count_a1, f32::NAN);
    assert!(allclose(&ref_val.view(), &val.view(), 1e-08, true));

    let val_f64 = read(file, true, count_a1, f64::NAN);
    assert!(allclose(&ref_val.view(), &val_f64.view(), 1e-08, true));

    let val2 = read(file, true, count_a1, -127);
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
    assert!(val1.shape() == val2.shape());
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
    let (iid_count, sid_count) = counts(filename);
    let iid_index_full = (0..iid_count).collect::<Vec<usize>>();
    let ref_val_float = reference_val(true);

    let val = read(filename, true, true, f32::NAN);
    assert!(allclose(&ref_val_float.view(), &val.view(), 1e-08, true));

    let val = read_with_indexes(filename, &iid_index_full, &[2], true, true, f32::NAN);
    assert!(allclose(
        &(ref_val_float.slice(nd::s![.., 2..3])),
        &val.view(),
        1e-08,
        true
    ));

    let val = read_with_indexes(filename, &[1usize], &[2usize], true, true, f32::NAN);
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
    );

    // !!!cmk look for nd::ndarray macro/function to stack columns
    let col1 = ref_val_float.slice(nd::s![.., 2..3]);
    let col2 = ref_val_float.slice(nd::s![.., -2..-1]);
    assert!(allclose(&col1, &val.slice(nd::s![.., 0..1]), 1e-08, true));
    assert!(allclose(&col2, &val.slice(nd::s![.., 1..2]), 1e-08, true));

    // Python has more tests, but they aren't needed in Rust until it gets fancy indexing
}

#[test]
fn writer() {
    let filename = "bed_reader/tests/data/some_missing.bed";
    let val = read(filename, false, true, -127);

    let filename2 = r"m:\deldir\rusttest\rust_bed_reader_writer_test.bed";
    write(filename2, &val.view(), true, (false, -127)).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }
    println!("{}", filename2);

    let val2 = read(filename2, false, true, -127);

    println!("{}", val);
    println!("{}", val2);
    assert!(allclose(&val.view(), &val2.view(), 0, true));

    let val = read(filename, false, true, f64::NAN);

    let filename2 = r"m:\deldir\rusttest\rust_bed_reader_writer_testf64.bed";
    write(filename2, &val.view(), true, (true, f64::NAN)).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }
    println!("{}", filename2);

    let val2 = read(filename2, false, true, f64::NAN);

    println!("{}", val);
    println!("{}", val2);
    assert!(allclose(&val.view(), &val2.view(), 1e-8, true));
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
        let mut val = read(filename, *output_is_order_f_ptr, true, f64::NAN);
        let mut stats = nd::Array2::<f64>::zeros((val.shape()[1], 2));

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
        );
        let mut stats = nd::Array2::<f64>::zeros((val.shape()[1], 2));
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
fn standardize_beta() {
    for output_is_order_f_ptr in [true, false].iter() {
        let mut val = read(
            r"D:\OneDrive\programs\fastlmm\fastlmm\feature_selection\examples\toydata.5chrom.bed",
            *output_is_order_f_ptr,
            false,
            f64::NAN,
        );
        let mut stats = nd::Array2::<f64>::zeros((val.shape()[1], 2));
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

// fn c_doer<T: Default + Copy + Sync + Send + Float + ToPrimitive + FromPrimitive>(
//     c: &dyn Fn(T) -> T,
// ) {
//     let input = T::from_f64(3.0).unwrap();
//     println!("{:?}", c(input).to_f64().unwrap());
// }

// #[test]
// fn pass_closure_to_function() {
//     pass_closure_to_function_internal::<f64>();
// }
// fn pass_closure_to_function_internal<
//     T: Default + Copy + Sync + Send + Float + ToPrimitive + FromPrimitive,
// >() {
//     let b = false;
//     let a = T::from_f64(2.0).unwrap();

//     let c0 = |x: T| T::from_f64(21.0).unwrap();
//     let c1 = move |x: T| a + x;

//     let c: &dyn Fn(T) -> T;
//     if b {
//         c = &c0
//     } else {
//         c = &c1
//     }

//     c_doer(&c);
//}

// cmk What does  pyo3::Python::with_gil mean?
// // !!!cmk add more tests from the python/c++ project
// // 2nd half of test_bed_int8
// // test_respect_read_inputs
// // What happens if index is out of range (too big? too small?)
// // Coverage
