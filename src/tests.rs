// https://stackoverflow.com/questions/32900809/how-to-suppress-function-is-never-used-warning-for-a-function-used-by-tests

#[cfg(test)]
use crate::file_aat_piece;
#[cfg(test)]
use crate::file_ata_piece;
#[cfg(test)]
use crate::file_b_less_aatbx;
#[cfg(test)]
use crate::read_into_f64;
#[cfg(test)]
use crate::try_div_4;
#[cfg(test)]
use crate::Dist;
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
use std::f32;
#[cfg(test)]
use std::f64;
#[cfg(test)]
use std::io::{LineWriter, Write};
#[cfg(test)]
use std::ops::Range;
#[cfg(test)]
use std::path::Path;
#[cfg(test)]
use std::path::PathBuf;
#[cfg(test)]
use std::{fs::File, io::BufReader};
#[cfg(test)]
use temp_testdir::TempDir;

#[test]
fn best_int8() {
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
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed(_))) => (),
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
    // Could be run in parallel
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

    let col0 = ref_val_float.slice(nd::s![.., 2]);
    let col1 = ref_val_float.slice(nd::s![.., -2]);
    let expected = nd::stack![nd::Axis(1), col0, col1];
    assert!(allclose(&expected.view(), &val.view(), 1e-08, true));

    let result = read_with_indexes(filename, &[usize::MAX], &[2], true, true, f32::NAN);
    match result {
        Err(BedErrorPlus::BedError(BedError::IidIndexTooBig(_))) => (),
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
    match result {
        Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(_))) => (),
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
    match result {
        Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(_))) => (),
        _ => panic!("test failure"),
    };

    let result = read_with_indexes(filename, &[2], &[usize::MAX], true, true, f32::NAN);
    match result {
        Err(BedErrorPlus::BedError(BedError::SidIndexTooBig(_))) => (),
        _ => panic!("test failure"),
    };

    let mut ignore_val = nd::Array2::zeros((1, 1));
    let buf_reader = BufReader::new(File::open("bed_reader/tests/data/small_no_bim.bed").unwrap());
    let result = internal_read_no_alloc(
        buf_reader,
        &"ignore",
        usize::MAX,
        usize::MAX,
        true,
        &[usize::MAX - 1],
        &[usize::MAX - 1],
        f64::NAN,
        &mut ignore_val.view_mut(),
    );
    match result {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };

    let result = read(
        &"bed_reader/tests/data/no_such_file.nsf",
        true,
        true,
        f64::NAN,
    );
    match result {
        Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(_))) => (),
        _ => panic!("test failure"),
    };
}

// Python has more tests, but they aren't needed in Rust until it gets fancy indexing

#[test]
fn writer() {
    let filename = "bed_reader/tests/data/some_missing.bed";
    let val = read(filename, false, true, -127).unwrap();

    let temp = TempDir::default();
    let path2 = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_test.bed");
    let filename2 = path2.as_os_str().to_str().unwrap();

    write(filename2, &val.view(), true, -127).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }

    let val2 = read(filename2, false, true, -127).unwrap();
    assert!(allclose(&val.view(), &val2.view(), 0, true));

    let val = read(filename, false, true, f64::NAN).unwrap();

    let path2 = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_testf64.bed");
    let filename2 = path2.as_os_str().to_str().unwrap();

    write(filename2, &val.view(), true, f64::NAN).unwrap();
    for ext in ["fam", "bim"].iter() {
        let from = Path::new(filename).with_extension(ext);
        let to = Path::new(filename2).with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }

    let val2 = read(filename2, false, true, f64::NAN).unwrap();

    assert!(allclose(&val.view(), &val2.view(), 1e-8, true));

    let mut val = read(filename, false, true, f64::NAN).unwrap();
    val[(0, 0)] = 5.0;
    let path = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_testf64_5.bed");
    let filename = path.as_os_str().to_str().unwrap();
    let result = write(filename, &val.view(), true, f64::NAN);
    match result {
        Err(BedErrorPlus::BedError(BedError::BadValue(_))) => (),
        _ => panic!("test failure"),
    };

    let val = nd::Array2::zeros((0, 0));
    let path = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_testf64_0s.bed");
    let filename = path.as_os_str().to_str().unwrap();
    write(filename, &val.view(), true, f64::NAN).unwrap();

    let val: nd::Array2<i8> = nd::Array2::zeros((3, 0));
    let path = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_testf64_3_0.bed");
    let filename = path.as_os_str().to_str().unwrap();
    write(filename, &val.view(), true, -127).unwrap();
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

    matrix_subset_no_alloc(
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

    matrix_subset_no_alloc(
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

    let result = matrix_subset_no_alloc(&in_val2.view(), &[0], &[], &mut out_val2.view_mut());
    match result {
        Err(BedErrorPlus::BedError(BedError::SubsetMismatch(_, _, _, _))) => (),
        _ => panic!("test failure"),
    }
}

#[test]
fn fill_in() {
    let filename = "bed_reader/tests/data/some_missing.bed";

    for output_is_order_f_ptr in [false, true].iter() {
        let mut val = read(filename, *output_is_order_f_ptr, true, f64::NAN).unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));

        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Unit,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();
        assert!((val[(0, 0)] - 0.16783627165933704).abs() < 1e-8);

        nd::Array2::fill(&mut val, f64::NAN);
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Unit,
            true,
            false,
            &mut stats.view_mut(),
        );
        match result {
            Err(BedErrorPlus::BedError(BedError::NoIndividuals)) => (),
            _ => panic!("test failure"),
        }

        let mut val = read(filename, *output_is_order_f_ptr, true, f64::NAN).unwrap();
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Beta { a: -10.0, b: 0.0 },
            true,
            false,
            &mut stats.view_mut(),
        );
        match result {
            Err(BedErrorPlus::BedError(BedError::CannotCreateBetaDist(_, _))) => (),
            _ => panic!("test failure"),
        }

        nd::Array2::fill(&mut val, 3.0);
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Beta { a: 0.5, b: 0.5 },
            true,
            false,
            &mut stats.view_mut(),
        );
        match result {
            Err(BedErrorPlus::BedError(BedError::IllegalSnpMean)) => (),
            _ => panic!("test failure"),
        }

        nd::Array2::fill(&mut val, 1.0);
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Beta { a: 0.5, b: 0.5 },
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();
    }
}

#[test]
fn standardize_unit() {
    for output_is_order_f_ptr in [true, false].iter() {
        let mut val = read(
            r"bed_reader/tests/data/toydata.5chrom.bed",
            *output_is_order_f_ptr,
            false,
            f64::NAN,
        )
        .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Unit,
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
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };
    match try_div_4(0, 256, 0u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 10, 5u8) {
        Ok(tup) => assert_eq!(tup, (25usize, 25u8)),
        Err(_) => panic!("test failure"),
    };

    match try_div_4(25 * 4 + 1, 10, 5u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 11, 5u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };

    match try_div_4(25 * 4, 10, 6u8) {
        Err(BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn standardize_beta() {
    for output_is_order_f_ptr in [true, false].iter() {
        let mut val = read(
            r"bed_reader/tests/data/toydata.5chrom.bed",
            *output_is_order_f_ptr,
            false,
            f64::NAN,
        )
        .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            Dist::Beta { a: 1.0, b: 25.0 },
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
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed(_))) => (),
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
    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed(_))) => (),
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
    match result {
        Err(BedErrorPlus::BedError(BedError::BadMode(_))) => (),
        _ => panic!("test failure"),
    };
}

#[cfg(test)]
fn write_fake_metadata(filename: &str, iid_count: usize, sid_count: usize) {
    for (ext, count) in ["fam", "bim"].iter().zip([iid_count, sid_count].iter()) {
        let meta_name = std::path::Path::new(filename).with_extension(ext);
        let file = std::fs::File::create(meta_name).unwrap();
        let mut file = LineWriter::new(file);
        for _ in 0..*count {
            file.write_all(b"\n").unwrap();
        }
    }
}

#[test]
fn zeros() {
    let filename = "bed_reader/tests/data/some_missing.bed";
    let (iid_count, sid_count) = counts(filename).unwrap();
    let iid_index_full = (0..iid_count).collect::<Vec<usize>>();
    let sid_index_full = (0..sid_count).collect::<Vec<usize>>();
    let ref_val_float = reference_val(true);

    // Test read on zero length indexes
    let val = read(filename, true, true, f32::NAN).unwrap();
    assert!(allclose(&ref_val_float.view(), &val.view(), 1e-08, true));

    let out_val10 =
        read_with_indexes(filename, &iid_index_full, &[], true, true, f64::NAN).unwrap();
    assert!(out_val10.shape() == [iid_count, 0]);

    let out_val01 =
        read_with_indexes(filename, &[], &sid_index_full, true, true, f64::NAN).unwrap();
    assert!(out_val01.shape() == [0, sid_count]);

    let out_val00 = read_with_indexes(filename, &[], &[], true, true, f64::NAN).unwrap();
    assert!(out_val00.shape() == [0, 0]);

    // Test subset on zero length indexes

    let shape = (ref_val_float.shape()[0], ref_val_float.shape()[1], 1usize);
    let in_val = ref_val_float.into_shape(shape).unwrap();

    let mut out_val = nd::Array3::<f64>::zeros((iid_count, 0, 1));
    matrix_subset_no_alloc(
        &(in_val.view()),
        &iid_index_full,
        &[],
        &mut out_val.view_mut(),
    )
    .unwrap();

    let mut out_val = nd::Array3::<f64>::zeros((0, sid_count, 1));
    matrix_subset_no_alloc(
        &(in_val.view()),
        &[],
        &sid_index_full,
        &mut out_val.view_mut(),
    )
    .unwrap();

    let mut out_val = nd::Array3::<f64>::zeros((0, 0, 1));
    matrix_subset_no_alloc(&(in_val.view()), &[], &[], &mut out_val.view_mut()).unwrap();

    // Writing zero length vals
    let temp = TempDir::default();
    let path = PathBuf::from(temp.as_ref()).join("rust_bed_reader_writer_zeros.bed");
    let filename = path.as_os_str().to_str().unwrap();

    write(filename, &out_val01.view(), true, f64::NAN).unwrap();
    write_fake_metadata(filename, 0, sid_count);
    let result = read(filename, true, true, f64::NAN);
    let in_val01 = result.unwrap();
    assert!(in_val01.shape() == [0, sid_count]);
    assert!(allclose(&in_val01.view(), &out_val01.view(), 1e-08, true));

    write(filename, &out_val10.view(), true, f64::NAN).unwrap();
    write_fake_metadata(filename, iid_count, 0);
    let in_val10 = read(filename, true, true, f64::NAN).unwrap();
    assert!(in_val10.shape() == [iid_count, 0]);
    assert!(allclose(&in_val10.view(), &out_val10.view(), 1e-08, true));

    write(filename, &out_val00.view(), true, f64::NAN).unwrap();
    write_fake_metadata(filename, 0, 0);
    let in_val00 = read(filename, true, true, f64::NAN).unwrap();
    assert!(in_val00.shape() == [0, 0]);
    assert!(allclose(&in_val00.view(), &out_val00.view(), 1e-08, true));
}
#[test]
fn file_ata_small() {
    let filename = "bed_reader/tests/data/small_array.memmap";
    let mut out_val = nd::Array2::<f64>::from_elem((3, 3), f64::NAN);
    file_ata(filename, 0, 2, 3, 2, &mut out_val.view_mut()).unwrap();
    println!("{:?}", out_val);

    let expected = nd::arr2(&[[17., 22., 27.], [22., 29., 36.], [27., 36., 45.]]);
    println!("{:?}", expected);
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[cfg(test)]
fn file_ata(
    filename: &str,
    offset: u64,
    iid_count: usize,
    sid_count: usize,
    sid_step: usize,
    val: &mut nd::ArrayViewMut2<'_, f64>,
) -> Result<(), BedErrorPlus> {
    for sid_start in (0..sid_count).step_by(sid_step) {
        let sid_range_len = sid_step.min(sid_count - sid_start);
        let mut ata_piece =
            nd::Array2::<f64>::from_elem((sid_count - sid_start, sid_range_len), f64::NAN);
        file_ata_piece(
            filename,
            offset,
            iid_count,
            sid_count,
            sid_start,
            &mut ata_piece.view_mut(),
            sid_range_len,
            read_into_f64,
        )?;
        insert_piece(
            Range {
                start: sid_start,
                end: sid_start + sid_range_len,
            },
            ata_piece,
            val,
        );
    }
    return Ok(());
}

#[cfg(test)]
fn insert_piece(sid_range: Range<usize>, piece: nd::Array2<f64>, val: &mut nd::ArrayViewMut2<f64>) {
    for range_index in sid_range.clone() {
        for j in range_index - sid_range.start..piece.shape()[0] {
            // this is the inner loop, so pre-computing indexes would speed it up
            val[(range_index, j + sid_range.start)] = piece[(j, range_index - sid_range.start)];
            val[(j + sid_range.start, range_index)] = val[(range_index, j + sid_range.start)];
        }
    }
}

#[test]
fn file_b_less_aatbx_medium() {
    let filename = "bed_reader/tests/data/500x400_o640_array.memmap";
    let iid_count = 500usize;

    let b_shape = ShapeBuilder::set_f((iid_count, 100usize), true);
    let mut b1 = nd::Array2::<f64>::from_elem(b_shape, 2.0);
    let atb_shape = ShapeBuilder::set_f((400usize, 100usize), true);
    let mut atb = nd::Array2::<f64>::zeros(atb_shape);
    let mut aatb = b1.clone();

    file_b_less_aatbx(
        filename,
        640,
        iid_count,
        &mut b1.view_mut(),
        &mut aatb.view_mut(),
        &mut atb.view_mut(),
        10,
    )
    .unwrap();
    // println!("{:?}", atb);
    // println!("{:?}", aatb);

    println!("{:?}", atb[(1, 1)]);
    assert!(abs(atb[(1, 1)] - 499.00749503747534) < 1e-11);

    println!("{:?}", aatb[(1, 1)]);
    assert!(abs(aatb[(1, 1)] - -597.6363313483225) < 1e-11);
}

#[test]
fn file_aat_small() {
    let filename = "bed_reader/tests/data/small_array.memmap";
    let mut out_val = nd::Array2::<f64>::from_elem((2, 2), f64::NAN);
    file_aat(filename, 0, 2, 3, 1, &mut out_val.view_mut()).unwrap();
    println!("{:?}", out_val);

    let expected = nd::arr2(&[[14.0, 32.0], [32.0, 77.0]]);
    println!("{:?}", expected);
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[test]
fn file_aat_small2() {
    let filename = "bed_reader/tests/data/small2_array.memmap";
    let mut out_val = nd::Array2::<f64>::from_elem((3, 3), f64::NAN);
    file_aat(filename, 0, 3, 4, 2, &mut out_val.view_mut()).unwrap();
    println!("{:?}", out_val);

    let expected = nd::arr2(&[
        [30.0, 70.0, 110.0],
        [70.0, 174.0, 278.0],
        [110.0, 278.0, 446.0],
    ]);
    println!("{:?}", expected);
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[cfg(test)]
fn file_aat(
    filename: &str,
    offset: u64,
    iid_count: usize,
    sid_count: usize,
    iid_step: usize,
    val: &mut nd::ArrayViewMut2<'_, f64>,
) -> Result<(), BedErrorPlus> {
    let (nrows, ncols) = val.dim();
    assert!(nrows == iid_count && ncols == iid_count); // real assert
    for iid_start in (0..iid_count).step_by(iid_step) {
        let iid_range_len = iid_step.min(iid_count - iid_start);
        let mut aat_piece =
            nd::Array2::<f64>::from_elem((iid_count - iid_start, iid_range_len), f64::NAN);
        file_aat_piece(
            filename,
            offset,
            iid_count,
            sid_count,
            iid_start,
            &mut aat_piece.view_mut(),
            iid_range_len,
            read_into_f64,
        )?;
        println!("piece:\n{:?}", aat_piece);

        for range0_index in 0..iid_count - iid_start {
            for range1_index in 0..iid_range_len {
                let val00 = aat_piece[(range0_index, range1_index)];
                val[(range0_index + iid_start, range1_index + iid_start)] = val00;
                if range0_index > range1_index {
                    val[(range1_index + iid_start, range0_index + iid_start)] = val00;
                }
            }
        }
        println!("val:\n{:?}", val);
    }
    return Ok(());
}
