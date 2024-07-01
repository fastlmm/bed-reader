// https://stackoverflow.com/questions/32900809/how-to-suppress-function-is-never-used-warning-for-a-function-used-by-tests

#![cfg(test)]

use crate::allclose;
use crate::assert_eq_nan;
use crate::assert_error_variant;
use crate::encode1;
use crate::file_aat_piece;
use crate::file_ata_piece;
use crate::file_b_less_aatbx;
use crate::read_into_f64;
use crate::sample_bed_file;
use crate::sample_file;
use crate::sample_files;
use crate::try_div_4;
use crate::Bed;
use crate::Dist;
use crate::Index;
use crate::Metadata;
use crate::ReadOptions;
use crate::SliceInfo1;
use crate::WriteOptions;
use crate::{impute_and_zero_mean_snps, matrix_subset_no_alloc};
use crate::{internal_read_no_alloc, read_no_alloc, BedError, BedErrorPlus};
use anyinput::anyinput;
use nd::s;
use ndarray as nd;
use ndarray::ShapeBuilder;
use ndarray_npy::read_npy;
use num_traits::abs;
use std::f32;
use std::f64;
use std::io::BufReader;
use std::ops::Range;
use std::ops::RangeInclusive;
use std::path::Path;
use std::path::PathBuf;
use temp_testdir::TempDir;

#[test]
fn best_int8() {
    let filename = sample_bed_file("some_missing.bed").unwrap();

    for output_order_is_f in &[true, false] {
        let mut bed = Bed::new(&filename).unwrap();
        let val = ReadOptions::builder()
            .is_f(*output_order_is_f)
            .i8()
            .read(&mut bed)
            .unwrap();
        let ref_val_i8 = reference_val_i8(true);
        assert_eq!(val, ref_val_i8);
    }
}

fn reference_val_i8(is_a1_counted: bool) -> nd::Array2<i8> {
    let ref_val = reference_val(is_a1_counted);

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
    ref_val_i8
}

#[test]
fn read_test() {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed").unwrap();
    let mut bed = Bed::new(file).unwrap();
    assert!(bed.iid_count().unwrap() == 10);
    assert!(bed.sid_count().unwrap() == 100);
    let val: nd::Array2<i8> = bed.read().unwrap();
    let val_f64 = val.mapv(|elem| elem as f64);
    let mean_ = val_f64.mean().unwrap();
    assert!((mean_ - -13.142).abs() < 1e-8); // really shouldn't do mean on data where -127 represents missing

    let mut bed2 = Bed::new(sample_bed_file("small_too_short.bed").unwrap()).unwrap();
    let result = bed2.read::<i8>();
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));
}

#[test]
fn rest_reader_bed() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("some_missing.bed").unwrap();
    let is_a1_counted = false;

    let ref_val = reference_val(is_a1_counted);
    let ref_val_i8 = reference_val_i8(is_a1_counted);

    let mut bed = Bed::new(file).unwrap();
    let val = ReadOptions::builder()
        .is_a1_counted(is_a1_counted)
        .f32()
        .read(&mut bed)?;
    assert!(allclose(&ref_val.view(), &val.view(), 1e-08, true));

    let val_f64 = ReadOptions::builder()
        .is_a1_counted(is_a1_counted)
        .f64()
        .read(&mut bed)?;
    assert!(allclose(&ref_val.view(), &val_f64.view(), 1e-08, true));

    let val2 = ReadOptions::builder()
        .is_a1_counted(is_a1_counted)
        .i8()
        .read(&mut bed)?;
    assert_eq!(val2, ref_val_i8);

    Ok(())
}

fn reference_val(is_a1_counted: bool) -> nd::Array2<f64> {
    let file = sample_file("some_missing.val.npy").unwrap();

    let mut val: nd::Array2<f64> = read_npy(file).unwrap();
    if !is_a1_counted {
        val = val * -1.0 + 2.0;
    }

    val
}

#[test]
fn index() {
    let filename = sample_bed_file("some_missing.bed").unwrap();
    let ref_val_float = reference_val(true);

    let val: nd::Array2<f64> = Bed::new(&filename).unwrap().read().unwrap();
    assert!(allclose(&ref_val_float.view(), &val.view(), 1e-08, true));

    let mut bed = Bed::new(&filename).unwrap();
    let val: nd::Array2<f32> = ReadOptions::builder()
        .sid_index(2)
        .f32()
        .read(&mut bed)
        .unwrap();

    assert!(allclose(
        &(ref_val_float.slice(nd::s![.., 2usize..3])),
        &val.view(),
        1e-08,
        true
    ));

    let val = ReadOptions::builder()
        .iid_index(1)
        .sid_index(2)
        .f32()
        .read(&mut bed)
        .unwrap();
    assert!(allclose(
        &ref_val_float.slice(nd::s![1usize..2, 2usize..3]),
        &val.view(),
        1e-08,
        true
    ));

    // val = bed.read([2, -2])
    let val = ReadOptions::builder()
        .sid_index([2, -2])
        .f32()
        .read(&mut bed)
        .unwrap();

    let col0 = ref_val_float.slice(nd::s![.., 2]);
    let col1 = ref_val_float.slice(nd::s![.., -2]);
    let expected = nd::stack![nd::Axis(1), col0, col1];
    assert!(allclose(&expected.view(), &val.view(), 1e-08, true));

    let result = ReadOptions::builder()
        .iid_index(isize::MAX)
        .sid_index(2)
        .f32()
        .read(&mut bed);
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IidIndexTooBig(_)));

    let bed_bim = sample_files(["small_no_fam.bed", "small_no_fam.bim"]).unwrap();
    let mut bed = Bed::new(&bed_bim[0]).unwrap();
    let result2 = ReadOptions::builder()
        .iid_index(0)
        .sid_index(0)
        .f32()
        .read(&mut bed);
    assert_error_variant!(result2, BedErrorPlus::IOError(_));

    let bed_fam = sample_files(["small_no_bim.bed", "small_no_bim.fam"]).unwrap();
    let mut bed = Bed::new(&bed_fam[0]).unwrap();
    let result3 = ReadOptions::builder()
        .iid_index(0)
        .sid_index(0)
        .f32()
        .read(&mut bed);
    assert_error_variant!(result3, BedErrorPlus::IOError(_));

    let mut bed = Bed::new(filename).unwrap();
    let result4 = ReadOptions::builder()
        .iid_index(2)
        .sid_index(isize::MAX)
        .f32()
        .read(&mut bed);
    assert_error_variant!(result4, BedErrorPlus::BedError(BedError::SidIndexTooBig(_)));

    let mut ignore_val = nd::Array2::zeros((1, 1));
    let buf_reader = BufReader::new(std::fs::File::open(&bed_fam[0]).unwrap());
    let result5 = internal_read_no_alloc(
        buf_reader,
        "ignore",
        usize::MAX,
        usize::MAX,
        true,
        &[isize::MAX - 1],
        &[isize::MAX - 1],
        f64::NAN,
        &mut ignore_val.view_mut(),
    );
    assert_error_variant!(
        result5,
        BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))
    );

    let result6 = Bed::new("no_such_file.nsf");
    assert_error_variant!(result6, BedErrorPlus::IOError(_));
}

#[test]
fn writer() {
    let path = sample_bed_file("some_missing.bed").unwrap();

    let mut bed = Bed::new(&path).unwrap();
    let val = ReadOptions::builder().c().i8().read(&mut bed).unwrap();

    let output_folder = TempDir::default();
    let path2 = output_folder.join("rust_bed_reader_writer_test.bed");

    Bed::write(&val, &path2).unwrap();

    for ext in &["fam", "bim"] {
        let from = path.with_extension(ext);
        let to = path2.with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }

    let mut bed2 = Bed::new(path2).unwrap();
    let val2 = ReadOptions::builder().c().i8().read(&mut bed2).unwrap();

    assert!(allclose(&val.view(), &val2.view(), 0, true));

    let mut bed = Bed::new(&path).unwrap();
    let val = ReadOptions::builder().c().f64().read(&mut bed).unwrap();

    let path2 = output_folder.join("rust_bed_reader_writer_testf64.bed");

    Bed::write(&val, &path2).unwrap();

    for ext in &["fam", "bim"] {
        let from = path.with_extension(ext);
        let to = path2.with_extension(ext);
        std::fs::copy(from, to).unwrap();
    }

    let mut bed2 = Bed::new(path2).unwrap();
    let val2 = ReadOptions::builder().c().f64().read(&mut bed2).unwrap();
    assert!(allclose(&val.view(), &val2.view(), 1e-8, true));

    let mut val = ReadOptions::builder().c().f64().read(&mut bed).unwrap();
    val[(0, 0)] = 5.0;
    let path = output_folder.join("rust_bed_reader_writer_testf64_5.bed");

    let result = Bed::write(&val, &path);
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::BadValue(_)));
    assert!(!path.exists(), "file should not exist");

    // let val = nd::Array2::zeros((0, 0));
    let val = nd::Array2::<f64>::zeros((0, 0));
    let path = output_folder.join("rust_bed_reader_writer_testf64_0s.bed");
    Bed::write(&val, &path).unwrap();

    let val: nd::Array2<i8> = nd::Array2::zeros((3, 0));
    let path = output_folder.join("rust_bed_reader_writer_testf64_3_0.bed");
    Bed::write(&val, &path).unwrap();
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
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::SubsetMismatch(_, _, _, _))
    );
}

#[test]
fn fill_in() {
    let filename = sample_bed_file("some_missing.bed").unwrap();

    for output_is_orderf_ptr in &[false, true] {
        let mut bed = Bed::builder(&filename).build().unwrap();
        let mut val = ReadOptions::builder()
            .is_f(*output_is_orderf_ptr)
            .f64()
            .read(&mut bed)
            .unwrap();

        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));

        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Unit,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();
        assert!((val[(0, 0)] - 0.167_836_271_659_337_04).abs() < 1e-8);

        nd::Array2::fill(&mut val, f64::NAN);
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Unit,
            true,
            false,
            &mut stats.view_mut(),
        );
        assert_error_variant!(result, BedErrorPlus::BedError(BedError::NoIndividuals));

        let mut bed = Bed::builder(&filename).build().unwrap();
        let mut val = ReadOptions::builder()
            .is_f(*output_is_orderf_ptr)
            .f64()
            .read(&mut bed)
            .unwrap();
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Beta { a: -10.0, b: 0.0 },
            true,
            false,
            &mut stats.view_mut(),
        );
        assert_error_variant!(
            result,
            BedErrorPlus::BedError(BedError::CannotCreateBetaDist(_, _))
        );

        nd::Array2::fill(&mut val, 3.0);
        let result = impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Beta { a: 0.5, b: 0.5 },
            true,
            false,
            &mut stats.view_mut(),
        );
        assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllegalSnpMean));

        nd::Array2::fill(&mut val, 1.0);
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Beta { a: 0.5, b: 0.5 },
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();
    }
}

#[test]
fn standardize_unit() {
    for output_is_orderf_ptr in &[true, false] {
        let mut bed = Bed::new(sample_bed_file("toydata.5chrom.bed").unwrap()).unwrap();
        let mut val = ReadOptions::builder()
            .count_a2()
            .is_f(*output_is_orderf_ptr)
            .f64()
            .read(&mut bed)
            .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Unit,
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();

        assert!((val[(0, 0)] - -0.305_026_183_261_766_8).abs() < 1e-8);
    }
}

#[test]
fn div_4() {
    assert_error_variant!(
        try_div_4(usize::MAX, usize::MAX),
        BedErrorPlus::BedError(BedError::IndexesTooBigForFiles(_, _))
    );
}

#[test]
fn standardize_beta() {
    for output_is_orderf_ptr in &[true, false] {
        let mut bed = Bed::new(sample_bed_file("toydata.5chrom.bed").unwrap()).unwrap();
        let mut val = ReadOptions::builder()
            .count_a2()
            .is_f(*output_is_orderf_ptr)
            .f64()
            .read(&mut bed)
            .unwrap();
        let mut stats = nd::Array2::<f64>::zeros((val.dim().1, 2));
        impute_and_zero_mean_snps(
            &mut val.view_mut(),
            &Dist::Beta { a: 1.0, b: 25.0 },
            true,
            false,
            &mut stats.view_mut(),
        )
        .unwrap();

        assert!((val[(0, 0)] - -0.000_031_887_380_905_091_765).abs() < 1e-8);
    }
}

#[test]
fn read_errors() {
    let iid_count = 100usize;
    let sid_count = 200;
    let iid_index = (0..iid_count as isize).collect::<Vec<isize>>();
    let sid_index = (0..iid_count as isize).collect::<Vec<isize>>();
    let output_is_orderf = true;
    let shape = ShapeBuilder::set_f((iid_index.len(), sid_index.len()), output_is_orderf);
    let mut val = nd::Array2::<f64>::default(shape);

    let result0 = read_no_alloc(
        "no_such_file.nsf",
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        1,
        &mut val.view_mut(),
    );
    assert_error_variant!(result0, BedErrorPlus::IOError(_));

    let result = read_no_alloc(
        sample_file("some_missing.fam").unwrap(),
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        1,
        &mut val.view_mut(),
    );
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));

    let result = read_no_alloc(
        sample_file("empty.bed").unwrap(),
        iid_count,
        sid_count,
        true,
        &iid_index,
        &sid_index,
        f64::NAN,
        1,
        &mut val.view_mut(),
    );
    assert_error_variant!(result, BedErrorPlus::IOError(_));
}

#[test]
fn read_modes() -> Result<(), Box<BedErrorPlus>> {
    let filename = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(filename)?;
    let iid_count_s1 = bed.iid_count()?;
    let sid_count_s1 = bed.sid_count()?;

    let mut val_small_mode_1 = nd::Array2::<i8>::default((iid_count_s1, sid_count_s1));
    bed.read_and_fill(&mut val_small_mode_1.view_mut())?;

    let bed_fam_bim = sample_files(["small_too_short.bed", "small.fam", "small.bim"])?;
    let mut bed_too_short = Bed::builder(&bed_fam_bim[0])
        .fam_path(&bed_fam_bim[1])
        .bim_path(&bed_fam_bim[2])
        .build()?;
    let result = bed_too_short.read_and_fill(&mut val_small_mode_1.view_mut());
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));

    let mut val_small_mode_0 = nd::Array2::<i8>::default((sid_count_s1, iid_count_s1));
    let mut bed_mode0 = Bed::new(sample_bed_file("smallmode0.bed")?)?;
    bed_mode0.read_and_fill(&mut val_small_mode_0.view_mut())?;
    assert_eq!(val_small_mode_0.t(), val_small_mode_1);

    let bed_fam_bim = sample_files(["smallmodebad.bed", "small.fam", "small.bim"])?;
    let mut bed_small_mode_bad = Bed::builder(&bed_fam_bim[0])
        .fam_path(&bed_fam_bim[1])
        .bim_path(&bed_fam_bim[2])
        .build()?;
    let result = bed_small_mode_bad.read_and_fill(&mut val_small_mode_1.view_mut());
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::BadMode(_)));

    Ok(())
}

#[test]
fn zeros() -> Result<(), Box<BedErrorPlus>> {
    let filename = sample_bed_file("some_missing.bed")?;
    let mut bed = Bed::new(&filename).unwrap();
    let iid_count = bed.iid_count().unwrap();
    let sid_count = bed.sid_count().unwrap();
    let iid_index_full = (0..iid_count).collect::<Vec<usize>>();
    let sid_index_full = (0..sid_count).collect::<Vec<usize>>();
    let ref_val_float = reference_val(true);

    // Test read on zero length indexes
    let mut bed = Bed::new(&filename).unwrap();
    let val: nd::Array2<f32> = bed.read().unwrap();
    assert!(allclose(&ref_val_float.view(), &val.view(), 1e-08, true));

    let out_val10 = ReadOptions::builder()
        .sid_index([0; 0])
        .f64()
        .read(&mut bed)
        .unwrap();
    assert!(out_val10.dim() == (iid_count, 0));

    let out_val01 = ReadOptions::builder()
        .iid_index([0; 0])
        .f64()
        .read(&mut bed)
        .unwrap();
    assert!(out_val01.dim() == (0, sid_count));

    let out_val00 = ReadOptions::builder()
        .iid_index([0; 0])
        .sid_index([0; 0])
        .f64()
        .read(&mut bed)
        .unwrap();
    assert!(out_val00.dim() == (0, 0));

    // Test subset on zero length indexes

    let shape = (ref_val_float.dim().0, ref_val_float.dim().1, 1usize);
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
    let output_folder = TempDir::default();
    let path = output_folder.join("rust_bed_reader_writer_zeros.bed");

    Bed::write(&out_val01, &path).unwrap();
    let in_val01 = Bed::new(&path).unwrap().read::<f64>().unwrap();
    assert!(in_val01.dim() == (0, sid_count));
    assert!(allclose(&in_val01.view(), &out_val01.view(), 1e-08, true));

    Bed::write(&out_val10, &path).unwrap();
    let in_val10 = Bed::new(&path).unwrap().read::<f64>().unwrap();
    assert!(in_val10.dim() == (iid_count, 0));
    assert!(allclose(&in_val10.view(), &out_val10.view(), 1e-08, true));

    Bed::write(&out_val00, &path).unwrap();
    let in_val00 = Bed::new(&path).unwrap().read::<f64>().unwrap();
    assert!(in_val00.dim() == (0, 0));
    assert!(allclose(&in_val00.view(), &out_val00.view(), 1e-08, true));

    Ok(())
}
#[test]
fn file_ata_small() {
    let filename = sample_file("small_array.memmap").unwrap();
    let mut out_val = nd::Array2::<f64>::from_elem((3, 3), f64::NAN);
    file_ata(filename, 0, 2, 3, 2, &mut out_val.view_mut()).unwrap();
    println!("{out_val:?}");

    let expected = nd::arr2(&[[17., 22., 27.], [22., 29., 36.], [27., 36., 45.]]);
    println!("{expected:?}");
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[anyinput]
fn file_ata(
    path: AnyPath,
    offset: u64,
    iid_count: usize,
    sid_count: usize,
    sid_step: usize,
    val: &mut nd::ArrayViewMut2<'_, f64>,
) -> Result<(), Box<BedErrorPlus>> {
    for sid_start in (0..sid_count).step_by(sid_step) {
        let sid_range_len = sid_step.min(sid_count - sid_start);
        let mut ata_piece =
            nd::Array2::<f64>::from_elem((sid_count - sid_start, sid_range_len), f64::NAN);
        file_ata_piece(
            path,
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
            &ata_piece,
            val,
        );
    }
    Ok(())
}

fn insert_piece(
    sid_range: Range<usize>,
    piece: &nd::Array2<f64>,
    val: &mut nd::ArrayViewMut2<f64>,
) {
    for range_index in sid_range.clone() {
        for j in range_index - sid_range.start..piece.dim().0 {
            // this is the inner loop, so pre-computing indexes would speed it up
            val[(range_index, j + sid_range.start)] = piece[(j, range_index - sid_range.start)];
            val[(j + sid_range.start, range_index)] = val[(range_index, j + sid_range.start)];
        }
    }
}

#[test]
fn file_b_less_aatbx_medium() {
    let filename = sample_file("500x400_o640_array.memmap").unwrap();
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

    println!("{:?}", atb[(1, 1)]);
    assert!(abs(atb[(1, 1)] - 499.007_495_037_475_34) < 1e-11);

    println!("{:?}", aatb[(1, 1)]);
    assert!(abs(aatb[(1, 1)] - -597.636_331_348_322_5) < 1e-11);
}

#[test]
fn file_aat_small() {
    let filename = sample_file("small_array.memmap").unwrap();
    let mut out_val = nd::Array2::<f64>::from_elem((2, 2), f64::NAN);
    file_aat(filename, 0, 2, 3, 1, &mut out_val.view_mut()).unwrap();
    println!("{out_val:?}");

    let expected = nd::arr2(&[[14.0, 32.0], [32.0, 77.0]]);
    println!("{expected:?}");
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[test]
fn file_aat_small2() {
    let filename = sample_file("small2_array.memmap").unwrap();
    let mut out_val = nd::Array2::<f64>::from_elem((3, 3), f64::NAN);
    file_aat(filename, 0, 3, 4, 2, &mut out_val.view_mut()).unwrap();
    println!("{out_val:?}");

    let expected = nd::arr2(&[
        [30.0, 70.0, 110.0],
        [70.0, 174.0, 278.0],
        [110.0, 278.0, 446.0],
    ]);
    println!("{expected:?}");
    assert!(allclose(&expected.view(), &out_val.view(), 1e-08, true));
}

#[anyinput]
fn file_aat(
    path: AnyPath,
    offset: u64,
    iid_count: usize,
    sid_count: usize,
    iid_step: usize,
    val: &mut nd::ArrayViewMut2<'_, f64>,
) -> Result<(), Box<BedErrorPlus>> {
    let (nrows, ncols) = val.dim();
    assert!(nrows == iid_count && ncols == iid_count); // real assert
    for iid_start in (0..iid_count).step_by(iid_step) {
        let iid_range_len = iid_step.min(iid_count - iid_start);
        let mut aat_piece =
            nd::Array2::<f64>::from_elem((iid_count - iid_start, iid_range_len), f64::NAN);
        file_aat_piece(
            path,
            offset,
            iid_count,
            sid_count,
            iid_start,
            &mut aat_piece.view_mut(),
            iid_range_len,
            read_into_f64,
        )?;
        println!("piece:\n{aat_piece:?}");

        for range0_index in 0..iid_count - iid_start {
            for range1_index in 0..iid_range_len {
                let val00 = aat_piece[(range0_index, range1_index)];
                val[(range0_index + iid_start, range1_index + iid_start)] = val00;
                if range0_index > range1_index {
                    val[(range1_index + iid_start, range0_index + iid_start)] = val00;
                }
            }
        }
        println!("val:\n{val:?}");
    }
    Ok(())
}

#[test]
fn test_allclose() -> Result<(), Box<BedErrorPlus>> {
    let val1 = nd::arr2(&[[1.0, 2.000_000_000_001], [3.0, f64::NAN]]);
    let val2 = nd::arr2(&[[1.0, 2.0], [3.0, f64::NAN]]);
    assert!(allclose(&val1.view(), &val2.view(), 1e-08, true));

    let val1 = nd::arr2(&[[1.0, 2.0], [3.0, f64::NAN]]);
    assert_eq_nan(&val1, &val2);

    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.bed");
    let val = nd::array![
        [1.0, 0.0, f64::NAN, 0.0],
        [2.0, 0.0, f64::NAN, 2.0],
        [0.0, 1.0, 2.0, 0.0]
    ];
    WriteOptions::builder(output_file).write(&val)?;

    Ok(())
}

fn expected_len(index: &Index, count: usize, len: usize) -> Result<(), Box<BedErrorPlus>> {
    assert!(index.to_vec(count)?.len() == len);
    assert!(index.len(count)? == len);
    assert!(index.is_empty(count)? == (len == 0));

    Ok(())
}

#[test]
fn index_len_is_empty() -> Result<(), Box<BedErrorPlus>> {
    expected_len(&s![0..0;-2].into(), 0, 0)?;
    expected_len(&s![0..;-2].into(), 4, 2)?;

    expected_len(&Index::All, 0, 0)?;
    expected_len(&Index::All, 2, 2)?;

    expected_len(&(-1).into(), 2, 1)?;

    expected_len(&(vec![] as Vec<isize>).into(), 0, 0)?;
    expected_len(&vec![2, -1].into(), 4, 2)?;

    expected_len(&(nd::array![] as nd::Array1<isize>).into(), 0, 0)?;
    expected_len(&nd::array![2, -1].into(), 4, 2)?;

    let empty_isize = nd::array![] as nd::Array1<isize>;
    expected_len(&(empty_isize.view()).into(), 0, 0)?;
    expected_len(&(&(empty_isize.view())).into(), 0, 0)?;
    expected_len(&(&nd::array![2, -1].view()).into(), 4, 2)?;
    expected_len(&(nd::array![2, -1].view()).into(), 4, 2)?;

    expected_len(&(vec![] as Vec<bool>).into(), 0, 0)?;
    expected_len(&vec![false, false, true, true].into(), 4, 2)?;

    let empty_bool = nd::array![] as nd::Array1<bool>;
    expected_len(&(empty_bool.view()).into(), 0, 0)?;
    expected_len(&(&empty_bool.view()).into(), 0, 0)?;
    expected_len(&(nd::array![] as nd::Array1<bool>).into(), 0, 0)?;
    expected_len(&nd::array![false, false, true, true].into(), 4, 2)?;

    expected_len(&(0..).into(), 0, 0)?;
    expected_len(&(0..).into(), 2, 2)?;
    Ok(())
}

#[test]
fn test_sample_file() -> Result<(), Box<BedErrorPlus>> {
    let filename = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(filename)?;
    println!("{}", bed.iid_count()?);
    println!("{}", bed.sid_count()?);

    let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    let mut bed = Bed::builder(&deb_maf_mib[0])
        .fam_path(&deb_maf_mib[1])
        .bim_path(&deb_maf_mib[2])
        .build()?;
    println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    println!("{:?}", bed.sid()?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]

    Ok(())
}

#[test]
#[allow(clippy::needless_borrow)]
#[allow(clippy::needless_borrows_for_generic_args)]
fn demo_path() -> Result<(), Box<BedErrorPlus>> {
    let path: &str = "bed_reader/tests/data/small.bed";
    let _ = Bed::new(&path)?; // borrow a &str
    let _ = Bed::new(path)?; // move a &str
    let path: String = "bed_reader/tests/data/small.bed".to_string();
    let _ = Bed::new(&path)?; // borrow a String
    let path2: &String = &path;
    let _ = Bed::new(&path2)?; // borrow a &String
    let _ = Bed::new(path2)?; // move a &String
    let _ = Bed::new(path)?; // move a String
    let path: &Path = Path::new("bed_reader/tests/data/small.bed");
    let _ = Bed::new(&path)?; // borrow a Path
    let _ = Bed::new(path)?; // move a Path
    let path: PathBuf = PathBuf::from("bed_reader/tests/data/small.bed");
    let _ = Bed::new(&path)?; // borrow a PathBuf
    let path2: &PathBuf = &path;
    let _ = Bed::new(&path2)?; // borrow a &PathBuf
    let _ = Bed::new(path2)?; // move a &PathBuf
    let _ = Bed::new(path)?; // move a PathBuf
    Ok(())
}

#[allow(clippy::single_char_pattern)]
#[allow(clippy::needless_borrow)]
#[allow(clippy::needless_borrows_for_generic_args)]
#[test]
fn demo_iter() -> Result<(), Box<BedErrorPlus>> {
    let list: [&str; 3] = ["i1", "i2", "i3"];
    let _ = Metadata::builder().iid(&list).build()?; // borrow fixed-size array
    let _ = Metadata::builder().iid(list).build()?; // move fixed-size array
    let list: [String; 3] = ["i1".to_string(), "i2".to_string(), "i3".to_string()];
    let _ = Metadata::builder().iid(&list).build()?; // borrow fixed-size array of String
    let _ = Metadata::builder().iid(list).build()?; // move fixed-size array of String
    let list: Vec<&str> = vec!["i1", "i2", "i3"];
    let _ = Metadata::builder().iid(&list).build()?; // borrow Vec<&str>
    let list2 = &list[..]; // borrowed slice
    let _ = Metadata::builder().iid(list2).build()?; // borrow slice
    let _ = Metadata::builder().iid(list).build()?; // move Vec<&str>
    let list = nd::array!["i1", "i2", "i3"];
    let view = list.view();
    let _ = Metadata::builder().iid(&view).build()?; // borrow nd view
    let _ = Metadata::builder().iid(view).build()?; // move nd view
    let _ = Metadata::builder().iid(&list).build()?; // borrow ndarray
    let _ = Metadata::builder().iid(list).build()?; // move ndarray
    let list: std::str::Split<&str> = "i1,i2,i3".split(",");
    let _ = Metadata::builder().iid(list).build()?; // move iterator
    Ok(())
}

#[test]
#[allow(clippy::needless_borrow)]
#[allow(clippy::needless_borrows_for_generic_args)]
fn demo_index() -> Result<(), Box<BedErrorPlus>> {
    #[allow(clippy::let_unit_value)]
    let index: () = ();
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: isize = 2;
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index = ..;
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: RangeInclusive<usize> = 0..=3;
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: SliceInfo1 = s![..;2];
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: [isize; 2] = [2, 5];
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: Vec<isize> = vec![2, 5];
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: &[isize] = &vec![2, 5][..];
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: nd::Array1<isize> = nd::array![2, 5];
    let view = index.view();
    let _ = ReadOptions::builder().iid_index(&view).i8().build()?;
    let _ = ReadOptions::builder().iid_index(view).i8().build()?;
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: [bool; 3] = [false, false, true];
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: Vec<bool> = vec![false, false, true];
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: &[bool] = &vec![false, false, true][..];
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    let index: nd::Array1<bool> = nd::array![false, false, true];
    let view = index.view();
    let _ = ReadOptions::builder().iid_index(&view).i8().build()?;
    let _ = ReadOptions::builder().iid_index(view).i8().build()?;
    let _ = ReadOptions::builder().iid_index(&index).i8().build()?;
    let _ = ReadOptions::builder().iid_index(index).i8().build()?;

    Ok(())
}

#[test]
fn demo_index2() -> Result<(), Box<BedErrorPlus>> {
    let _ = ReadOptions::builder().iid_index(()).i8().build()?;
    let _ = ReadOptions::builder().iid_index(2).i8().build()?;
    let _ = ReadOptions::builder().iid_index(..).i8().build()?;
    let _ = ReadOptions::builder().iid_index(0..=3).i8().build()?;
    let _ = ReadOptions::builder().iid_index(s![..;2]).i8().build()?;
    let _ = ReadOptions::builder().iid_index([2, 5]).i8().build()?;
    let _ = ReadOptions::builder().iid_index(vec![2, 5]).i8().build()?;
    let _ = ReadOptions::builder()
        .iid_index(&vec![2, 5][..])
        .i8()
        .build()?;
    let _ = ReadOptions::builder()
        .iid_index(nd::array![2, 5])
        .i8()
        .build()?;
    let _ = ReadOptions::builder()
        .iid_index([false, false, true])
        .i8()
        .build()?;
    let _ = ReadOptions::builder()
        .iid_index(vec![false, false, true])
        .i8()
        .build()?;
    let _ = ReadOptions::builder()
        .iid_index(nd::array![false, false, true])
        .i8()
        .build()?;
    Ok(())
}

#[test]
#[allow(clippy::needless_borrow)]
#[allow(clippy::needless_borrows_for_generic_args)]
fn use_index() -> Result<(), Box<BedErrorPlus>> {
    fn len100(index: impl Into<Index>) -> Result<usize, Box<BedErrorPlus>> {
        let index = index.into();
        let len = index.len(100)?;
        Ok(len)
    }

    #[allow(clippy::let_unit_value)]
    let index: () = ();
    let _ = len100(index)?;

    let index = 2;
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = ..;
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = 0..=3;
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = s![..;2];
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = [2, 5];
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = vec![2, 5];
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = &vec![2, 5][..];
    let _ = len100(index)?;

    let index = nd::array![2, 5];
    let view = index.view();
    let _ = len100(&view)?;
    let _ = len100(view)?;
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = [false, false, true];
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = vec![false, false, true];
    let _ = len100(&index)?;
    let _ = len100(index)?;

    let index = &vec![false, false, true][..];
    let _ = len100(index)?;

    let index = nd::array![false, false, true];
    let view = index.view();
    let _ = len100(&view)?;
    let _ = len100(view)?;
    let _ = len100(&index)?;
    let _ = len100(index)?;

    Ok(())
}

#[test]
fn another_bed_read_example() -> Result<(), Box<BedErrorPlus>> {
    let filename = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(filename)?;
    let val = ReadOptions::builder()
        .sid_index(..3)
        .i8()
        .c()
        .num_threads(1)
        .read(&mut bed)?;
    println!("{:?}", val.dim());
    Ok(())
}

#[test]
fn read_encode1() -> Result<(), Box<BedErrorPlus>> {
    let val = nd::Array1::from(vec![2i8, 2, 2, 2, -127, -127, 2, 2, 1, -127, -127, -127]);
    let mut buffer = vec![0u8; 3];
    encode1(&val.view(), &mut buffer, true, -127)?;
    assert_eq!(buffer, vec![0, 5, 86]);
    Ok(())
}

#[test]
fn zero_length_encode1() -> Result<(), Box<BedErrorPlus>> {
    let val = nd::Array1::from(vec![0i8; 0]); // zero length
    let mut buffer = vec![0u8; 0]; // zero length
    encode1(&val.view(), &mut buffer, true, -127)?;
    assert_eq!(buffer, vec![0u8; 0]); // zero length
    Ok(())
}
