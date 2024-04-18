"""Read and write the PLINK BED format, simply and efficiently."""

# flake8: noqa: F401

# See: https://towardsdatascience.com/nine-rules-for-writing-python-extensions-in-rust-d35ea3a4ec29?sk=f8d808d5f414154fdb811e4137011437
# for an article on how this project uses Rust to create a Python extension.

from bed_reader._open_bed import get_num_threads, open_bed
from bed_reader._sample_data import sample_file, sample_url, tmp_path
from bed_reader._to_bed import create_bed, to_bed

from .bed_reader import (
    encode1_i8,
    file_aat_piece_f32_orderf,
    file_aat_piece_f64_orderf,
    file_ata_piece_f32_orderf,
    file_ata_piece_f64_orderf,
    file_b_less_aatbx,
    file_dot_piece,
    read_cloud_i8,
    read_f32,
    read_f64,
    read_i8,
    standardize_f32,
    standardize_f64,
    subset_f32_f32,
    subset_f32_f64,
    subset_f64_f64,
    write_f32,
    write_f64,
    write_i8,
)
