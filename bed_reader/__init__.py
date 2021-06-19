"""Read and write the PLINK BED format, simply and efficiently."""

from bed_reader._open_bed import get_num_threads, open_bed  # noqa
from bed_reader._sample_data import sample_file, tmp_path  # noqa
from bed_reader._to_bed import to_bed  # noqa

from .bed_reader import file_ata_piece_f64  # !!!cmk make a f32, too?
from .bed_reader import (
    file_aat_piece_f64,
    file_b_less_aatbx,
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
