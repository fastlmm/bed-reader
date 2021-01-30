"""Read and write the PLINK BED format, simply and efficiently."""
__version__ = "0.1.3"

from bed_reader._open_bed import open_bed, get_num_threads
from bed_reader._sample_data import sample_file, tmp_path
from bed_reader._to_bed import to_bed
from .bed_reader import write_i8, write_f64, write_f32
from .bed_reader import standardize_f64, standardize_f32
from .bed_reader import subset_f64_f64, subset_f32_f64, subset_f32_f32

