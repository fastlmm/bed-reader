from bed_reader import open_bed, sample_file

file_name = sample_file("sparse.bed")
with open_bed(file_name) as bed:
    val_sparse = bed.read_sparse(dtype="int8")
