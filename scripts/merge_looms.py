import os
import loompy

####################
# GLOBAL VARIABLES #
####################
input_files = snakemake.input.input_list #["file1.loom","file2.loom", ... ]
output_filename = snakemake.output[0]

#combine loom files
loompy.combine(files = input_files, output_file = str(output_filename), key = )