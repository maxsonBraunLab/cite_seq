import os
import loompy
import glob

####################
# GLOBAL VARIABLES #
####################
input_files = snakemake.input.input_list #["file1.loom","file2.loom", ... ]
output_filename = snakemake.output[0]

#loom_files = [str(file) for file in input_files]

#input_dir = os.path.dirname(input_file)
#loom_files =  glob.glob("{}/*.loom".format(input_dir)) #ie ["file1.loom","file2.loom", ... ]
#loom_files =  input_files #ie ["file1.loom","file2.loom", ... ]

#combine loom files
loompy.combine(input_files, str(output_filename))