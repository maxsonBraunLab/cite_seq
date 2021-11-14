""" Process CITE-Seq data """

""" 1 - preprocess will plot important QC metrics for data filtration """
""" 2 - integrate will find anchors and integrate all samples into an integrated seurat object """
""" 3 - cluster will cluster individual samples plus integrated object """
""" 4 - differential will test for DE genes in seurat objects """

import os
import sys
import time
import glob

configfile: "config.yaml"

# SAMPLES, = glob_wildcards("data/raw/{sample}")
SAMPLES = [ os.path.basename(i) for i in glob.glob("data/raw/*")]

def message(m):
	sys.stderr.write("|--- {} \n".format(m))

for i in SAMPLES:
	message("Samples to process: {}".format(i))

if config["date"]:
	time_stamp = config["date"]
else:
	time_stamp = time.strftime("%Y-%m-%d")

message("time_stamp: {}".format(time_stamp))

p = config["project_name"]

if not os.path.isdir("data/reports"):
	os.mkdir("data/reports")
if not os.path.isdir("data/rda"):
	os.mkdir("data/rda")

# snakemake -j 8 --use-conda --cluster-config cluster.yaml --profile slurm

rule all:
	input:
		# 1 - preprocess
		"data/reports/1-preprocess.html",
		"data/rda/preprocess.{project_name}.{date}.rds".format(project_name = p, date = time_stamp),
		# 2 - integrate
		"data/reports/2-integrate.html",
		"data/rda/integrate.{project_name}.{date}.rds".format(project_name = p, date = time_stamp),
		# 3 - cluster
		"data/reports/3-cluster.html",
		"data/rda/cluster.{project_name}.{date}.rds".format(project_name = p, date = time_stamp),
		# 4 - differential
		"data/reports/4-differential.html"

rule preprocess:
	input:
		rmd = "scripts/1-preprocess.Rmd",
		samples = expand("data/raw/{sample}/outs/filtered_feature_bc_matrix", sample = SAMPLES)
	output:
		report = "data/reports/1-preprocess.html",
		rds = "data/rda/preprocess.{project_name}.{date}.rds".format(project_name = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	log:
		"logs/1-preprocess.log"
	shell:
		"""
		Rscript -e 'rmarkdown::render( here::here("{input.rmd}"), output_file = here::here("{output.report}"), knit_root_dir = here::here(), envir = new.env(), params = list(input_samples = "{input.samples}", output_rds = "{output.rds}" ))' > {log} 2>&1
		"""

# knit rmarkdown report with multiple outputs workaround: https://github.com/snakemake/snakemake/issues/178
# how to use rmd params to specify input/output files:
# https://stackoverflow.com/questions/32479130/passing-parameters-to-r-markdown
# here::here() will turn a relative path absolute, and is required for rmarkdown I/O

rule integrate:
	input:
		rmd = "scripts/2-integrate.Rmd",
		rds = rules.preprocess.output.rds
	output:
		report = "data/reports/2-integrate.html",
		rds = "data/rda/integrate.{project_name}.{date}.rds".format(project_name = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	log:
		"logs/2-integrate.log"
	shell:
		"""
		Rscript -e 'rmarkdown::render( here::here("{input.rmd}"), output_file = here::here("{output.report}"), knit_root_dir = here::here(), envir = new.env(), params = list(input_rds = "{input.rds}", output_rds = "{output.rds}" ))' > {log} 2>&1
		"""

rule cluster:
	input:
		rmd = "scripts/3-cluster.Rmd",
		rds = rules.integrate.output.rds
	output:
		report = "data/reports/3-cluster.html",
		rds = "data/rda/cluster.{project_name}.{date}.rds".format(project_name = p, date = time_stamp),
		umaps = directory("data/umaps")
	conda:
		"envs/seurat.yaml"
	log:
		"logs/3-cluster.log"
	shell:
		"""
		Rscript -e 'rmarkdown::render( here::here("{input.rmd}"), output_file = here::here("{output.report}"), knit_root_dir = here::here(), envir = new.env(), params = list(input_rds = "{input.rds}", output_rds = "{output.rds}" ))' > {log} 2>&1
		"""

rule differential:
	input:
		rmd = "scripts/4-differential.Rmd",
		rds = rules.cluster.output.rds
	output:
		# directory("data/markers"),
		report = "data/reports/4-differential.html"
	conda:
		"envs/seurat_go.yaml"
	log:
		"logs/4-differential.log"
	shell:
		"""
		Rscript -e 'rmarkdown::render( here::here("{input.rmd}"), output_file = here::here("{output.report}"), knit_root_dir = here::here(), envir = new.env(), params = list(input_rds = "{input.rds}" ))' > {log} 2>&1
		"""