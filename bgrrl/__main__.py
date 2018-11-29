#!/usr/bin/env python3

import os
from os.path import join
import csv
import sys
import shutil
import yaml
from argparse import ArgumentParser
from textwrap import dedent
import copy

from snakemake.utils import min_version
min_version("4.0")


from . import PipelineStep, __version__, BGRRLModuleRunner, BGRRLRunner

from .snakemake_helper import make_exeenv_arg_group, ExecutionEnvironment, NOW

VALID_ASSEMBLERS = ["unicycler", "velvet"]
VALID_ANNOTATORS = ["prokka", "ratt", "both"]

	
def main():
	print("Starting EI BGRRL V " + __version__)
	print()

	parser = ArgumentParser("The Earlham Institute Bacterial Genome Reconstruction & Recognition Pipeline (BGRR|)",
							description="""This program controls the various Snakemake pipelines making up the EI-BGRR| pipeline.""")

	parser.add_argument(
		"input", 
		help="""The (EI-PAP) samplesheet to process. This is a comma-separated file, 
				containing location and meta-information for each sample to be processed."""
	)

	parser.add_argument(
		"-o", "--output_dir", 
		default="BGRRL_<timestamp>",
		help="If specified BGRRL will output data to this directory."
	)

	parser.add_argument(
		"--project-prefix", 
		default=os.path.basename(os.getcwd()), 
		help="Reports and resultfiles/-folders will be prefixed by this."
	)

	parser.add_argument(
		"-f", "--force", 
		action="store_true", 
		help="Force overwriting existing output directory, causes pipeline to be restarted. (disabled)"
	)

	parser.add_argument(
		"--bgrrl_config", 
		help="""Configuration file for BGRRL. This file specifies details for accessing services and commands 
				to be executed prior to running each pipeline tool."""
	) 

	parser.add_argument(
		"-m", "--module", 
		choices=[ps.name.lower() for ps in PipelineStep], 
		default="read_qc",
		help="""This option controls which part of the pipeline to run. Due to the task at hand 
				being highly dependent on input data quality,
				it is recommended to run the 4 major steps of EI-BGRRL manually, carrying out manual/visual 
				data checks, and additional data manipulation steps in between steps.
				\"read_qc\" - This step preprocesses the read sets and assesses the assemblability 
				of individual samples. It will produce an input-samplesheet for the assembly-step.
				\"assembly\" - This step performs genome assembly and post-assembly quality control for 
				each sample.
				\"annotation\" - This step performs genome annotation and post-annotation quality control 
				for each assembly.
				\"report_and_package\" - This step performs report generation and data packaging.
				\"attempt_full\" - If you know the data and what you're doing, this option attempts to run the 
				individual pipeline steps consecutively without breaks (not implemented yet)"""
	)

	parser.add_argument(
		"--report-only", 
		action="store_true", 
		help="Only runs reporting modules, no snakemake pipelines [False]"
	)

	parser.add_argument(
		"--no-normalization", 
		action="store_true", 
		help="Use non-normalized reads in asm module (kills fallback-mode!) [False]"
	)

	parser.add_argument(
		"--package-mode", 
		choices=["ann", "asm"],
		help="""By default, packaging runs automatically after the assembly/annotation steps. 
				In case packages need to be (re-)generated, the finalize step can be ran for each module separately."""
	)

	parser.add_argument(
		"--no-packaging",
		action="store_true",
		help="""Disable automatic packaging. [False]"""

	parser.add_argument(
		"--enterobase-groups", 
		type=str, 
		help="""Comma-separated list of Enterobase microbial organisms. 
				The set of assemblies is tested against organism-specific criteria and assemblies are 
				packaged according to their species. [NEEDS REWORDING!]. 
				By default, the enterobase mode is disabled."""
	)

	parser.add_argument(
		"--assembler", 
		type=str, 
		choices=VALID_ASSEMBLERS,
		default="unicycler", 
		help="""Assembly software to use for genome assembly. [unicycler]"""
	)

	parser.add_argument(
		"--contig-minlen", 
		type=int, 
		default=0, 
		help="Minimum length [bp] of contigs retained in filtering step [0]."
	)

	parser.add_argument(
		"--annotation", 
		type=str, 
		choices=VALID_ANNOTATORS, 
		default="prokka",
		help="""Annotation software to use for genome annotation. [prokka]"""
	)

	parser.add_argument(
		"--prokka-package-style",
		type=str,
		choices=["by_sample", "all_in_one"],
		default="by_sample",
		help="""Should the prokka annotation be packaged into one directory per sample (by_sample) or into one single directory (all_in_one)? [by_sample]"""
	)

	parser.add_argument(
		"--ratt-reference-dir", 
		type=str, 
		help="Path to reference data for ratt", 
		default=""
	)

	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False)	# Add in cluster and DRMAA options
	args = parser.parse_args()

	bgrrl_runner = BGRRLRunner(args)
	run_result = bgrrl_runner.run()


if __name__ == "__main__":
	main()
