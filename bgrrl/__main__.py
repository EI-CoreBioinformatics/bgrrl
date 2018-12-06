import os
from os.path import join
import csv
import sys
import argparse
from argparse import ArgumentParser

from snakemake.utils import min_version
min_version("4.0")

from . import __version__, BGRRLRunner
from .snakemake_helper import make_exeenv_arg_group

# __version__ = "0.4.5"

BGSURVEY_DESC = "This is the quality control stage at which samples are tested for their assemble-ability."
BGASM_DESC = "This stage performs (short-read) assembly of samples passing the survey stage."
BGANN_DESC = "This stage performs de-novo (and on demand reference-based) gene/functional annotation."
BGPACKAGE_DESC = "This stage (re-)generates data packages on demand. By default, data packages are automatically generated after each stage."


def add_default_options(parser):

	common_group = parser.add_argument_group(
		"bgrr| options"
	)
	common_group.add_argument(                                                                                                                                   	
		"input_sheet", 
		help="""The samplesheet to process. This is a comma-separated file, containing location and meta-information for each sample to be processed."""
	)
	
	common_group.add_argument(
		"-o", 
		"--output-dir",
		default="Analysis",
		help="""Output will be written to the specified directory."""
	)

	common_group.add_argument(
		"--project-prefix",
		default="bgrrl_" + os.path.basename(os.getcwd()),
		help="""Reports and data packages will be prefixed by this."""
	)

	common_group.add_argument(
		"--config", 
		default="Analysis/config/bgrrl_config.yaml", #!!!
		help="""Configuration file. This file specifies details for accessing services and commands to be executed prior to running each pipeline tool."""
	) 	
	
	common_group.add_argument(
		"--report-only", 
		action="store_true", 
		help="""Only runs reporting modules, no snakemake pipelines [False]"""
	)

	common_group.add_argument(
		"-f", 
		"--force", 
		action="store_true", 
		help="""Force overwriting existing output directory, causes pipeline to be restarted. [DISABLED]"""
	)

	common_group.add_argument(
		"--enterobase-groups", 
		type=str, 
		help="""Comma-separated list of Enterobase microbial organisms. 
				The set of assemblies is tested against organism-specific criteria and assemblies are 
				packaged according to their species. [NEEDS REWORDING!]. 
				By default, the enterobase mode is disabled."""
	)

	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)


def add_survey_parser(subparsers):
	survey_parser = subparsers.add_parser(
		"survey",
		help=BGSURVEY_DESC,
		description=BGSURVEY_DESC
	)

	survey_parser.add_argument(
		"--no-normalization", 
		action="store_true", 
		help="""Disable read normalization. [False]""" # !TODO
	)

	survey_parser.add_argument(
		"--no-packaging",
		action="store_true",
		help="""Disable automatic packaging. [False]"""
	)

	survey_parser.add_argument(
		"--full-qaa-analysis",
		action="store_true",
		help="""Perform full qaa-analysis on survey assemblies. [False]"""
	)	

	add_default_options(survey_parser)
	survey_parser.set_defaults(runmode="survey")


def add_asm_parser(subparsers):
	asm_parser = subparsers.add_parser(
		"assemble",
		help=BGASM_DESC,
		description=BGASM_DESC
	)

	asm_parser.add_argument(
		"--assembler", 
		type=str, 
		choices=["unicycler", "velvet"],
		default="unicycler", 
		help="""Assembly software to use for genome assembly. [unicycler]""")

	asm_parser.add_argument(
		"--contig-minlen", 
		type=int, 
		default=0, 
		help="""Minimum length [bp] of contigs retained in filtering step [0]."""
	)

	asm_parser.add_argument(
        "--no-normalization",
        action="store_true",
        help="""Use non-normalized reads in asm module [False]"""
    )

	asm_parser.add_argument(
		"--no-packaging",
		action="store_true",
		help="""Disable automatic packaging. [False]"""
	)

	add_default_options(asm_parser)
	asm_parser.set_defaults(runmode="assemble")
	

def add_ann_parser(subparsers):
	ann_parser = subparsers.add_parser(
		"annotate", 
		help=BGANN_DESC,
		description=BGANN_DESC
	)

	ann_parser.add_argument(
		"--ratt-reference",
		type=str,
		help="Path to reference data for ratt annotation transfer"
	)

	ann_parser.add_argument(
		"--no-packaging",
		action="store_true",
		help="""Disable automatic packaging. [False]"""
	)

	ann_parser.add_argument(
		"--prokka-package-style",
		type=str,
		choices=["by_sample", "all_in_one"],
		default="by_sample",
		help="""Should the prokka annotation be packaged into one directory per sample (by_sample) or into one single directory (all_in_one)? [by_sample]"""
	)

	add_default_options(ann_parser)
	ann_parser.set_defaults(runmode="annotate")


def add_package_parser(subparsers):
	package_parser = subparsers.add_parser(
		"package", 
		help=BGPACKAGE_DESC,                  	
		description=BGPACKAGE_DESC
	)

	package_parser.add_argument(
		"package-mode", 
		# choices=["asm", "ann", "ann"],
		type=str, # comma-separated list ["asm", "ann"]
		help="""By default, packaging runs automatically after the assembly/annotation stages. 
				In case packages need to be (re-)generated, the packaging stage can be run for each module separately."""
	)

	package_parser.add_argument(
		"--prokka-package-style",
		type=str,
		choices=["by_sample", "all_in_one"],
		default="by_sample",
		help="""Should the prokka annotation be packaged into one directory per sample (by_sample) or into one single directory (all_in_one)? [by_sample]"""
	)

	package_parser.add_argument(
		"--make-ratt-data-tarballs",
		action="store_true",
		help="""ratt-based annotation transfers may generate large amounts of results. This option will generate a set of tarballs containing the transferred annotations. [False]"""
	)


	add_default_options(package_parser)
	package_parser.set_defaults(runmode="package")




def main():
	print("Starting EI bgrr| V " + __version__)
	print()

	if len(sys.argv) == 1:
		sys.argv.append("-h")

	bgrrl_parser = ArgumentParser(
		"bgrrl",
		description="""The Earlham Institute Bacterial Genome Reconstruction & Recognition Pipeline (bgrr|). 
					This utility is the central control program of bgrr|."""
	)

	stage_parsers = bgrrl_parser.add_subparsers(
		help="""bgrr| comprises three main stages: *survey*, *assemble*, and *annotate* and an auxiliary stage *package*.
				Reporting and quality assessment of the results of each main stage is performed using qaa."""
	)

	add_survey_parser(stage_parsers)
	add_asm_parser(stage_parsers)
	add_ann_parser(stage_parsers)
	add_package_parser(stage_parsers)
	

	args = bgrrl_parser.parse_args()
	print("ARGS", args)

	print(vars(args))

	BGRRLRunner(args).run()



if __name__ == "__main__":
	main()	
