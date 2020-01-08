import os
from os.path import join
import csv
import sys
import argparse
from argparse import ArgumentParser

from snakemake.utils import min_version
min_version("4.0")

from . import __version__
from bgrrl.runners import BgrrlRunner
from eicore.snakemake_helper import make_exeenv_arg_group

BGSURVEY_DESC = "This is the quality control stage at which samples are tested for their assemble-ability."
BGASM_DESC = "This stage performs (short-read) assembly of samples passing the survey stage."
BGANN_DESC = "This stage performs de-novo (and on demand reference-based) gene/functional annotation."
BGPACKAGE_DESC = "This stage (re-)generates data packages on demand. By default, data packages are automatically generated after each stage."
BGAUTO_DESC = "This mode runs the full pipeline (survey + assembly, with optional annotation and packaging) in one go."


def add_default_options(parser):

	common_group = parser.add_argument_group(
		"bgrr| options"
	)
	
	common_group.add_argument(
		"input_sheet", 
		help="""The samplesheet to process. This is a comma-separated file, 
				containing location and meta-information for each sample to be processed."""
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
		"--no-normalization",                                 	
		action="store_true", 
		help="""Disable read normalization. [False]""" # !TODO
	)

	common_group.add_argument(
		"--enterobase-groups", 
		type=str, 
		help="""Comma-separated list of Enterobase microbial organisms. 
				The set of assemblies is tested against organism-specific criteria and assemblies are 
				packaged according to their species. [NEEDS REWORDING!]. 
				By default, the enterobase mode is disabled."""
	)

	common_group.add_argument(
		"--single-cell",
		action="store_true",
		help="""Support for single-cell data"""
	)



def add_survey_group(parser, shared_only=False):
	
	group = parser.add_argument_group(
        "survey options"
    )
	
	group.add_argument(
		"--full-qaa-analysis",
		action="store_true",
		help="""Perform full qaa-analysis on survey assemblies. [False]"""
	)	

	group.add_argument(
		"--minimum-survey-assembly-size",
		type=int,
		default=1000000,
		help="""Minimum size (in bp) for tadpole assembly to pass survey stage [1Mbp] This allows plasmids to be processed."""
	)


def add_asm_group(parser, shared_only=False):

	group = parser.add_argument_group(
        "assembly options"
    )

	group.add_argument(
		"--assembler", 
		type=str, 
		choices=["unicycler", "velvet"],
		default="unicycler", 
		help="""Assembly software to use for genome assembly. [unicycler]""")

	group.add_argument(
		"--contig-minlen", 
		type=int, 
		default=0, 
		help="""Minimum length [bp] of contigs retained in filtering step [0]."""
	)

	group.add_argument(
		"--run-annotation",
		action="store_true",
		help="""Run annotation on assembly. If set, de novo annotation with prokka will be run. 
				Additionally, you may enable annotation transfer by specifying a path to a reference annotation with --ratt-reference. [False]"""
	)

	group.add_argument(
		"--is-final-step",
		action="store_true",
		help="""If set, analysis packaging will take place after the assembly stage.
				Otherwise, assume that an annotation stage will follow, which will then take care of analysis packaging. [False]"""
	)

def add_ann_group(parser, shared_only=False):

	group = parser.add_argument_group(
        "annotation options"
    )

	group.add_argument(
		"--custom-prokka-proteins",
		type=str,
		default="",
		help="""If you have a custom protein database that you would like prokka to use (--proteins option), then specify the path to it here. [n/a]"""
	)

	group.add_argument(
		"--ratt-reference",
		type=str,
		help="Path to reference data for ratt annotation transfer"
	)


def add_package_group(parser, is_packaging_module=False):

	group = parser.add_argument_group(
        "packaging options"
    )

	if is_packaging_module:
		group.add_argument(
			"package-mode", 
			type=str, 
			help="""By default, packaging runs automatically after the assembly/annotation stages. 
					In case packages need to be (re-)generated, the packaging stage can be run for each module separately."""
		)
	else:
		group.add_argument(
			"--no-packaging",
			action="store_true",
			help="""Disable automatic packaging. [False]"""
		)

	group.add_argument(
		"--prokka-package-style",
		type=str,
		choices=["by_sample", "all_in_one"],
		default="by_sample",
		help="""Should the prokka annotation be packaged into one directory per sample (by_sample) or into one single directory (all_in_one)? [by_sample]"""
	)

	group.add_argument(
		"--make-ratt-data-tarballs",
		action="store_true",
		help="""ratt-based annotation transfers may generate large amounts of results. This option will generate a set of tarballs containing the transferred annotations. [False]"""
	)


def add_survey_parser(subparsers):
	parser = subparsers.add_parser(
		"survey",
		help=BGSURVEY_DESC,
		description=BGSURVEY_DESC
	)

	add_default_options(parser)
	add_survey_group(parser, shared_only=False)
	add_package_group(parser)
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
	parser.set_defaults(runmode="survey")


def add_asm_parser(subparsers):
	parser = subparsers.add_parser(
		"assemble",
		help=BGASM_DESC,
		description=BGASM_DESC
	)

	add_default_options(parser)
	add_asm_group(parser, shared_only=False)
	add_ann_group(parser, shared_only=False)
	add_package_group(parser)
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
	parser.set_defaults(runmode="assemble")


def add_ann_parser(subparsers):
	parser = subparsers.add_parser(
		"annotate", 
		help=BGANN_DESC,
		description=BGANN_DESC
	)

	add_default_options(parser)
	add_ann_group(parser, shared_only=False)
	add_package_group(parser)
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
	parser.set_defaults(runmode="annotate")


def add_package_parser(subparsers):
	parser = subparsers.add_parser(
		"package", 
		help=BGPACKAGE_DESC,                  	
		description=BGPACKAGE_DESC
	)

	add_default_options(parser)
	add_package_group(parser, is_packaging_module=True)
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
	parser.set_defaults(runmode="package")


def add_auto_parser(subparsers):
	parser = subparsers.add_parser(
		"auto",
		help=BGAUTO_DESC,
		description=BGAUTO_DESC
	)
	add_default_options(parser)
	add_survey_group(parser)
	add_asm_group(parser)
	add_ann_group(parser)
	add_package_group(parser)
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
	parser.set_defaults(runmode="auto")


def main():
	print("Starting EI bgrr| V " + __version__)
	from qaa import QAA_ID
	print("Using qaa " + QAA_ID)
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

	add_auto_parser(stage_parsers)

	args = bgrrl_parser.parse_args()
	BgrrlRunner(args).run()



if __name__ == "__main__":
	main()	
