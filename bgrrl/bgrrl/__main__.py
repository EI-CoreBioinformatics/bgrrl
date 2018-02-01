#!/usr/bin/env python3

import os
import sys
import shutil
import yaml
# https://stackoverflow.com/a/600612
import pathlib 
from argparse import ArgumentParser
from textwrap import dedent

from snakemake.utils import min_version

from . import NOW, DEFAULT_HPC_CONFIG_FILE, DEFAULT_BGRRL_CONFIG_FILE, PipelineStep, RunMode, __version__, ExecutionEnvironment, make_exeenv_arg_group
from bgrrl.bgrrl import run_qc, run_asm, compileQUASTReport, ENTERO_FILTER, TAX_FILTER, run_fin, compileBUSCO, run_ann
from bgrrl.bin.qc_eval import main as qc_eval_main
from bgrrl.bin.asm_report import main as asm_report_main

from qaa import TIME_CMD as tcmd
print("QAA_TIME_CMD="+tcmd) 


# min_version("4.0")


def main():
	print("Starting EI BGRRL V" + __version__)
	print()

	# bgrrl_config = yaml.load(open(DEFAULT_BGRRL_CONFIG_FILE))
	bgrrl_config_file = DEFAULT_BGRRL_CONFIG_FILE

	parser = ArgumentParser("The Earlham Institute Bacterial Genome Reconstruction & Recognition Pipeline (BGGR|)",
	                        description="""This program controls the EI-BGRRL pipeline via a number of Snakemake pipelines.""")

	parser.add_argument("input", help="""The (EI-PAP) samplesheet to process. This is a comma-separated file, containing location and meta-information for each sample to be processed.""")

	parser.add_argument("-o", "--output_dir", default="BGRRL_<timestamp>",
	                    help="If specified BGRRL will output data to this directory.")

	parser.add_argument("--project-prefix", default=os.path.basename(os.getcwd()), help="Reports and resultfiles/-folders will be prefixed by this.")

	parser.add_argument("-f", "--force", action="store_true", help="Force overwriting existing output directory, causes pipeline to be restarted.")

	parser.add_argument("-m", "--mode", choices=[ps.name.lower() for ps in PipelineStep], default="read_qc",
						help=dedent("""This option controls which part of the pipeline to run. Due to the task at hand being highly dependent on input data quality,
						it is recommended to run the 4 major steps of EI-BGRRL manually, with data checks, and/or additional manual data manipulation, in between.
						\"read_qc\" - This step preprocesses the read sets and assesses the assemblability of individual samples. It will produce an input-samplesheet for the assembly-step.
						\"assembly\" - This step performs genome assembly and post-assembly quality control for each read set.
						\"annotation\" - This step performs genome annotation and post-annotation quality control for each assembly.
						\"report_and_package\" - This step performs report generation and data packaging.
						\"attempt_full\" - If you know the data and what you're doing, you can attempt to run the 4 pipeline steps without breaks (not implemented yet)
						"""))

	parser.add_argument("--unlock", action='store_true', default=False,
	                    help="If snakemake is not running because it is reporting that the directory is locked, then you can unlock it using this option.  Please make sure that there are no other snakemake jobs running in this directory before using this option!")
	parser.add_argument("--bgrrl-config", help="Configuration file for BGRRL. This file specifies details for accessing services and commands to be executed prior to running each pipeline tool.  Default config file is: " + DEFAULT_BGRRL_CONFIG_FILE)

	parser.add_argument("--contig-minlen", type=int, default=0)
	parser.add_argument("--enterobase-groups", type=str, default="", help="Comma-separated list of Enterobase microbial organisms. The set of assemblies is tested against organism-specific criteria and assemblies are packaged according to their species. [NEEDS REWORDING!]. By default, the enterobase mode is disabled.")

	make_exeenv_arg_group(parser)	# Add in cluster and DRMAA options

	args = parser.parse_args()

	# Set run mode
	run_mode = PipelineStep[args.mode.upper()]

	# Establish a valid cluster configuration... may throw if invalid
	print("Configuring execution environment ... ", end="", flush=True)
	exe_env = ExecutionEnvironment(args, NOW, job_suffix=args.input + "_" + args.output_dir, log_dir=os.path.join(args.output_dir, "hpc_logs"))
	print("done.")
	print(str(exe_env))
	print()

	if args.bgrrl_config:
		print("Custom BGRRL configuration file specifed, overriding defaults ... ", end="", flush=True)
		bgrrl_config = yaml.load(open(args.bgrrl_config))
		bgrrl_config_file = args.bgrrl_config
		print("done.")
		print()
	else:
		print("Loading default BGRRL configuration ...", end="", flush=True)
		bgrrl_config = yaml.load(open(bgrrl_config_file))
		print("done.")
		print()

	if os.path.exists(args.output_dir):
		if args.force:
			print("Output directory already exists and overwrite was requested (-f option).  Deleting directory contents ... ",
			      end='',
			      flush=True)
			shutil.rmtree(args.output_dir)
			os.makedirs(args.output_dir)
			print("done.")
		else:
			print("Output already exists, attempting to resume.", flush=True)
	else:
		print("Output directory doesn't exist creating ... ", end="", flush=True)
		os.makedirs(args.output_dir)
		print("done.")

	logs_dir = os.path.join(args.output_dir, "hpc_logs")
	if not os.path.exists(logs_dir) and exe_env.use_scheduler:
		print("HPC log dir doesn't exist.  Creating " + logs_dir + " now ... ", end="", flush=True)
		os.makedirs(logs_dir)
		print("done.")

	print()

	"""
        run_dir = None
	if args.run_dir:
		run_dir = args.run_dir

	if not run_dir:
		raise ValueError("Illegal configuration: No run directory established.")
        """
	print(args.mode.upper())
	if run_mode == PipelineStep.READ_QC:
		run_result = run_qc(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)	
		qc_eval_main([args.input, args.output_dir])
	elif run_mode == PipelineStep.ASSEMBLY:
		# run_result = True
		run_result = run_asm(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
		asm_report_main([args.output_dir, args.enterobase_groups])
		"""
		report_dir = os.path.join(args.output_dir, "reports")
		pathlib.Path(report_dir).mkdir(parents=True, exist_ok=True)
		with open(os.path.join(report_dir, "quast_report.tsv"), "w") as qout, open(os.path.join(report_dir, "quast_report.enterobase.tsv"), "w") as qeout, open(os.path.join(report_dir, "blobtools_taxonomy_report.tsv"), "w") as teout, open(os.path.join(report_dir, "Salmonella_EB_samples.txt"), "w") as vout:
			entero_pass_assembly = set(ENTERO_FILTER(compileQUASTReport(os.path.join(args.output_dir, "qa", "quast"), out=qout), organism="Salmonella", out=qeout))
			entero_pass_taxonomy = set(TAX_FILTER(os.path.join(args.output_dir, "qa", "blobtools", "blob"), organism="Salmonella", out=teout))
			print(*sorted(entero_pass_assembly.intersection(entero_pass_taxonomy)), sep="\n", file=vout)
		"""
				
	elif run_mode == PipelineStep.ANNOTATION:
		run_result = run_ann(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
	elif run_mode == PipelineStep.FINALIZE:
		run_result = run_fin(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)

		with open(os.path.join(report_dir, "busco_report.tsv"), "w") as bout:
			compileBUSCO(os.path.join(args.output_dir, "qa", "busco"), out=bout)
	else:
		print("Wrong runmode: (ATTEMPT_FULL is not implemented yet)", run_mode)
		exit(1)

	# result = eipap.illumina.pap.do_illumina_pap(jira, out_dir, run_dir, args, run_mode, exe_env, pap_config,
	#	                                            pap_per_lane=args.pap_per_lane)
	# result = True

	print()
	if run_result:
		print("BGRRL completed successfully.")
	else:
		print("BGRRL failed.  Please consult logs to debug run.")
		exit(1)


if __name__ == "__main__":
	main()
