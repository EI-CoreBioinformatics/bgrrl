import pkg_resources

__title__ = "bgrr|"
__author__ = "Christian Schudoma"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2017-2018 Earlham Institute"
__version__ = pkg_resources.require("bgrrl")[0].version

import os
from os.path import join, dirname, basename, exists
import sys
from enum import Enum, unique
from collections import namedtuple, Counter
import yaml
import csv
import shutil
# https://stackoverflow.com/a/600612
import pathlib
from copy import copy
import glob

from .snakemake_helper import *
from .workflow_runner import WorkflowRunner

from bgrrl.enterobase_helpers import validateEnterobaseInput, loadEnterobaseCriteria 
from bgrrl.bin.qc_eval import main as qc_eval_main
from bgrrl.bin.asm_report import main as asm_report_main
from bgrrl.bin.ann_report import main as ann_report_main
from bgrrl.bin.asm_stage_report import main as asm_stage_report_main
from bgrrl.bin.annocmp import main as annocmp_main
from bgrrl.samplesheet import verifySamplesheet, Samplesheet, BaseSample, ASM_Sample

from .qaa_helpers import QAA_ArgumentManager

from qaa import QAA_Runner, QAA_ID
print("QAA_ID="+QAA_ID)


TIME_CMD = " /usr/bin/time -v"

@unique
class PipelineStep(Enum):
	READ_QC = 0
	ASSEMBLY = 1
	ANNOTATION = 2
	DATA_QA = 3
	FINALIZE = 4
	ATTEMPT_FULL = 5


class BGRRLModuleRunner(object):
	def __init__(self, module, args, exe_env, hpc_config_file, config=dict()):
		print(config)
		self.config = dict(config)
		self.module = module
		self.outdir = args.output_dir
		self.unlock = args.unlock
		self.exe_env = exe_env
		self.config["hpc_config"] = hpc_config_file

		try:
			sample_sheet = Samplesheet(args.input, sampletype=ASM_Sample)
		except:
			print("This samplesheet {} is not an ASM_Sample. Trying to parse it as BaseSample.")
			sample_sheet = Samplesheet(args.input, sampletype=BaseSample)

		if sample_sheet.verifySampleData(fields=["R1", "R2"]):
			self.config["samplesheet"] = args.input

		self.config["out_dir"] = self.outdir

		for k, v in args._get_kwargs():
			# need to find a solution here how to not override already sanitised input, e.g. enterobase-groups
			# this should do it for the moment, but the whole args-structure needs to be changed
			if k not in self.config:
				print("WRITING {}: {} TO CONFIG".format(k, v))
				self.config[k] = v

		self.config_file = os.path.join(self.outdir, module + ".conf.yaml")
		with open(self.config_file, "w") as conf_out:
			yaml.dump(self.config, conf_out, default_flow_style=False)

	def run(self):
 
		print("Running " + self.module)
		snake = os.path.join(os.path.dirname(__file__), "zzz", self.module + ".smk.py")
		run_result = run_snakemake(snake, self.outdir, self.config_file, self.exe_env, dryrun=False, unlock=self.unlock)

		return run_result


class BGRRLRunner(WorkflowRunner):

	def __init__(self, args, **kwargs):

		self.args = copy(args)
	
		args.alt_hpc_config_warning = "Please run bginit or provide a valid HPC configuration file with --hpc_config."
		args.alt_config_warning = "Please run bginit or provide a valid configuration file with --bgrrl_config/--config."
		args.config = args.bgrrl_config

		super().__init__(args)

		# Set run mode
		print(args.module.upper())
		self.run_mode = PipelineStep[args.module.upper()]

	def __run_qc(self): 
		readtype = "bbduk" if self.args.no_normalization else "bbnorm"

		if self.args.report_only:
			run_result = qc_eval_main(["--readtype", readtype, self.args.input, self.args.output_dir])
		else:
			run_result = BGRRLModuleRunner("bgrrl-qc", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
			if run_result:
				qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="qc_survey")
				qaa_run = QAA_Runner(qaa_args).run()					
				if qaa_run:
					run_result = qc_eval_main(["--readtype", readtype, self.args.input, self.args.output_dir])
					if run_result:
						self.args.input = join(self.args.output_dir, "reports", "samplesheets", "samplesheet.qc_pass.tsv")
						qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="qc_report")
						qaa_run = QAA_Runner(qaa_args).run()

		return run_result


	def __run_asm(self):

		self.config["etc"] = os.path.join(os.path.dirname(__file__), "etc")
		self.config["cwd"] = os.getcwd()
		self.config["reapr_correction"] = False
	
		if self.args.contig_minlen: # this is really, really, really bad! 
			self.config["use_asm_lengthfilter"] = True
			self.config["asm_lengthfilter_contig_minlen"] = self.args.contig_minlen
		else:
			self.config["use_asm_lengthfilter"] = False
			self.config["asm_lengthfilter_contig_minlen"] = 0

		eb_criteria = self.config.get("enterobase_criteria", "")

		if self.args.report_only:
			run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
			if self.args.enterobase_groups: # needs validation?
				run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
		else:
			run_result = BGRRLModuleRunner("bgrrl-asm", self.args, self.exe_env, self.hpc_config_file, config=self.config).run() 
			if run_result:
				run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
				if run_result:
						qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="asm")
						qaa_run = QAA_Runner(qaa_args).run()
						if qaa_run:
							self.args.package_mode = "asm"
							if self.args.enterobase_groups:
								run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
							if run_result and not self.args.no_packaging:
								run_result = self.__run_fin() 

		return run_result

	def __run_ann(self):
		self.config["etc"] = os.path.join(os.path.dirname(__file__), "etc")
		self.config["cwd"] = os.getcwd()
	
		self.config["run_ratt"] = self.args.annotation == "both"
		self.config["run_prokka"] = self.args.annotation in ("both", "prokka")
		self.config["ratt_reference"] = self.args.ratt_reference_dir
		self.config["prokka_package_style"] = self.args.prokka_package_style 

		assert not self.config["run_ratt"] or os.path.exists(self.config["ratt_reference"]), "Missing reference data for ratt. Please make sure to use the --ratt-reference-dir parameter."
		run_result = False
		print(self.args)

		if self.args.report_only:
			run_result = False
			if self.args.annotation == "both":
				run_result = ann_report_main(["--ref-dir", self.args.ratt_reference_dir, join(self.args.output_dir, "annotation", "ratt")])
				annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
		else:
			if self.args.annotation in ("prokka", "both"):
				print("WARNING: Prokka annotation selected. If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails).")
				run_result = BGRRLModuleRunner("bgrrl-ann", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
				if not run_result:
					print("ANNOTATION RUN FAILED?")
				if run_result:

					with open(join(self.args.output_dir, "reports", "ann_run_report.txt"), "at") as run_report:
						nodes = set()
						for f in glob.glob(join(self.args.output_dir, "annotation", "prokka", "*", "PROKKA_FAILED")):
							node = open(f).read().strip()
							nodes.add(node)
							print(basename(dirname(f)), node, sep="\t", file=run_report)
						print("Failed prokka jobs were executed on nodes: {}. Try to exclude those nodes from being used for rule ann_prokka.".format(sorted(list(nodes))), file=run_report)


					qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="ann")
					qaa_run = QAA_Runner(qaa_args).run()
	
					if qaa_run and self.args.annotation == "both":
						ann_report_main(["--ref-dir", self.args.ratt_reference_dir, join(self.args.output_dir, "annotation", "ratt")])
						annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
					else:
						open(join(self.args.output_dir, "reports", "annotation_report.tsv"), "at")
					if qaa_run and not self.args.no_packaging:
						self.args.package_mode = "ann"
						run_result = self.__run_fin() 
		return run_result

	def __run_fin(self):
		self.config["package_dir"] = os.path.join(os.path.dirname(self.args.output_dir), "Data_Package")
		self.config["etc"] = os.path.join(os.path.dirname(__file__), "etc")
		self.config["cwd"] = os.getcwd()	
		self.config["run_ratt"] = self.args.annotation == "both"
		self.config["run_prokka"] = self.args.annotation in ("both", "prokka")
		self.config["prokka_package_style"] = self.args.prokka_package_style 

		if "enterobase_groups" in self.args: 
			try:
				eb_criteria = loadEnterobaseCriteria(self.config["enterobase_criteria"])
			except:
				print("Enterobase-groups selected but missing enterobase_criteria entry in config, please add path to enterobase criteria.", file=sys.stderr)
				sys.exit(1)
			self.config["enterobase_groups"] = validateEnterobaseInput(self.args.enterobase_groups, eb_criteria)
		else:
			self.config["enterobase_groups"] = list()
		if "project_prefix" in self.args and self.args.project_prefix:
			self.config["misc"]["project"] = self.args.project_prefix
		print("_FIN_CONFIG")
		print(self.config)

		run_result = BGRRLModuleRunner("bgrrl-fin", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
		return run_result

	def __run_all(self):
		raise ValueError("Wrong runmode: (ATTEMPT_FULL is not implemented yet)")

	def run(self):
	
		if self.run_mode == PipelineStep.READ_QC:
			run_result = self.__run_qc()
		elif self.run_mode == PipelineStep.ASSEMBLY:
			run_result = self.__run_asm() 
		elif self.run_mode == PipelineStep.ANNOTATION:
			run_result = self.__run_ann()
		elif self.run_mode == PipelineStep.FINALIZE:
			run_result = self.__run_fin() 
		else:
			run_result = self.__run_all()		
	
		print()
		if run_result:
			print("BGRR| completed successfully.")
		else:
			raise ValueError("BGRR| failed.  Please consult logs to debug run.")

		return True


