import pkg_resources

__title__ = "bgrr|"
__author__ = "Christian Schudoma (cschu)"
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
			sample_sheet = Samplesheet(args.input_sheet, sampletype=ASM_Sample)
		except:
			print("This samplesheet {} is not an ASM_Sample. Trying to parse it as BaseSample.")
			sample_sheet = Samplesheet(args.input_sheet, sampletype=BaseSample)

		if sample_sheet.verifySampleData(fields=["R1", "R2"]):
			self.config["samplesheet"] = args.input_sheet

		self.config["out_dir"] = self.outdir

		for k, v in sorted(vars(args).items()):
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

class BGRRLConfigManager(object):
	@staticmethod
	def manageArgs(args, config): 
		config["etc"] = join(dirname(__file__), "etc")
		config["cwd"] = os.getcwd()
		config["reapr_correction"] = False

		if hasattr(args, "contig_minlen"):
			config["asm_lengthfilter_contig_minlen"] = max(0, args.contig_minlen)
			config["use_asm_lengthfilter"] = config["asm_lengthfilter_contig_minlen"] > 0 # probably not needed anymore

		config["run_prokka"] = True
		config["run_ratt"] = False
		if hasattr(args, "ratt_reference") and args.ratt_reference is not None:
			config["run_ratt"] = True
			config["ratt_reference"] = args.ratt_reference
			if not os.path.exists(config["ratt_reference"]):
				raise ValueError("ratt reference location is invalid: " + config["ratt_reference"])
			if hasattr(args, "make_ratt_data_tarballs"):
				config["make_ratt_data_tarballs"] = args.make_ratt_data_tarballs

		if hasattr(args, "prokka_package_style"):
			config["prokka_package_style"] = args.prokka_package_style


		config["package_dir"] = join(dirname(args.output_dir), "Data_Package")
		if hasattr(args, "project_prefix") and args.project_prefix:
			config["misc"]["project"] = args.project_prefix

		pass



class BGRRLRunner(WorkflowRunner):

	def __init__(self, args):

		self.args = copy(args)
	
		args.alt_hpc_config_warning = "Please run bginit or provide a valid HPC configuration file with --hpc_config."
		args.alt_config_warning = "Please run bginit or provide a valid configuration file with --bgrrl_config/--config."

		super().__init__(args)

		BGRRLConfigManager.manageArgs(args, self.config)


	def __run_survey(self): 
		readtype = "bbduk" if self.args.no_normalization else "bbnorm"
		min_tadpole_size = str(int(self.args.minimum_survey_assembly_size))

		if self.args.report_only:
			run_result = qc_eval_main(["--readtype", readtype, "--min_tadpolesize", min_tadpole_size, self.args.input_sheet, self.args.output_dir])
		else:
			run_result = BGRRLModuleRunner("bgsurvey", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
			if run_result:
				qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="qc_survey")
				qaa_run = QAA_Runner(qaa_args).run()					
				if qaa_run:
					run_result = qc_eval_main(["--readtype", readtype, "--min_tadpolesize", min_tadpole_size, self.args.input_sheet, self.args.output_dir])
					
					self.args.input_sheet = join(self.args.output_dir, "reports", "samplesheets", "samplesheet.qc_pass.tsv")

					if run_result and self.args.full_qaa_analysis:
						qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="qc_report")
						run_result = QAA_Runner(qaa_args).run()

					if run_result and not self.args.no_packaging:
						self.args.package_mode = "processed_reads"
						run_result = self.__run_package()

		return run_result


	def __run_asm(self):

		eb_criteria = self.config.get("enterobase_criteria", "")

		if self.args.report_only:
			run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
			if self.args.enterobase_groups: # needs validation?
				run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
		else:
			run_result = BGRRLModuleRunner("bgasm", self.args, self.exe_env, self.hpc_config_file, config=self.config).run() 
			if run_result:
				run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
				if run_result:
						qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="asm")
						run_result = QAA_Runner(qaa_args).run()
						if run_result:
							if self.args.enterobase_groups:
								run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
							if run_result and not self.args.no_packaging:
								self.args.package_mode = "asm"
								run_result = self.__run_package() 
								if run_result and self.args.is_final_step:
									self.args.package_mode = "analysis"
									run_result = self.__run_package()

		return run_result

	def __run_ann(self):
		run_result = False

		if self.args.report_only:
			run_result = False
			if self.args.annotation == "both":
				run_result = ann_report_main(["--ref-dir", self.args.ratt_reference, join(self.args.output_dir, "annotation", "ratt")])
				annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
		else:
			if True: #self.args.annotation in ("prokka", "both"):
				print("WARNING: Prokka annotation selected. If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails).")
				run_result = BGRRLModuleRunner("bgann", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
				if not run_result:
					print("ANNOTATION RUN FAILED?")
				else:
					with open(join(self.args.output_dir, "reports", "ann_run_report.txt"), "at") as run_report:
						nodes = set()
						for f in glob.glob(join(self.args.output_dir, "annotation", "prokka", "*", "PROKKA_FAILED")):
							node = open(f).read().strip()
							nodes.add(node)
							print(basename(dirname(f)), node, sep="\t", file=run_report)
						print("Failed prokka jobs were executed on nodes: {}. Try to exclude those nodes from being used for rule ann_prokka.".format(sorted(list(nodes))), file=run_report)

					qaa_args = QAA_ArgumentManager.get_qaa_args(self.args, self.config_file, self.hpc_config_file, stage="ann")
					qaa_run = QAA_Runner(qaa_args).run()
	
					if qaa_run and hasattr(self.args, "ratt_reference") and self.args.ratt_reference is not None:
						ann_report_main(["--ref-dir", self.args.ratt_reference, join(self.args.output_dir, "annotation", "ratt")])
						annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
					else:
						open(join(self.args.output_dir, "reports", "annotation_report.tsv"), "at").close()

					if qaa_run and not self.args.no_packaging:
						self.args.package_mode = "ann"
						run_result = self.__run_package()
						if run_result:
							self.args.package_mode = "analysis"
							run_result = self.__run_package()

		return run_result

	def __run_package(self):
	
		req_pmodes = set(self.args.package_mode.split(","))
		req_valid_pmodes = req_pmodes.intersection({"ann", "asm", "analysis", "processed_reads"})
		if req_pmodes != req_valid_pmodes:
			print("Warning: Dropping invalid modes {} from requested package modes.".format(req_modes.difference(req_valid_pmodes)))
			if not req_valid_pmodes:
				raise ValueError("No valid package modes requested: {}.".format(req_pmodes))
		self.args.package_mode = tuple(req_valid_pmodes)

		if "enterobase_groups" in self.args: 
			try:
				eb_criteria = loadEnterobaseCriteria(self.config["enterobase_criteria"])
			except:
				raise ValueError("Enterobase-groups selected but missing enterobase_criteria entry in config, please add path to enterobase criteria.")
			self.config["enterobase_groups"] = validateEnterobaseInput(self.args.enterobase_groups, eb_criteria)
		else:
			self.config["enterobase_groups"] = list()

		run_result = BGRRLModuleRunner("bgpackage", self.args, self.exe_env, self.hpc_config_file, config=self.config).run()
		return run_result

	def __run_all(self):
		raise ValueError("Wrong runmode: (ATTEMPT_FULL is not implemented yet)")

	def run(self):
		if self.args.runmode == "survey":
			run_result = self.__run_survey()
		elif self.args.runmode == "assemble":
			run_result = self.__run_asm() 
		elif self.args.runmode == "annotate":
			run_result = self.__run_ann()
		elif self.args.runmode == "package":
			run_result = self.__run_package() 
		else:
			run_result = self.__run_all()		
	
		print()
		if run_result:
			print("BGRR| completed successfully.")
		else:
			raise ValueError("BGRR| failed.  Please consult logs to debug run.")

		return True


