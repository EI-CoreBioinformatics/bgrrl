import os
from os.path import join, dirname, basename, exists
import sys
from enum import Enum, unique
from collections import namedtuple, Counter
import glob

from eicore.snakemake_helper import *

from bgrrl.enterobase_helpers import validateEnterobaseInput, loadEnterobaseCriteria 
from bgrrl.bin.asm_report import main as asm_report_main
from bgrrl.bin.ann_report import main as ann_report_main
from bgrrl.bin.annocmp import main as annocmp_main
from bgrrl.samplesheet import verifySamplesheet, Samplesheet, BaseSample, ASM_Sample, ANN_Sample
from bgrrl.bgrrl_config import BgrrlConfigurationManager

from qaa.runners import QaaRunner

class BgrrlModuleRunner:
	def __init__(self, module, config_manager):
		self.snake = join(dirname(__file__), "zzz", module + ".smk.py")
		self.module = module
		self.config_manager = config_manager
		self.cols_to_verify = list()
		self.sampletype = BaseSample

	def verify_data(self, sheet, cols_to_verify, sampletype):
		if cols_to_verify:
			return Samplesheet(sheet, sampletype=sampletype).verifySampleData(fields=cols_to_verify)
		return True

	def run_module(self):
		self.verify_data(self.config_manager.input_sheet, self.cols_to_verify, self.sampletype)
		self.config_manager._config["samplesheet"] = self.config_manager.input_sheet
		self.config_manager.module = self.module

		config_file = self.config_manager.generate_config_file(self.module)
		print("Run configuration file: " + config_file)

		print("Running " + self.module)
		return run_snakemake(
			self.snake,
			self.config_manager.output_dir,
			config_file,
			self.config_manager.exe_env,
			dryrun=False,
			unlock=self.config_manager.unlock
		)

class BgrrlSurveyRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super().__init__(module, config_manager)
		self.cols_to_verify = ["R1", "R2", "S"]
		self.sampletype = BaseSample

	def run(self):
		run_result = self.run_module()
		if run_result:
			#qaa_args = self.config_manager.create_qaa_args(stage="qc_survey")
			self.config_manager.input_sheet = join(
				self.config_manager.output_dir, 
				"reports", 
				"samplesheets", 
				"samplesheet.survey_pass.yaml"
			)

			if self.config_manager.full_qaa_analysis:
				qaa_args = self.config_manager.create_qaa_args(stage="qc_report")
				run_result = QaaRunner(qaa_args).run()

			if run_result and not self.config_manager.no_packaging:
				self.config_manager.package_mode = "processed_reads"
				run_result = BgrrlPackageRunner("bgpackage", self.config_manager).run_module()						


		return run_result


class BgrrlAssemblyRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super().__init__(module, config_manager)
		self.cols_to_verify = list()
		self.sampletype = ASM_Sample

	def run(self):
		asm_report_args = [
			self.config_manager.output_dir, 
			self.config_manager.enterobase_groups, 
			self.config_manager._config.get("enterobase_criteria", "")
		]

		run_result = self.run_module()
		qaa_stage = "asm"
		if run_result:
			if self.config_manager.run_annotation:
				print(
					"WARNING: Prokka annotation selected\n" + \
					"If your jobs fail, you might have to update tbl2asn and/or exclude nodes " + \
					"(hmmscan/GNU parallel fails)."
				)
				qaa_stage = "asm,ann" if self.config_manager._config["run_ratt"] else "asm"
				
			qaa_args = self.config_manager.create_qaa_args(stage=qaa_stage)
			run_result = QaaRunner(qaa_args).run()
			if run_result:
				if self.config_manager.enterobase_groups:
					run_result = asm_report_main(asm_report_args)
				if run_result and not self.config_manager.no_packaging:
					package_mode = qaa_stage + (",analysis" if self.config_manager.is_final_step or self.config_manager.run_annotation else "")
					self.config_manager.package_mode = package_mode # "asm"
					run_result = BgrrlPackageRunner("bgpackage", self.config_manager).run_module()

		return run_result


class BgrrlAnnotationRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super().__init__(module, config_manager)
		self.cols_to_verify = ["Assembly"]
		self.sampletype = ANN_Sample

	@staticmethod
	def run_ratt_report(ratt_reference, ratt_dir, prokka_dir, report_dir):
		try:
			ann_report_main(["--ref-dir", ratt_reference, ratt_dir])
			annocmp_main([prokka_dir, ratt_dir, report_dir])
		except:
			open(join(report_dir, "annotation_report.tsv"), "at").close()

	def run(self):
		run_result = False

		print(
			"WARNING: Prokka annotation selected\n" + \
			"If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails)."
		)
		run_result = self.run_module()
		if run_result:
			qaa_args = self.config_manager.create_qaa_args(stage="ann")
			run_result = QaaRunner(qaa_args).run()

			if self.config_manager._config["run_ratt"]: 
				BgrrlAnnotationRunner.run_ratt_report(
					self.config_manager.ratt_reference,                          				
					join(self.config_manager.output_dir, "annotation", "ratt"),
					join(self.config_manager.output_dir, "annotation", "prokka"),
					join(self.config_manager.output_dir, "reports")
				)

			if run_result and not self.config_manager.no_packaging:
				self.config_manager.package_mode = "ann"
				run_result = BgrrlPackageRunner("bgpackage", self.config_manager).run_module()
				if run_result:
					self.config_manager.package_mode = "analysis,ann"
					run_result = BgrrlPackageRunner("bgpackage", self.config_manager).run_module()

		return run_result


class BgrrlPackageRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super().__init__(module, config_manager)
		if self.config_manager.package_mode in ("processed_reads", "asm", "asm,analysis", "asm,ann,analysis"):
			self.cols_to_check = list()
			self.sampletype = ASM_Sample
		elif self.config_manager.package_mode in ("ann", "ann,analysis", "analysis,ann"):
			self.cols_to_check = ["Assembly"]
			self.sampletype = ANN_Sample


	def run(self):
		req_pmodes = set(self.config_manager.package_mode.split(","))
		req_valid_pmodes = req_pmodes.intersection({"ann", "asm", "analysis", "processed_reads"})
		if req_pmodes != req_valid_pmodes:
			invalid_modes = req_modes.difference(req_valid_pmodes)
			print("Warning: Dropping invalid modes {} from requested package modes.".format(invalid_modes))
			if not req_valid_pmodes:
				raise ValueError("No valid package modes requested: {}.".format(req_pmodes))
		self.config_manager.package_mode = tuple(req_valid_pmodes)

		if hasattr(self.config_manager, "enterobase_groups"):
			try:
				eb_criteria = loadEnterobaseCriteria(self.config_manager._config["enterobase_criteria"])
			except:
				raise ValueError(
					"""Enterobase-groups selected but missing enterobase_criteria entry in config, 
					please add path to enterobase criteria."""
				)
			self.config_manager._config["enterobase_groups"] = validateEnterobaseInput(
				self.config_manager.enterobase_groups, 
				eb_criteria
			)
		else:
			self.config_manager._config["enterobase_groups"] = list()

		run_result = self.run_module()
		return run_result

class BgrrlRunner:
	modules = {
		"survey": ("bgsurvey", BgrrlSurveyRunner),
		"assemble": ("bgasm", BgrrlAssemblyRunner),
		"annotate": ("bgann", BgrrlAnnotationRunner),
		"package": ("bgpackage", BgrrlPackageRunner),
	}
	def __init__(self, args):
		self.config_manager = BgrrlConfigurationManager(args)
	def run(self):
		runner = self.modules.get(self.config_manager.runmode, None)
		if runner is None and self.config_manager.runmode != "auto":
			raise ValueError("Not a valid runmode: {}".format(self.config_manager.runmode))
		elif self.config_manager.runmode == "auto":
			module, runner = self.modules["survey"]
			run_result = runner(module, self.config_manager).run()
			if run_result:
				module, runner = self.modules["assemble"]
				run_result = runner(module, self.config_manager).run()
		else:
			module, runner = runner
			run_result = runner(module, self.config_manager).run()

		print()
		if run_result:
			print("Bgrr| completed successfully.")
		else:
			raise ValueError("Bgrr| failed in module:{}. Please consult logs to troubleshoot run.".format(module))

		return True
