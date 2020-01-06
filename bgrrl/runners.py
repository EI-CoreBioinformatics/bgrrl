import os
from os.path import join, dirname, basename, exists
import sys
from enum import Enum, unique
from collections import namedtuple, Counter
import glob

from eicore.snakemake_helper import *

from bgrrl.enterobase_helpers import validateEnterobaseInput, loadEnterobaseCriteria 
from bgrrl.bin.qc_eval import main as qc_eval_main
from bgrrl.bin.survey_stage_evaluation import main as survey_eval_main
from bgrrl.bin.asm_report import main as asm_report_main
from bgrrl.bin.ann_report import main as ann_report_main
from bgrrl.bin.asm_stage_report import main as asm_stage_report_main
from bgrrl.bin.annocmp import main as annocmp_main
from bgrrl.samplesheet import verifySamplesheet, Samplesheet, BaseSample, ASM_Sample, ANN_Sample
from bgrrl.bgrrl_config import BgrrlConfigurationManager

from qaa.runners import QaaRunner

class BgrrlModuleRunner:
	def __init__(self, module, config_manager):
		self.module = module
		self.config_manager = config_manager
		self.cols_to_verify = list()
		self.sampletype = BaseSample

	def verify_data(self, sheet, cols_to_verify, sampletype):
		if cols_to_verify:
			return Samplesheet(sheet, sampletype=sampletype).verifySampleData(fields=cols_to_verify)
		return True

	def run_module(self):

		def get_smk_file(module):
			return join(dirname(__file__), "zzz", module + ".smk.py")

		self.verify_data(self.config_manager.input_sheet, self.cols_to_verify, self.sampletype)
		self.config_manager._config["samplesheet"] = self.config_manager.input_sheet
		self.config_manager.module = self.module

		config_file = self.config_manager.generate_config_file(self.module)
		print("Run configuration file: " + config_file)

		print("Running " + self.module)
		snake = get_smk_file(self.module)
		
		return run_snakemake(
			snake,
			self.config_manager.output_dir,
			config_file,
			self.config_manager.exe_env,
			dryrun=False,
			unlock=self.config_manager.unlock
		)

class BGSurveyRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super(BGSurveyRunner, self).__init__(module, config_manager)
		self.cols_to_verify = ["R1", "R2", "S"]
		self.sampletype = BaseSample

	def run(self):
		readtype = "bbduk" if self.config_manager.no_normalization else "bbnorm"
		min_tadpole_size = str(int(self.config_manager.minimum_survey_assembly_size))

		qc_eval_args = [
			"--readtype", readtype,
			"--min_tadpole_size", min_tadpole_size,
			self.config_manager.output_dir
		]
		if self.config_manager.single_cell:
			qc_eval_args.insert(0, "--single-cell")

		if self.config_manager.report_only:
			run_result = survey_eval_main(qc_eval_args)
		else:
			run_result = self.run_module()
			if run_result:
				#qaa_args = self.config_manager.create_qaa_args(stage="qc_survey")
				self.config_manager.input_sheet = join(
					self.config_manager.output_dir, 
					"reports", 
					"samplesheets", 
					"samplesheet.survey_pass.yaml"
				)

				if run_result and self.config_manager.full_qaa_analysis:
					qaa_args = self.config_manager.create_qaa_args(stage="qc_report")
					run_result = QaaRunner(qaa_args).run()

				if run_result and not self.config_manager.no_packaging:
					self.config_manager.package_mode = "processed_reads"
					run_result = BGPackageRunner("bgpackage", self.config_manager).run_module()						


		return run_result


class BGAssemblyRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super(BGAssemblyRunner, self).__init__(module, config_manager)
		self.cols_to_verify = list()
		self.sampletype = ASM_Sample

	def run(self):
		eb_criteria = self.config_manager._config.get("enterobase_criteria", "")
		asm_stage_report_args = [
			self.config_manager.output_dir,
			join(self.config_manager.output_dir, "reports")	
		]
		asm_report_args = [
			self.config_manager.output_dir, 
			self.config_manager.enterobase_groups, 
			eb_criteria
		]

		if self.config_manager.report_only:
			run_result = asm_stage_report_main(asm_stage_report_args)
			if self.config_manager.enterobase_groups: # needs validation?
				run_result = asm_report_main(asm_report_args)
		else:
			run_result = self.run_module()
			if run_result:
				run_result = asm_stage_report_main(asm_stage_report_args)
				if run_result:

					if self.config_manager.run_annotation:
						print(
								"WARNING: Prokka annotation selected\n" + \
								"If your jobs fail, you might have to update tbl2asn and/or exclude nodes " + \
								"(hmmscan/GNU parallel fails)."
						)

						BGAnnotationRunner.check_prokka_nodes(
							join(self.config_manager.output_dir, "annotation", "prokka"),
							join(self.config_manager.output_dir, "reports", "ann_run_report.txt")
						)
						qaa_stage = "asm,ann"
						
						if self.config_manager._config["run_ratt"]: 
							BGAnnotationRunner.run_ratt_report(
								self.config_manager.ratt_reference,
								join(self.config_manager.output_dir, "annotation", "ratt"),
								join(self.config_manager.output_dir, "annotation", "prokka"),
                        	   	join(self.config_manager.output_dir, "reports")
							)
					else:
						qaa_stage = "asm"
				
					qaa_args = self.config_manager.create_qaa_args(stage=qaa_stage)
					run_result = QaaRunner(qaa_args).run()
					if run_result:
						if self.config_manager.enterobase_groups:
							run_result = asm_report_main(asm_report_args)
						if run_result and not self.config_manager.no_packaging:
							package_mode = qaa_stage + (",analysis" if self.config_manager.is_final_step or self.config_manager.run_annotation else "")
							self.config_manager.package_mode = package_mode # "asm"
							run_result = BGPackageRunner("bgpackage", self.config_manager).run_module()

		return run_result


class BGAnnotationRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super(BGAnnotationRunner, self).__init__(module, config_manager)
		self.cols_to_verify = ["Assembly"]
		self.sampletype = ANN_Sample

	@staticmethod
	def check_prokka_nodes(prokka_path, prokka_run_report):
		nodes = set()
		for f in glob.glob(join(prokka_path, "*", "PROKKA_FAILED")):
			nodes.add(open(f).read().strip())
		if nodes:
			with open(prokka_run_report, "at") as run_report:
				for node in sorted(nodes):
					print(node, file=run_report)
			print(
				"Failed prokka jobs were executed on nodes: {}.\n" + \
				"Try to exclude those nodes from being used for rule ann_prokka.".format(sorted(list(nodes))), 
			)                                                                                                  			

	@staticmethod
	def run_ratt_report(ratt_reference, ratt_dir, prokka_dir, report_dir):
		try:
			ann_report_main(["--ref-dir", ratt_reference, ratt_dir])
			annocmp_main([prokka_dir, ratt_dir, report_dir])
		except:
			open(join(report_dir, "annotation_report.tsv"), "at").close()

	def run(self):
		run_result = False

		if self.config_manager.report_only:
			run_result = False
			if self.config_manager._config["run_ratt"]: 
				BGAnnotationRunner.run_ratt_report(
					self.config_manager.ratt_reference,
					join(self.config_manager.output_dir, "annotation", "ratt"),
					join(self.config_manager.output_dir, "annotation", "prokka"),
					join(self.config_manager.output_dir, "reports")
				)
		else:
			print(
				"WARNING: Prokka annotation selected\n" + \
				"If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails)."
			)
			run_result = self.run_module()
			if not run_result:
				print("ANNOTATION RUN FAILED?")
			else:
				BGAnnotationRunner.check_prokka_nodes(
					join(self.config_manager.output_dir, "annotation", "prokka"),
					join(self.config_manager.output_dir, "reports", "ann_run_report.txt")
				)

				qaa_args = self.config_manager.create_qaa_args(stage="ann")
				run_result = QaaRunner(qaa_args).run()

				if self.config_manager._config["run_ratt"]: 
					BGAnnotationRunner.run_ratt_report(
						self.config_manager.ratt_reference,                          				
						join(self.config_manager.output_dir, "annotation", "ratt"),
						join(self.config_manager.output_dir, "annotation", "prokka"),
						join(self.config_manager.output_dir, "reports")
					)

				if run_result and not self.config_manager.no_packaging:
					self.config_manager.package_mode = "ann"
					run_result = BGPackageRunner("bgpackage", self.config_manager).run_module()
					if run_result:
						self.config_manager.package_mode = "analysis,ann"
						run_result = BGPackageRunner("bgpackage", self.config_manager).run_module()

		return run_result


class BGPackageRunner(BgrrlModuleRunner):
	def __init__(self, module, config_manager):
		super(BGPackageRunner, self).__init__(module, config_manager)
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
	def __init__(self, args):
		self.config_manager = BgrrlConfigurationManager(args)

	def __run_all(self):
		raise ValueError("Wrong runmode: (ATTEMPT_FULL is not implemented yet)")

	def run(self):

		modules = {
			"survey": ("bgsurvey", BGSurveyRunner),
			"assemble": ("bgasm", BGAssemblyRunner),
			"annotate": ("bgann", BGAnnotationRunner),
			"package": ("bgpackage", BGPackageRunner)
		}
		
		runner = modules.get(self.config_manager.runmode, None)
		if runner is None:
			raise ValueError("Not a valid runmode: " + self.config_manager.runmode)
		module, runner = runner
		run_result = runner(module, self.config_manager).run()
	
		print()
		if run_result:
			print("BGRR| completed successfully.")
		else:
			raise ValueError("BGRR| failed.  Please consult logs to debug run.")

		return True
