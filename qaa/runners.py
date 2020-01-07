import os
from os.path import join, dirname, basename
import csv
from collections import namedtuple
import yaml
import sys
from copy import copy
import pathlib

from eicore.snakemake_helper import *
from .qaa_environment import QaaEnvironment
from .qaa_config import QaaConfigurationManager

from qaa.reporting.busco_report import compileBUSCOReport
from qaa.reporting.quast_report import compileQUASTReport
from qaa.reporting.blobtools_report import compileBlobReport
from qaa.reporting.sixteen_s_reporter import SixteenSReporter, HEADER as sixteenS_header

class QaaRunner:
	def __init__(self, args):
		self.config_manager = QaaConfigurationManager(args)
		self.qaa_env = QaaEnvironment(self.config_manager.run_config)
		self.config = self.config_manager.run_config
	def __clean_blobtools_trash(self):
		import glob	
		for f in glob.glob("*.bam.cov"):
			try:
				if os.path.isfile(f):
					os.unlink(f)
			except Exception as e:
				print(e)

	def report(self):
		def report_func(qa_dir, report_out, rfunc):
			with open(report_out, "w") as rep_out:
				rfunc(qa_dir, out=rep_out)

		if self.config_manager.runmode == "survey":
			#report_func(
			#	self.qaa_env.quast_dir, 
			#	join(self.config_manager.report_dir, "quast_survey_report.tsv"), 
			#	compileQUASTReport
			#)
			if self.config["run_busco"]:
				report_func(
					self.qaa_env.busco_geno_dir, 
					join(self.config_manager.report_dir, "busco_survey_report.tsv"),
					compileBUSCOReport
				)
			if self.config["run_blobtools"]:
				report_func(
					self.qaa_env.blob_dir,
					join(self.config_manager.report_dir, "blobtools_survey_report.tsv"),
					compileBlobReport
				)
		else:
			if "asm" in self.config_manager.runmode:
				if self.config["run_genome_module"]:
					#report_func(
					#	self.qaa_env.quast_dir,
					#	join(self.config_manager.report_dir, "quast_report.tsv"),
					#	compileQUASTReport
					#)
					pass
				if self.config["run_busco"]:
					report_func(
						self.qaa_env.busco_geno_dir,
						join(self.config_manager.report_dir, "busco_genome_report.tsv"),
						compileBUSCOReport
					)
				if self.config["run_blobtools"]:
					report_func(
						self.qaa_env.blob_dir,
						join(self.config_manager.report_dir, "blobtools_report.tsv"),
						compileBlobReport
					)
			if "ann" in self.config_manager.runmode:
				if self.config["run_transcriptome_module"] or self.config["run_proteome_module"]:
					report_func(
						self.qaa_env.busco_dir,
						join(self.config_manager.report_dir, "busco_report.tsv"),
						compileBUSCOReport
					)
				with open(join(self.config_manager.report_dir, "16S_report.tsv"), "wt") as ostream:
					SixteenSReporter(
						self.qaa_env.prokka_dir, 
						header=sixteenS_header, 
						out=ostream
					).generateReport()
	def run(self):
		run_result = run_snakemake(
			join(dirname(__file__), "zzz", "qaa.smk.py"),
			self.config_manager.output_dir,
			join(self.config_manager.output_dir, "qaa_conf.yaml"), 
			self.config_manager.exe_env,
			dryrun=False,
			unlock=self.config_manager.unlock
		)
		self.report()
		self.__clean_blobtools_trash()
		return run_result
