import os
from os.path import join, dirname, basename
import csv
from collections import namedtuple
import yaml
import sys
from copy import copy
import pathlib

from eicore.snakemake_helper import *
from eicore.config_manager import ConfigurationManager

ETC_DIR = join(dirname(__file__), "..", "..", "..", "etc")

QAA_MODES = ["geno", "genome", "g", "tran", "transcriptome", "t", "prot", "proteome", "p", "all"]

class QaaConfigurationManager(ConfigurationManager):

	@staticmethod
	def calc_qualimap_memory_setting(hpc_config_file):
		""" from eimethyl: Extract qualimap mem setting from hpc_config, 
			convert to gigabytes and subtract a bit to account for memory above java heap """
		import json
		data = json.load(open(hpc_config_file))
		try:
			return "{}G".format(int(int(data["qaa_qualimap"]["memory"]) / 1000) - 2)
		except:
			raise ValueError("Could not find qaa_qualimap:memory record in {}.".format(hpc_config_file))

	def write_run_config(self):
		with open(join(self.output_dir, "qaa_conf.yaml"), "w") as conf_out:
			yaml.dump(self.run_config, conf_out, default_flow_style=False)

	
	def __manage(self):

		def check_qaa_modes(modes, is_survey_run=False):
			requested_modes = set(modes.split(","))
			valid_modes = requested_modes.intersection(QAA_MODES)
			if not valid_modes:
				raise ValueError(
					"No valid qaa-modes provided. Valid modes are {}.".format(",".join(QAA_MODES))
				)
			if "all" in requested_modes:
				requested_modes = {"all"}
			if requested_modes != valid_modes:
				invalid_modes = requested_modes.difference(QAA_MODES)
				print("Dropped invalid qaa-modes {}. Proceeding with {}.".format(
					",".join(invalid_modes),
					",".join(requested_modes)
				))
			return {
				"genome": is_survey_run or bool(
					{"geno", "g", "genome", "all"}.intersection(requested_modes)
				),
				"transcriptome": not is_survey_run and bool(
					{"tran", "t", "transcriptome", "all"}.intersection(requested_modes)
				),
				"proteome": not is_survey_run and bool(
					{"prot", "p", "proteome", "all"}.intersection(requested_modes)
				)
			}
		
		if not hasattr(self, "runmode"):
			self.runmode = "asm"

		qaa_modes = check_qaa_modes(self.qaa_mode, is_survey_run=self.runmode=="survey")

		print(*(vars(self).items()), sep="\n")

		align_reads = self.align_reads
		if type(align_reads) is str:
			align_reads = align_reads.lower() if align_reads.lower() != "no" else False

		self.run_config = {
			"samplesheet": dict(),
			"output_dir": self.output_dir,
			"etc_dir": ETC_DIR,
			"qualimap_mem": self.calc_qualimap_memory_setting(self.hpc_config_file),
			"align_reads": align_reads,
			"normalized": self.normalized if hasattr(self, "normalized") else False,
			"project_prefix": self.project_prefix if hasattr(self, "project_prefix") else "",
			"multiqc_dir": self.multiqc_dir if hasattr(self, "multiqc_dir") else ".",
			"multiqc_config": self.multiqc_config if hasattr(self, "multiqc_config") else ".", #?
			"survey_assembly": self.runmode == "survey",
			"quast_mincontiglen": self.quast_mincontiglen if hasattr(self, "quast_mincontiglen") else 0,
			"annotation": self.annotation if hasattr(self, "annotation") else None,
			"single_cell_mode": self.single_cell_mode if hasattr(self, "single_cell_mode") else False
		}
		self.run_config.update(self._config)

		for tool in ["multiqc", "blobtools", "qualimap", "busco"]:
			run_tool = "run_" + tool
			self.run_config[run_tool] = getattr(self, run_tool) if hasattr(self, run_tool) else True
		for mode, run in qaa_modes.items():
			run_mode = "run_{}_module".format(mode)
			self.run_config[run_mode] = run

		sheet = self.run_config["samplesheet"]
		annotation_path = join(self.output_dir, "annotation", "prokka") 
		if self.runmode == "survey" or self.runmode == "ann":
			for row in csv.reader(open(self.input), delimiter=","):
				sample = row[0]
				if self.runmode == "survey":
					sheet[sample] = {
						"assembly": join(self.output_dir, "qc", "tadpole", sample, sample + "_tadpole_contigs.fasta"),
						"r1": row[2],
						"r2": row[3]
					}
				else:
					sheet[sample] = {
						"transcripts": join(annotation_path, sample, sample + ".ffn"),
						"proteins": join(annotation_path, sample, sample + ".faa")
					}
				sheet[sample]["busco_db"] = self.busco_db
		else:
			try:
				sample_data = yaml.load(open(self.input_sheet), Loader=yaml.SafeLoader)
			except:
				sample_data = yaml.load(open(self.input), Loader=yaml.SafeLoader)

			for sample, data in sample_data.items():
				if self.runmode == "asm":
					sheet[sample] = {
						"assembly": join(self.output_dir, "assembly", sample, sample + ".assembly.fasta")
					}
				elif self.runmode == "asm,ann":
					sheet[sample] = {
						"assembly": join(annotation_path, sample, sample + ".fna"),
						"transcripts": join(annotation_path, sample, sample + ".ffn"),					
						"proteins": join(annotation_path, sample, sample + ".faa"),
						"busco_db": self.busco_db	
					}
				sheet[sample]["r1"] = data["bbduk_r1"]
				sheet[sample]["r2"] = data["bbduk_r2"]



	def __init__(self, ap_args):
		self.alt_hpc_config_warning = "Please provide a valid HPC configuration file with --hpc_config."
		self.alt_config_warning = "Please provide a valid configuration file with --config."

		self.job_suffix_prefix = "qaa"

		super().__init__(ap_args)

		self.__manage()
		self.write_run_config()
