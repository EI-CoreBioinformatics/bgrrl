from collections import OrderedDict, namedtuple
import os 
from os.path import exists, dirname, basename, join
import sys
import yaml
import pathlib
import glob

from eicore.snakemake_helper import *
from eicore.config_manager import ConfigurationManager


class BgrrlConfigurationManager(ConfigurationManager):
	
	def __manage(self):
		cfg_d = {
			"etc": join(dirname(__file__), "..", "etc"),
			"cwd": os.getcwd(),
			"reapr_correction": False,
			"run_prokka": self.runmode == "annotate" or (hasattr(self, "run_annotation") and self.run_annotation),
			"run_ratt": False,
			"package_dir": self.package_dir,
			"single_cell_mode": self.single_cell	
		}
		self._config.update(cfg_d)
		
		if hasattr(self, "project_prefix"):
			prefix = self.project_prefix if self.project_prefix is not None else ""
			self._config["project_prefix"] = prefix
		
		if hasattr(self, "prokka_package_style"):
			self._config["prokka_package_style"] = self.prokka_package_style

		if hasattr(self, "contig_minlen"):
			self._config["asm_lengthfilter_contig_minlen"] = max(0, self.contig_minlen)

		self._config["run_ratt"] = hasattr(self, "ratt_reference") and self.ratt_reference is not None
		if self._config["run_ratt"]:
			if not exists(self.ratt_reference):
				raise ValueError("Invalid ratt reference location: " + self.ratt_reference)
			self._config["ratt_reference"] = self.ratt_reference
			if hasattr(self, "make_ratt_data_tarballs"):
				self._config["make_ratt_data_tarballs"] = self.make_ratt_data_tarballs		


	def __init__(self, ap_args):
		self.alt_hpc_config_warning = "Please run bginit or provide a valid HPC configuration file with --hpc_config."
		self.alt_config_warning = "Please run bginit or provide a valid configuration file with --bgrrl_config/--config."
		self.alt_multiqc_config_warning = "Please run bginit to obtain a valid MultiQC configuration file template."
		self.job_suffix_prefix = "bgrrl"

		super().__init__(ap_args)

		# check for input
		if not hasattr(self, "input"):
			try:
				self.input = self.input_sheet
			except:
				raise ValueError("Configuration has neither 'input' nor 'input_sheet' attribute.")

		self.__manage()

	def create_qaa_args(self, stage="init"):
		from qaa.qaa_args import QAA_ArgumentsAdapter as QAA_Args
		from .qaa_helpers import STAGE_QAA_ARGS

		qaa_args = QAA_Args(**STAGE_QAA_ARGS["init"])
		qaa_args.update(
			quast_mincontiglen=1000,
			project_prefix=self.project_prefix,
			config=self.config_file,
			hpc_config=self.hpc_config_file,
			multiqc_config=self.multiqc_config_file,	
			normalized=not self.no_normalization if hasattr(self, "no_normalization") else True,
			multiqc_dir=join(self.report_dir, "multiqc", stage.split(",")[0]),
			single_cell_mode = self.single_cell 
		)

		if not stage == "init":

			qaa_args.update(**vars(self))

			if stage == "asm,ann":
				qaa_args.update(**STAGE_QAA_ARGS["asm"])
				qaa_args.update(**{
					"qaa_mode": "genome,transcriptome,proteome", 
					"runmode": "asm,ann"
				})
			else:			
				try:
					qaa_args.update(**STAGE_QAA_ARGS[stage])
				except:
					raise ValueError("Invalid stage '{}' in BCM::create_qaa_args().".format(stage))
		
		return qaa_args
