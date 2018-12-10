from collections import OrderedDict, namedtuple
import os 
from os.path import exists, dirname, basename, join
import sys
import yaml
import pathlib
import glob

from .snakemake_helper import *

ExecutionEnvironmentArguments = namedtuple(
	"ExecutionEnvironmentArguments",
	[
		"scheduler",
		"partition",
		"no_drmaa",
		"max_nodes",
		"max_cores",
		"hpc_config"
	]
)


class ConfigurationManager(OrderedDict):
	
	def __handle_config_file(self, config_file, config_type, file_pattern, warning=""):
		try:
			init_file = glob.glob(self.config_dir, file_pattern)[0]
		except:
			init_file = ""

		if hasattr(self, config_file) and exists(getattr(self, config_file)):
			print("Custom {} file specified, overriding defaults.".format(config_type))
			setattr(self, config_file + "_file", getattr(self, config_file))
		elif init_file:
			print("Found {} file at init location, using this.".format(config_type))
			setattr(self, config_file + "_file", init_file)
		else:
			raise ValueError(
				"No valid {} file specified.{}".format(
					config_type,
					("\n" + warning) if warning else ""
				)
			)
		


	def __handle_output_dir(self, output_dir, overwrite=False):
		outdir_exists = exists(output_dir)
		if outdir_exists:
			if overwrite:
				print(
					"Output directory already exists and overwrite was requested (-f option).  Deleting directory contents ... ",
					end="", flush=True
				)
				print("DEACTIVATED DUE TO TOO MANY ACCIDENTS.")
				# shutil.rmtree(output_dir)
				# os.makedirs(output_dir)
			else:
				print("Output directory already exists, attempting to resume.", flush=True)
		else:
			pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

		self.logs_dir = join(output_dir, "hpc_logs")
		if not exists(self.logs_dir):
			print("HPC log dir doesn't exist. Creating " + self.logs_dir + " now ... ", end="", flush=True)
			pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
		print("done.")

		self.config_dir = join(output_dir, "config")
		if not exists(self.config_dir):
			print("Config dir does not exist. Creating " + self.config_dir + " now ...", end="", flush=True)
		print("done.")

		self.package_dir = "Data_Package"
		if not exists(self.package_dir):
			print("Package dir does not exist. Creating " + self.package_dir + " now ...", end="", flush=True)
		print("done.")

		print()

		return outdir_exists

	def __make_exe_env_args(self):
		return ExecutionEnvironmentArguments(
			self.scheduler,
			self.partition,
			self.no_drmaa,
			self.max_nodes,
			self.max_cores,
			self.hpc_config
		)

	def __make_exe_env(self):
		print("Configuring execution environment ... ", end="", flush=True)
		self.exe_env = ExecutionEnvironment(
			self.__make_exe_env_args(),
			NOW,
			job_suffix=self.input + "_" + self.output_dir,
			log_dir=self.logs_dir
		)
		print("done.")
		print(str(self.exe_env))
	

	def __init__(self, ap_args):		

		# take all items from argparse args
		for k, v in vars(ap_args).items():
			setattr(self, k, v)

		# check for input
		if not hasattr(self, "input"):
			try:
				self.input = self.input_sheet
			except:
				raise ValueError("Configuration has neither 'input' nor 'input_sheet' attribute.")

		
		

		# make sure output-directory exists and create hpclog-directory  	
		self.__handle_output_dir(self.output_dir, overwrite=self.force)

		# handle hpc configuration
		self.__handle_config_file(
			"hpc_config", 
			"HPC configuration", 
			"hpc_config.json", 
			warning=self.alt_hpc_config_warning if hasattr(self, "alt_hpc_config_warning") else "")

		# hand this over to ExecutionEnvironment
		self.__make_exe_env()

		# handle main main configuration
		self.__handle_config_file(
			"config",
			"configuration",
			"*config.yaml",
			warning=self.alt_config_warning if hasattr(self, "alt_config_warning") else "")

		# Load configuration
		print("Loading configuration from {} ... ".format(self.config_file), end="", flush=True)
		self._config = yaml.load(open(self.config_file))
		print("done.")
		print()




	def __str__(self):
		return super(ConfigurationManager, self).__str__() + "\n" + str(self.exe_env)

	def setConfiguration(self):
		pass



class BGRRLConfigurationManager(ConfigurationManager):
	
	def __manage(self):
		cfg_d = {
			"etc": join(dirname(__file__), "etc"),
			"cwd": os.getcwd(),
			"reapr_correction": False,
			"run_prokka": True,
			"run_ratt": False,
			"package_dir": self.package_dir,			
		}
		self._config.update(cfg_d)
		
		if hasattr(self, "project_prefix"):
			self._config["project_prefix"] = self._config["misc"]["project"] = self.project_prefix if self.project_prefix is not None else ""
		
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

		super(BGRRLConfigurationManager, self).__init__(ap_args)

		self.__manage()




"""
no_normalization        False
no_packaging    False
full_qaa_analysis       False
minimum_survey_assembly_size    1000000.0
input_sheet     samplesheet.csv.20
output_dir      Analysis
project_prefix  bgrrl_test
config  Analysis/config/bgrrl_config.yaml
report_only     False
force   False
enterobase_groups       None
unlock  False
runmode survey
"""	
		
