from collections import OrderedDict, namedtuple
import os 
from os.path import exists, dirname, basename, join
import sys
import yaml
import pathlib
import glob

from eicore.snakemake_helper import *

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
			init_file = glob.glob(join(self.config_dir, file_pattern))[0]
		except:
			init_file = ""

		if hasattr(self, config_file) and exists(getattr(self, config_file)):
			print("Custom {} file specified {}.".format(config_type, getattr(self, config_file)))
			setattr(self, config_file + "_file", getattr(self, config_file))
		elif init_file:
			print("Found {} file at default location, using this.".format(config_type))
			setattr(self, config_file + "_file", init_file)
			setattr(self, config_file, init_file)
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


		def _create_subdir(subdir, name):		
			if not exists(subdir):
				print(name + " does not exist. Creating " + subdir + " now ... ", end="", flush=True)
				pathlib.Path(subdir).mkdir(parents=True, exist_ok=True)
				print("done.")
			return subdir

		self.logs_dir = _create_subdir(join(output_dir, "hpc_logs"), "HPC log dir")
		self.config_dir = _create_subdir(join(output_dir, "config"), "Config dir")
		self.report_dir = _create_subdir(join(output_dir, "reports"), "Report dir")
		self.package_dir = _create_subdir("Data_Package", "Package dir")
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
			job_suffix=self.job_suffix_prefix + "_" + self.output_dir,
			log_dir=self.logs_dir
		)
		print("done.")
		print(str(self.exe_env))


	def generate_config_file(self, module):
		config_file = join(self.config_dir, module + ".conf.yaml")

		with open(config_file, "wt") as cfg_out:
			config_d = OrderedDict(self._config)
			for k, v in sorted(vars(self).items()):
				if not k in self._config and k != "_config":
					config_d[k] = v
			print("Writing configuration to file {} ... ".format(config_file), end="", flush=True)
			yaml.dump(config_d, cfg_out, default_flow_style=False)
			print("done.")

		return config_file
	

	def __init__(self, ap_args):		

		for k, v in vars(ap_args).items():
			setattr(self, k, v)

		if not hasattr(self, "config"):
			raise ValueError("Cannot establish configuration as config argument is not present.")

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

		# handle main configuration
		self.__handle_config_file(
			"config",
			"configuration",
			#"bgrrl_config.yaml",
			self.config,
			warning=self.alt_config_warning if hasattr(self, "alt_config_warning") else "")

		# Load/edit configuration
		print("Loading configuration from {} ... ".format(self.config_file), end="", flush=True)
		self._config = yaml.load(open(self.config_file), Loader=yaml.SafeLoader)
		print("done.")
		print()

		self._config["out_dir"] = self.output_dir

		# get multiqc configuration for qaa
		self.__handle_config_file(
			"multiqc_config",
			"MultiQC configuration",
			"multiqc_config.yaml",
			warning=self.alt_multiqc_config_warning if hasattr(self, "alt_multiqc_config_warning") else "")
			
	def __str__(self):
		return super(ConfigurationManager, self).__str__() + "\n" + str(self.exe_env)

	def setConfiguration(self):
		pass
