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
			init_file = glob.glob(join(self.config_dir, file_pattern))[0]
		except:
			init_file = ""

		if hasattr(self, config_file) and exists(getattr(self, config_file)):
			print("Custom {} file specified {}, overriding defaults.".format(config_type, getattr(self, config_file)))
			setattr(self, config_file + "_file", getattr(self, config_file))
		elif init_file:
			print("Found {} file at init location, using this.".format(config_type))
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
		#if not exists(self.logs_dir):
		#	print("HPC log dir doesn't exist. Creating " + self.logs_dir + " now ... ", end="", flush=True)
		#	pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
		#	print("done.")

		#self.config_dir = join(output_dir, "config")
		#if not exists(self.config_dir):
		#	print("Config dir does not exist. Creating " + self.config_dir + " now ...", end="", flush=True)
		#	pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
		#	print("done.")

		#self.report_dir = join(output_dir, "reports")
		#if not exists(self.config_dir):
		#	print("Report dir does not exist. Creating " + self.report_dir + " now ...", end="", flush=True)
		#	pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
		#	print("done.")

		#self.package_dir = "Data_Package"
		#if not exists(self.package_dir):
		#	print("Package dir does not exist. Creating " + self.package_dir + " now ...", end="", flush=True)
		#	pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
		#	print("done.")

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


	def generate_config_file(self, module):
		config_file = join(self.config_dir, module + ".conf.yaml")

		#if not exists(dirname(config_file)):
		#	print("Could not find config-dir, creating ... ", end="", flush=True)
		#	pathlib.Path(dirname(config_file)).mkdir(exist_ok=True, parents=True)
		#	print("done.")	

		with open(config_file, "wt") as cfg_out:
			config_d = OrderedDict(self._config)
			for k, v in sorted(vars(self).items()):
				if not k in self._config and k != "_config":
					print("WRITING {}: {} -> CONFIG".format(k, v))
					config_d[k] = v
			print("Writing configuration to file {} ... ".format(config_file), end="", flush=True)
			yaml.dump(config_d, cfg_out, default_flow_style=False)
			print("done.")

		return config_file
	

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

		# handle main configuration
		self.__handle_config_file(
			"config",
			"configuration",
			"bgrrl_config.yaml",
			warning=self.alt_config_warning if hasattr(self, "alt_config_warning") else "")

		# Load/edit configuration
		print("Loading configuration from {} ... ".format(self.config_file), end="", flush=True)
		self._config = yaml.load(open(self.config_file))
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



class BGRRLConfigurationManager(ConfigurationManager):
	
	def __manage(self):
		cfg_d = {
			"etc": join(dirname(__file__), "..", "etc"),
			"cwd": os.getcwd(),
			"reapr_correction": False,
			"run_prokka": hasattr(self, "run_annotation") and self.run_annotation,
			"run_ratt": False,
			"package_dir": self.package_dir,			
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

		super(BGRRLConfigurationManager, self).__init__(ap_args)

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
			multiqc_dir=join(self.report_dir, "multiqc", stage.split(",")[0]) 
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
