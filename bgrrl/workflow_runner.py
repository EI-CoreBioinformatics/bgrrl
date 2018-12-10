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


class WorkflowRunner(object):
	def __init(self):
		pass


#	def __handle_output_dir(self, output_dir, overwrite=False):
#		outdir_exists = exists(output_dir)
#		if outdir_exists:
#			if overwrite:
#				print(
#					"Output directory already exists and overwrite was requested (-f option).  Deleting directory contents ... ",
#					end="", flush=True
#				)
#				print("DEACTIVATED DUE TO TOO MANY ACCIDENTS.")
#				# shutil.rmtree(output_dir)
#				# os.makedirs(output_dir)
#			else:
#				print("Output already exists, attempting to resume.", flush=True)
#		else:
#			pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
#
#		self.logs_dir = join(output_dir, "hpc_logs")
#		if not exists(self.logs_dir) and self.exe_env.use_scheduler:
#			print("HPC log dir doesn't exist.  Creating " + self.logs_dir + " now ... ", end="", flush=True)
#			pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)
#
#		print("done.")
#		print()
#
#		return outdir_exists
#
#	def __init__(self):
#		
#		if not hasattr(self.args, "input"):
#			try:
#				self.args.input = args.input_sheet
#			except:
#				raise ValueError("args does neither have input nor input_sheet attribute.")
#
#		# Establish a valid cluster configuration... may throw if invalid
#		print("Configuring execution environment ... ", end="", flush=True)
#		self.exe_env = ExecutionEnvironment(
#			self.args.exe_env_args, 
#			NOW, 
#			job_suffix=self.args.input_sheet + "_" + self.args.output_dir, 
#			log_dir=join(self.args.output_dir, "hpc_logs")
#		)
#		print("done.")
#		print(str(self.exe_env))
#
##		### make sure output-directory exists and create hpclog-directory
##		### outdir_exists = self.__handle_output_dir(self.args.output_dir, overwrite=self.args.force)
##
##		# see if there are config files at the init location 
##		try:
##			hpc_config_init = glob.glob(join(self.args.output_dir, "config", "hpc_config.json"))[0]
##		except:
##			hpc_config_init = ""
##		try:
##			config_init = glob.glob(join(self.args.output_dir, "config", "*config.yaml"))[0]
##		except:
##			config_init = ""
##		
##
##		if self.args.hpc_config and exists(self.args.hpc_config):
##			print("Custom HPC configuration file specified, overriding defaults")
##			self.hpc_config_file = self.args.hpc_config
##			self.args.hpc_config_file = self.args.hpc_config
##		elif hpc_config_init:
##			print("Found HPC configuration at init location ({}), using this.".format(hpc_config_init))
##			self.hpc_config_file = hpc_config_init
##			self.args.hpc_config_file = hpc_config_init
##		else:	   
##			raise ValueError("No valid HPC configuration specified ({}). {}".format(self.args.hpc_config, self.args.alt_hpc_config_warning if hasattr(self.args, "alt_hpc_config_warning") else ""))
##
##		print()
##
##		if self.args.config and exists(self.args.config):
##			print("Custom configuration file specified, overriding defaults")
##			self.config_file = self.args.config
##			self.args.config_file = self.args.config
##		elif config_init:
##			print("Found configuration at init location ({}), using this.".format(config_init))
##			self.config_file = config_init
##			self.args.config_file = config_init
##		else:
##			raise ValueError("No valid configuration specified ({}). {}".format(self.args.config, self.args.alt_config_warning if hasattr(self.args, "alt_config_warning") else ""))
##
##		print("Loading configuration from {} ... ".format(self.config_file), end="", flush=True)
##		self.config = yaml.load(open(self.config_file))
##		print("done.")
##		print()
##
##
