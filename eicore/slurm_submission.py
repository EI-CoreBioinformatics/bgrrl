#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import subprocess
from abc import ABC, abstractmethod

from eipap import create_eipap_argparser, DEFAULT_PAP_CONFIG_FILE, NOW
from eipap.illumina.transfer import create_eitransfer_argparser
from eipap.illumina.analysis import create_eianalysis_argparser

DEFAULT_PARTITION = "tgac-pap"
DEFAULT_EMAIL = "$(whoami)@nbi.ac.uk"

from contextlib import contextmanager

@contextmanager
# https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)



class SlurmSubmitter(ABC):

	@staticmethod
	def get_resume_string(output_dir, prefix):
		r = 1
		resume_infix = ".resume" if prefix else ""
		file_prefix = output_dir + "." + prefix + resume_infix 
		while os.path.exists(file_prefix + str(r) + ".out.log") or os.path.exists(file_prefix + str(r) + ".err.log"):
			r += 1
		
		return resume_infix + str(r)

	def __init__(self, config=):
		if not hasattr(self, "whoami"):
			raise ValueError("Error: Cannot run without specifying module.")
		config = 
		
class EISlurmSubmitter(SlurmSubmitter):
	

	def __init__(self): 

		pap_config = yaml.load(open(DEFAULT_PAP_CONFIG_FILE), Loader=yaml.SafeLoader)
		self.argparser = self.argparser_func(pap_config)
		self.argparser.add_argument("--no-submit", action="store_true")	
		args = self.argparser.parse_args()

		self.__compile_eipap_options(args)

		self.qc_dir = None
		if self.jira is not None:
			try:
				p = subprocess.Popen("eipap_cd qc {}".format(self.jira), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				self.qc_dir = p.communicate()[0].decode().strip().split("\n")[-1]
			except:
				raise ValueError("Error: Cannot determine qc directory for {}".format(self.jira))
			self.rundir = self.qc_dir
		elif self.whoami == "eipap":
			self.rundir = "." # args.run_dir if hasattr(args, "run_dir") else None
		elif self.whoami == "eianalysis":
			self.rundir = "." # args.input
		else:
			raise ValueError("Error: Looks like someone is trying to run eitransfer without jira. Impossible.")

		default_output_dir = "PAP_" + NOW.split("_")[0]
		if args.output_dir:
			self.output_dir = args.output_dir
			# This will be treated as local, but we cd into the PAP qc-dir anyway
			self.submission_script = os.path.basename(self.output_dir) + self.submission_script_suffix
		elif self.whoami != "eitransfer":
			self.output_dir = default_output_dir
			# This will be treated as local, but we cd into the PAP qc-dir anyway
			self.submission_script = os.path.basename(self.output_dir) + self.submission_script_suffix
		else:
			self.output_dir = ""
			self.submission_script = os.path.basename(args.input) + self.submission_script_suffix

		self.__compile_slurm_options(args)


	def __remove_submission_script(self):
			try:
				os.remove(self.submission_script)
			except:
				pass

	def __write_submission_script(self):
		# adding output_dir -o, as it was ignored in eipap_option parsing: easier to deal with custom output_dir 
		od_string = ("-o " + self.output_dir) if self.output_dir else ""
		try:	
			with open(self.submission_script, "w") as cmd_out:
				print("#!/bin/bash -e", file=cmd_out)
				print(*([self.whoami, od_string] + self.eipap_options), file=cmd_out)
		except:
			raise ValueError("Could not write submission script.")

	def __compile_eipap_options(self, args, ignore={"partition", "input", "output_dir", "no_submit"}):
		self.eipap_options = list()
		for k, v in vars(args).items():
			if k not in ignore and v is not None and v:			
				opt = "--" + k                                                  	
				if not type(v) is bool:
					opt += "=" + str(v)
				self.eipap_options.append(opt)			

		# eipap: jira = input, if input == "0" (jira-less): jira -> None
		# eianalysis: jira = args.jira (optional), if not args.jira: jira -> None
		# eitransfer: jira = args.jira (mandatory)
		try:
			self.jira = args.jira
		except:
			self.jira = args.input if (self.whoami == "eipap" and args.input != "0") else None

		self.no_submit = args.no_submit

		self.eipap_options.append(args.input)
		print(*self.eipap_options)

	def __compile_slurm_options(self, args):

		email = os.environ.get("SNAKEMAKE_EMAIL", DEFAULT_EMAIL)

		resume_str = get_resume_string(
			self.output_dir if self.whoami != "eitransfer" else args.input, 
			self.whoami.replace("ei", "").replace("transfer", "")
		) if ((hasattr(args, "mode") and args.mode.upper() == "RESUME") or self.whoami == "eitransfer") else ""

		partition = args.partition if hasattr(args, "partition") and args.partition else DEFAULT_PARTITION
	
		self.slurm_options = "-J {5}_{0}_{1} -o {1}.{6}{2}.out.log -e {1}.{6}{2}.err.log --mail-type=FAIL --mail-user={3} -p {4}"

		#out_prefix = args.input if self.whoami == "eitransfer" or not self.output_dir

		self.slurm_options = self.slurm_options.format(
			args.input if self.whoami != "eitransfer" else "",
			#self.output_dir if self.output_dir else args.input,
			self.output_dir if self.whoami != "eitransfer" else args.input,
			resume_str,
			email,
			partition,
			self.whoami,
			self.whoami.replace("ei", "")
		)

	def run(self, cleanup=True):

		with cd(self.rundir if self.rundir is not None else "."):
			if self.jira is not None:
				print("Changing into QC directory for {0} ".format(self.jira))
			else:
				print("Changing into rundir {} ".format(self.rundir))

			self.__remove_submission_script()
			self.__write_submission_script()

			print("Running submission command: sbatch {} --wrap=\"{} {}\"".format(
				self.whoami,
				self.slurm_options, 
				" ".join(map(str, self.eipap_options))
				)
			)

			if not self.no_submit:                                                                               
				p = subprocess.Popen(
					"sbatch {} {}".format(self.slurm_options, self.submission_script),
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					shell=True
				)

				o, e = p.communicate()

				print(o.decode())
				if e.decode():
					print(e.decode())

			if cleanup:
				self.__remove_submission_script()



class EipapSubmitter(SlurmSubmitter):
	def __init__(self):
		self.whoami = "eipap"
		self.argparser_func = create_eipap_argparser
		self.submission_script_suffix = ".pap_command.sh"
		super().__init__()

class EitransferSubmitter(SlurmSubmitter):
	def __init__(self):
		self.whoami = "eitransfer"
		self.argparser_func = create_eitransfer_argparser
		self.submission_script_suffix = ".transfer_command.sh"
		super().__init__()
		
		if not hasattr(self, "jira") or not self.jira:
			raise ValueError("Error: eitransfer requires the -j parameter")	

class EianalysisSubmitter(SlurmSubmitter):
	def __init__(self):
		self.whoami = "eianalysis"
		self.argparser_func = create_eianalysis_argparser
		self.submission_script_suffix = ".analysis_command.sh"
		super().__init__()


def eipap_main():
	print("EIPAP_MAIN")
	EipapSubmitter().run()
def eitransfer_main():
	print("EITRANSFER_MAIN")
	EitransferSubmitter().run()
def eianalysis_main():
	print("EIANALYSIS_MAIN")
	EianalysisSubmitter().run()


def main():
	pass
