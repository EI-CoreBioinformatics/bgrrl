import pkg_resources

__title__ = "ei:bact-grrl"
__author__ = "Christian Schudoma"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2017 Earlham Institute"
__version__ = "0.1" # pkg_resources.require("eipap")[0].version

import datetime
import os
import sys
import time
import urllib
from enum import Enum, unique
from io import StringIO
from snakemake import snakemake

import requests

NOW = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

DEFAULT_HPC_CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "etc", "hpc_config.json")
DEFAULT_BGRRL_CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "etc", "bgrrl_config.yaml")

PAP_CONFIG = None

@unique
class RunMode(Enum):
	NORMAL = 0
	RESTART = 1
	RESUME = 2
	VALIDATE = 3
	DRYRUN = 4


@unique
class PipelineStep(Enum):
	READ_QC = 0
	ASSEMBLY = 1
	ANNOTATION = 2
	DATA_QA = 3
	PACKAGING = 4
	ATTEMPT_FULL = 5


class Capturing(list):
	def __enter__(self):
		self._stdout = sys.stdout
		self._stderr = sys.stderr
		sys.stdout = self._stringioout = StringIO()
		sys.stderr = self._stringioerr = StringIO()
		return self

	def __exit__(self, *args):
		self.extend(self._stringioout.getvalue().splitlines())
		self.extend(self._stringioerr.getvalue().splitlines())
		del self._stringioout  # free up some memory
		del self._stringioerr  # free up some memory
		sys.stdout = self._stdout
		sys.stderr = self._stderr

def make_exeenv_arg_group(parser):
	hpc_group = parser.add_argument_group("HPC Options",
	                                      "Controls for how jobs should behave across the HPC resources.")

	hpc_group.add_argument("--partition", type=str,
	                       help="Will run all child jobs on this partition/queue, this setting overrides anything specified in the \"--hpc_config\" file.")
	hpc_group.add_argument("--scheduler", type=str,
	                       help="The job scheduler to use.  LSF, PBS and SLURM currently supported.  If running without a scheduler type NONE here. Assumes SLURM by default.")
	hpc_group.add_argument("--no_drmaa", action='store_true', default=False,
	                       help="Use this flag if DRMAA is not available")
	hpc_group.add_argument("-N", "--max_nodes", type=int, default=100,
	                       help="Maximum number of nodes to use concurrently")
	hpc_group.add_argument("-c", "--max_cores", type=int, default=200,
	                       help="Maximum number of cores to use concurrently")
	hpc_group.add_argument("--hpc_config", default=DEFAULT_HPC_CONFIG_FILE,
	                       help="Configuration file for the HPC.  Can be used to override what resources and partitions each job uses.")

class ExecutionEnvironment:
	def __init__(self, args=None, now="eipap", job_suffix=None, log_dir="logs"):
		self.use_drmaa = False
		self.use_scheduler = False
		self.partition = ""
		self.max_nodes = 1
		self.max_cores = 4
		self.sub_cmd = ""
		self.res_cmd = ""
		self.hpc_config = ""
		if args:
			scheduler = args.scheduler if args.scheduler and args.scheduler != '' else "SLURM"
			self.partition = args.partition if args.partition else "{cluster.partition}"
			self.hpc_config = args.hpc_config
			self.use_drmaa = not args.no_drmaa
			self.max_nodes = args.max_nodes
			self.max_cores = args.max_cores
			log_prefix = os.path.join(log_dir, now + "_{rule}_%j")
			job_name = "{rule}_" + job_suffix if job_suffix and job_suffix != "" else "{rule}"
			if scheduler.upper() == "LSF":
				self.sub_cmd = "bsub"
				self.res_cmd = " -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads} -q " + self.partition + " -J " + job_name + " -oo " + log_prefix + ".lsf.log"
				self.use_scheduler = True
			elif scheduler.upper() == "PBS":
				self.sub_cmd = "qsub"
				self.res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q " + self.partition + " -N " + job_name + " -o " + log_prefix + ".pbs.stdout -e " + log_prefix + ".pbs.stderr"
				self.use_scheduler = True
			elif scheduler.upper() == "SLURM":
				self.sub_cmd = "sbatch"
				self.res_cmd = 	" -c {threads}" + \
								" -p " + self.partition + \
								" --exclude={cluster.exclude}" + \
								" --mem={cluster.memory}" + \
								" -J " + job_name + \
								" -o " + log_prefix + ".slurm.log" + \
								" --time={cluster.time}"
				self.use_scheduler = True
			elif scheduler == "" or scheduler.upper() == "NONE":
				pass
			else:
				raise ValueError("Unexpected scheduler configuration.  Check settings.")

	def __str__(self):
		return '\n'.join([
			"Use scheduler: " + str(self.use_scheduler),
			"Use DRMAA: " + str(self.use_drmaa),
			"Submission command: " + self.sub_cmd,
			"Resource command: " + self.res_cmd,
			"Partition: " + self.partition,
			"Max nodes: " + str(self.max_nodes),
			"Max cores: " + str(self.max_cores),
			"HPC configuration file: " + self.hpc_config
		])

def loadPreCmd(command):
	'''
	Used to prefix a shell command that utilises some external software with another command used to load that software
	'''
	if command:
		cc = command.strip()
		if cc != "":
			return "set +u && {} &&".format(cc)

	return ""

def run_snakemake(snakefile, out_dir, cfg_file, exe_env, dryrun=False, unlock=False):
	res = False
	if dryrun:
		print("Dry run requested.  Will not execute tasks.")
		print()
		with Capturing() as output:
			res = snakemake(snakefile,
			                cores=exe_env.max_cores,
			                nodes=exe_env.max_nodes,
			                configfile=cfg_file,
			                workdir=".",
			                unlock=unlock,
			                # force_incomplete=args.force_incomplete,
			                # detailed_summary=args.detailed_summary,
			                # list_resources=args.list_resources,
			                latency_wait=60 if exe_env.use_scheduler else 1,
			                printdag=dryrun,
			                dryrun=False,
			                forceall=dryrun,
			                # allowed_rules=args.allowed_rules
			                )
		if res:
			dag_file = os.path.join(out_dir,
			                        "eibgrrl-" + os.path.basename(snakefile).split('.')[0] + "_dag-" + NOW + ".dot")
			print("Saving DAG of pipeline to " + dag_file + " ... ", end="", flush=True)
			with open(dag_file, 'w') as df:
				print("\n".join(output), file=df)
			print("done.")
		else:
			print("Error occured processing EI-BGRRL:\n" + "\n".join(output))

	else:

		cluster_cfg = exe_env.hpc_config if exe_env.use_scheduler else None
		cluster = (
			exe_env.sub_cmd + exe_env.res_cmd if not exe_env.use_drmaa else None) if exe_env.use_scheduler else None
		drmaa = exe_env.res_cmd if exe_env.use_drmaa else None
		res = snakemake(snakefile,
		                cores=exe_env.max_cores,
		                local_cores=exe_env.max_cores,
		                nodes=exe_env.max_nodes,
		                configfile=cfg_file,
		                workdir=".",
		                cluster_config=cluster_cfg,
		                cluster=cluster,
		                drmaa=drmaa,
		                unlock=unlock,
		                printshellcmds=True,
		                printreason=True,
		                stats=os.path.join(out_dir, os.path.basename(snakefile) + "-" + NOW + ".stats"),
		                jobname="bgrrl.{rulename}.{jobid}",
		                force_incomplete=True,
		                # detailed_summary=args.detailed_summary,
		                # list_resources=True,
		                latency_wait=60 if exe_env.use_scheduler or exe_env.use_drmaa else 1,
		                printdag=dryrun,
		                dryrun=dryrun,
		                forceall=dryrun,
		                verbose=True
		                # allowed_rules=args.allowed_rules
		                )
	return res
