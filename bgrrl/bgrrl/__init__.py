import pkg_resources

__title__ = "ei:bact-grrl"
__author__ = "Christian Schudoma"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2017 Earlham Institute"
__version__ = pkg_resources.require("bgrrl")[0].version

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

TIME_CMD = " /usr/bin/time -v"

@unique
class PipelineStep(Enum):
	READ_QC = 0
	ASSEMBLY = 1
	ANNOTATION = 2
	DATA_QA = 3
	FINALIZE = 4
	ATTEMPT_FULL = 5

"""
def loadPreCmd(command, is_dependency=True):
	'''
	Used to prefix a shell command that utilises some external software with another command used to load that software
	'''
	if command:
		cc = command.strip()
		if cc != "":
			if is_dependency:
				return "set +u && {} &&".format(cc)
			else:
				return " {} ".format(cc)

	return ""
"""
