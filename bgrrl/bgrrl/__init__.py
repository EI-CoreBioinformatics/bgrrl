import pkg_resources

__title__ = "ei:bgrr|"
__author__ = "Christian Schudoma"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2017-2018 Earlham Institute"
__version__ = pkg_resources.require("bgrrl")[0].version

import datetime
import os
import sys
import time
import urllib
from enum import Enum, unique
from io import StringIO
from snakemake import snakemake
from collections import namedtuple, Counter
import yaml
import requests

from eicore.external_process.snakemake_helper import run_snakemake

NOW = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

DEFAULT_HPC_CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "etc", "hpc_config.json")
DEFAULT_BGRRL_CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "etc", "bgrrl_config.yaml")

TIME_CMD = " /usr/bin/time -v"

@unique
class PipelineStep(Enum):
    READ_QC = 0
    ASSEMBLY = 1
    ANNOTATION = 2
    DATA_QA = 3
    FINALIZE = 4
    ATTEMPT_FULL = 5

Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
"""
def run_qc(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config # dict()
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["no_normalization"] = args.no_normalization

    config_file = os.path.join(out_dir, "bg-qc.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-QC")
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-qc.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock)


    # write samplesheet into configfile
    # test if sample data exists
    # run bg-qc on samples
    # upon return, parse bg-qc's output and annotate samples
    # generate samplesheet for bg-asm
    return res

"""

def readSamplesheet(fn, delimiter=","):
    import csv
    with open(fn) as fi:
        for row in csv.reader(fi, delimiter=delimiter):
            # print(row)
            sample = Sample(*row)
            # print(sample)
            # print((sample.R1 and exists(sample.R1) and sample.R2 and exists(sample.R2)))
            assert sample.sampleID
            # assert (sample.R1 and exists(sample.R1) and sample.R2 and exists(sample.R2)) or (sample.S and exists(sample.S))
            if not sample.customerSampleID:
                row[1] = sample.sampleID
            yield (sample.sampleID, Sample(*row))

def verifySamplesheet(fn, delimiter=","):
    for sample_id, sample in readSamplesheet(fn, delimiter=delimiter):
        r1_exists, r2_exists, s_exists = map(os.path.exists, (sample.R1, sample.R2, sample.S))
        if not (r1_exists and r2_exists):
            raise ValueError("Cannot find R1/R2 data at R1={}, R2={}.".format(sample.R1, sample.R2))
        return True


class BGRRLModuleRunner(object):
    def __init__(self, module, args, exe_env, config=dict()):
        print(config)
        self.config = dict(config)
        self.module = module
        self.outdir = args.output_dir
        self.unlock = args.unlock
        self.exe_env = exe_env


        if verifySamplesheet(args.input):
            self.config["samplesheet"] = args.input
        self.config["out_dir"] = self.outdir

        for k, v in args._get_kwargs():
            self.config[k] = v

        self.config_file = os.path.join(self.outdir, module + ".conf.yaml")
        with open(self.config_file, "w") as conf_out:
            yaml.dump(self.config, conf_out, default_flow_style=False)

    def run(self):
 
        print("Running " + self.module)
        snake = os.path.join(os.path.dirname(__file__), "zzz", self.module + ".smk.py")
        run_result = run_snakemake(snake, self.outdir, self.config_file, self.exe_env, dryrun=False, unlock=self.unlock)

        return run_result


