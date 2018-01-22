#!/usr/bin/env python
import sys
import os
import glob
import csv
import argparse
from collections import namedtuple
from os.path import exists
import yaml

from bgrrl import run_snakemake

Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))

def readSamplesheet(fn, delimiter=","):
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



def run_qc(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config # dict()
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    
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

def run_asm(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
    config["cwd"] = os.getcwd()
    # print(args)
    
    if args.contig_minlen:
        config["use_asm_lengthfilter"] = True
        config["asm_lengthfilter_contig_minlen"] = args.contig_minlen
        # print(config)

    config_file = os.path.join(out_dir, "bg-asm.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-ASM")
    # print("THIS:", __file__)
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-asm.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock)
    return res



if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('samplesheet', type=str)

    args = ap.parse_args()

    assert args.samplesheet and os.path.exists(args.samplesheet)
    samples = dict(readSamplesheet(args.samplesheet))

    for sample in samples:
        print(*samples[sample])
