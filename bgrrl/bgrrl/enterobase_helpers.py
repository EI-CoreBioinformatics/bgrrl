#!/usr/bin/env python
import sys
import os
import glob
import csv
import argparse
from collections import namedtuple, Counter
from os.path import exists
import yaml

ECriteria = namedtuple("ECriteria", "minsize maxsize n50 ncontigs ncount spcount".split(" "))

def loadEnterobaseCriteria(criteria_file):
    with open(criteria_file) as crit_in:
        d = yaml.load(crit_in)
    eb_crit = dict()
    for k in d:                
        eb_crit[k] = ECriteria(int(d[k]["min_genome_size"]),
                               int(d[k]["max_genome_size"]),
                               int(d[k]["n50"]),
                               int(d[k]["ncontigs"]),
                               float(d[k]["N_fraction"]),
                               float(d[k]["genus_fraction"]))
    return eb_crit

                     



def validateEnterobaseInput(eb_input_str, eb_criteria):
    print("Determining downstream processing mode...", end="", flush=True)
    if not eb_input_str:
        print(" DEFAULT")
        valid_sample_groups = list()
    else:
        if eb_input_str == "all":
            print(" EB:all")
            check_sample_groups = list(eb_criteria.keys())
        else:
            check_sample_groups = eb_input_str.strip().split(",")
       
        valid_sample_groups = list(filter(lambda g:g in eb_criteria, check_sample_groups))
        invalid_sample_groups = set(check_sample_groups).difference(valid_sample_groups)

        if not valid_sample_groups:
            print(" DEFAULT:FALLBACK")
            print("No valid Enterobase groups found. Proceeding without checking Enterobase criteria.")
            pass
        elif not eb_input_str == "all":
            print(" EB")
            print("Moving forward with Enterobase groups: {}.".format(",".join(sorted(valid_sample_groups))))
        if invalid_sample_groups:
            print("Found invalid Enterobase groups: {}.".format(",".join(sorted(invalid_sample_groups))))

    return valid_sample_groups

