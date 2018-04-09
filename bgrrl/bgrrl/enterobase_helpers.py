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

# https://bitbucket.org/enterobase/enterobase-web/wiki/EnteroBase%20Backend%20Pipeline%3A%20QA%20evaluation
ENTERO_CRITERIA = { "Salmonella": ECriteria(4000000, 5800000, 20000, 600, 0.03, 0.7),
                    "Escherichia": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Shigella": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Yersinia": ECriteria(3700000, 5500000, 15000, 600, 0.03, 0.65),
                    "Moraxella": ECriteria(1800000, 2600000, 20000, 600, 0.03, 0.65)
 }

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

