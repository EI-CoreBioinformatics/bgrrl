import sys
import csv
import os
from os.path import join, dirname, basename
import subprocess
import pathlib

from collections import namedtuple, OrderedDict, Counter
from argparse import ArgumentParser


GFFeature = namedtuple("GFFeature", "seqname source feature start end score strand frame attribute".split(" "))

PROKKA_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA"]

HEADER = ["Sample"] + PROKKA_FEATURES
# HEADER.extend(map(lambda s:s + "(Prokka)", PROKKA_FEATURES))

def readProkkaGFF(_in):
    data = OrderedDict() 
    for row in csv.reader(_in, delimiter="\t"):
        if row and not row[0].startswith("#"):
            attr = OrderedDict((item.split("=")[0], item.split("=")[1]) for item in row[8].split(";"))
            key = attr.get("ID", None)
            if key is not None:
                # print(row, file=sys.stderr)
                data[key] = GFFeature(row[0], row[1], row[2], int(row[3]), int(row[4]), row[5], row[6], row[7], attr)
    return data


def main(args=sys.argv[1:]):
    ap = ArgumentParser()
    ap.add_argument("prokka_dir", type=str)
    ap.add_argument("report_dir", type=str)
    args = ap.parse_args(args)

    delta_dir = join(args.prokka_dir, "..", "delta")
    pathlib.Path(delta_dir).mkdir(parents=True, exist_ok=True) 

    samples = next(os.walk(args.prokka_dir))[1]

    with open(join(args.report_dir, "annotation_report.tsv"), "w") as prokka_out:
        print(*HEADER, sep="\t", file=prokka_out)
        for sample in samples:
            prokka_gff = join(args.prokka_dir, sample, sample + ".prokka.gff")
            with open(prokka_gff) as _in:
                prokka_full = readProkkaGFF(_in)
                prokka_fcount = Counter(feature.feature for feature in prokka_full.values())    


            row = [sample]
            row.extend(prokka_fcount[feature] for feature in PROKKA_FEATURES)            

            print(*row, sep="\t", file=prokka_out)

    return True

if __name__ == "__main__":
    main()
