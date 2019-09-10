import os
import sys
import glob
import re
import argparse


def main(argv=sys.argv[1:]):

    ap = argparse.ArgumentParser()
    ap.add_argument("readpath", type=str)
    ap.add_argument("--lane", "-l", type=str, default="")
	ap.add_argument("--single-cell", action="store_true")
    args = ap.parse_args()

    lane = ("_" + args.lane) if args.lane else ""

    with open("samplesheet.csv", "w") as samplesheet:
        for f in sorted(glob.glob(os.path.join(args.readpath, "*{}_R1.fastq.gz".format(lane)))):
            s = os.path.basename(f).replace("{}_R1.fastq.gz".format(lane), "")
            print(s, s, f, f.replace("R1.fastq.gz", "R2.fastq.gz") if not args.single_cell else "", ",,,,,", sep=",", file=samplesheet)
    return True



if __name__ == "__main__":
    main()
