import os
import sys
import argparse


def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("--report-dir", type=str, default="")
    args = ap.parse_args(args_in)


    for




if __name__ == "__main__":
    main()
