#!/usr/bin/env python
import sys
import os
import glob
import re
import csv
import argparse

from collections import Counter, namedtuple

from bgrrl.bgrrl import readSamplesheet, Sample

TestResult = namedtuple("TestResult", "test status errmsg data".split(" "))

def test_fastqc_readcount(sample, min_reads=1000):
    def extractReadCount(fn):
        try:
                with open(fn) as fi:
                    for line in fi:
                        if line.startswith("Total Sequences"):
                            return int(line.strip().split("\t")[-1])
        except FileNotFoundError:
            return 0

    test = "FASTQC:READCOUNT"
    fastqc_dir = os.path.join(QCDIR, "fastqc", "bbnorm", sample)
    fastqc_r1 = os.path.join(fastqc_dir, sample + "_R1.bbnorm_fastqc", "fastqc_data.txt")
    fastqc_r2 = os.path.join(fastqc_dir, sample + "_R2.bbnorm_fastqc", "fastqc_data.txt")   
    if not (os.path.exists(fastqc_r1) and os.path.exists(fastqc_r2)):
        return TestResult(test, "FAIL", "MISSING", (0, 0)) 
    n1, n2 = map(extractReadCount, (fastqc_r1, fastqc_r2))
    if n1 != n2: 
        return TestResult(test, "FAIL", "INCONSISTENT", (n1, n2))
    if n1 < min_reads:
        return TestResult(test, "FAIL", "LOW", (n1, n2))
    return TestResult(test, "PASS", "", (n1, n2)) 

def checkKatPeaks(kat_log):
    def extractKatPeaks(fn):
        with open(fn) as _in:
            for line in _in:
                if line.startswith('Peaks in analysis'):
                    return int(line.strip().split(':')[1].strip())
        return 0

    if not os.path.exists(kat_log):
        return ("FAIL", "MISSING", (0,))
    kat_peaks = extractKatPeaks(kat_log)
    if kat_peaks == 1:
        return ("PASS", "", (1,))
    return ("FAIL", "MULTIMODAL" if kat_peaks else "NO_PEAK", (kat_peaks,))


def test_kat_hist(sample):    
    status, errmsg, data = checkKatPeaks(os.path.join(QCLOGDIR, sample + ".qc_kathist.log"))
    return TestResult("KAT:HIST", status, errmsg, data)
def test_kat_gcp(sample):
    status, errmsg, data = checkKatPeaks(os.path.join(QCLOGDIR, sample + ".qc_katgcp.log"))
    return TestResult("KAT:GCP", status, errmsg, data)
   
"""
Assembly	FD01543412_L006
# contigs (>= 0 bp)	1477
# contigs (>= 1000 bp)	944
# contigs (>= 5000 bp)	319
# contigs (>= 10000 bp)	90
# contigs (>= 25000 bp)	3
# contigs (>= 50000 bp)	0
Total length (>= 0 bp)	4727022
Total length (>= 1000 bp)	4504601
Total length (>= 5000 bp)	2908093
Total length (>= 10000 bp)	1314760
Total length (>= 25000 bp)	85557
Total length (>= 50000 bp)	0
# contigs	1477
Largest contig	30650
Total length	4727022
GC (%)	52.11
N50	6426
N75	3429
L50	224
L75	472
# N's per 100 kbp	0.00
"""
 
def test_tadpole_size(sample, min_size=1e6):
    def extractAssemblySize(qreport_in):
        for row in csv.reader(qreport_in, delimiter="\t"):
            if row[0].startswith("Total length"):
                return int(row[1])
        return 0
    test = "TADPOLE:SIZE"
    quast_report = os.path.join(QCDIR, "tadpole", sample, "quast", "report.tsv")   
    if os.path.exists(quast_report) and extractAssemblySize(open(quast_report)) < min_size:
        return TestResult(test, "FAIL", "TOO_SMALL", extractAssemblySize(open(quast_report)))
    return TestResult(test, "PASS", "", str(extractAssemblySize(open(quast_report))))

TESTS = [("FASTQC:READCOUNT", test_fastqc_readcount),
         # ("KAT:HIST", test_kat_hist),
         # ("KAT:GCP", test_kat_gcp),
         ("TADPOLE:SIZE", test_tadpole_size)]


def main(args_in=sys.argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=str)
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("--min_readcount", type=int, default=1000)
    ap.add_argument("--min_tadpolesize", type=int, default=1e6)
    args = ap.parse_args(args_in)
    
    print("Running qc:evaluation...", end="", flush=True)
    

    INPUTFILES = dict(readSamplesheet(args.input))
    global QCDIR
    # QCDIR  = os.path.join(args.indir, "Analysis", "qc")
    QCDIR = os.path.join(args.indir, "qc")

    #Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
    
    qc_eval_outf, asm_samplesheet_f = os.path.join(args.indir, "..", "qc_eval.tsv"), os.path.join(args.indir, "..", "samplesheet.qc_pass.tsv")
    with open(qc_eval_outf, "w") as qc_eval_out, open(asm_samplesheet_f, "w") as asm_samplesheet:
        print("SAMPLE", *tuple(test for test, testf in TESTS), sep="\t", file=qc_eval_out)
        for sample in sorted(INPUTFILES):
            results = list(testf(sample) for test, testf in TESTS)
            if all(result[1] == "PASS" for result in results):
                print(*INPUTFILES[sample], sep=",", file=asm_samplesheet)

            # p_results = tuple(",".join(result) for result in results)
            # print(p_results)

            p_result = tuple(",".join(map(str, result[1:])) for result in results)
          
            print(sample, *p_result, sep="\t", file=qc_eval_out)


    print(" Done.\n Generated qc_eval report in {} and asm-samplesheet in {}.".format(qc_eval_outf, asm_samplesheet_f))



if __name__ == "__main__":
    main()
