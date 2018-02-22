#!/usr/bin/env python
import sys
import os
from os.path import exists, join, basename, dirname
import glob
import re
import csv
import argparse
import pathlib 

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
    def checkReadFile(_file):
        return os.stat(_file).st_size if exists(_file) else "NA"

    test = "FASTQC:READCOUNT"
    fastqc_dir = join(QCDIR, "fastqc", "bbnorm", sample)
    fastqc_r1 = join(fastqc_dir, sample + "_R1.bbnorm_fastqc", "fastqc_data.txt")
    fastqc_r2 = join(fastqc_dir, sample + "_R2.bbnorm_fastqc", "fastqc_data.txt")
    bbnorm_dir = join(QCDIR, "bbnorm", sample)
    if not (exists(fastqc_r1) and exists(fastqc_r2)):
        # print(join(bbnorm_dir, sample + "_R1.bbnorm.fastq.gz"), file=sys.stderr)
        R1_size = checkReadFile(join(bbnorm_dir, sample + "_R1.bbnorm.fastq.gz"))
        R2_size = checkReadFile(join(bbnorm_dir, sample + "_R2.bbnorm.fastq.gz"))
        if R1_size > 20 or R2_size > 20:
            print("WARNING: One or both FASTQC-report(s) missing for sample {}, but read files seem not to be empty: R1_fastqc={},R1_read_file_size={} R2_fastqc={},R2_read_file_size={}.".format(sample, exists(fastqc_r1), R1_size, exists(fastqc_r2), R2_size), file=sys.stderr)
        return TestResult(test, "FAIL", "MISSING", (0, 0)) 
    n1, n2 = map(extractReadCount, (fastqc_r1, fastqc_r2))
    if n1 != n2: 
        return TestResult(test, "FAIL", "INCONSISTENT", (n1, n2))
    if n1 < min_reads:
        return TestResult(test, "FAIL", "LOW", (n1, n2))
    return TestResult(test, "PASS", "", (n1, n2)) 

def test_kat_peaks(sample, max_peak_volume_threshold=0.9):
    test = "KAT:PEAKS"
    kat_log = join(QCDIR, "kat", sample, sample + ".kat")
    status, errmsg, data = "PASS", "", tuple()
    kmer_peaks, gc_peaks, gsize, unit, kmer_freq, mean_gc = None, None, None, None, None, None
    if not exists(kat_log) or os.stat(kat_log).st_size == 0:
        stats, errmsg = "FAIL", "MISSING"
    else:
        with open(kat_log) as _in:
            kmer_peak_table, state, read_table = list(), "", False
            for line in _in:
                if line.startswith("K-mer frequency spectra statistics"):
                    state = "kmer_freq"
                    try:
                        [next(_in), next(_in), next(_in), next(_in)]
                    except:
                        print("PROBLEM:", kat_log, "truncated in K-mer frequency section", file=sys.stderr)
                        sys.exit(1)
                    read_table = True
                elif line.startswith("Calculating genome statistics"):
                    state = "genome_stats"
                elif line.startswith("GC distribution statistics"):
                    state = "gc_dist"
                else:
                    if state == "kmer_freq":
                        if line.startswith("Peaks in analysis: "):
                            try:
                                kmer_peaks = int(line.replace("Peaks in analysis: ", "").strip())
                            except:
                                kmer_peaks = None
                        elif line.startswith("Overall mean k-mer frequency: "):
                            kmer_freq = line.replace("Overall mean k-mer frequency: ", "").strip()
                        elif not line.strip() or line.startswith("K-value used: "):
                            read_table = False
                        elif read_table:
                            try:
                                kmer_peak_table.append(list(map(float, re.sub(" +", " ", line.strip()).split(" "))))
                            except:
                                kmer_peak_table.append(list(map(float, re.sub(" +", " ", line.strip()).split(" ")[:-1])))
                                # print("Error with kmer_peak_table in", kat_log, file=sys.stderr)
                                # sys.exit(1)
                                # volume_col = 6
                    elif state == "gc_dist":
                        try:
                            gc_peaks = int(line.replace("Peaks in analysis: ", "").strip())
                        except:
                            gc_peaks = None
                    elif state == "genome_stats" and line.startswith("Estimated genome size: "):
                        line = line.replace("Estimated genome size: ", "").strip()
                        try:
                            gsize, unit = line.split(" ")
                        except:
                            print("XXX", file=sys.stderr)
                            gsize, unit = line, "bp"
                        gsize = float(gsize)
                    elif state == "gc_dist" and line.startswith("Mean GC: "):
                        mean_gc = line.replace("Mean GC: ", "").strip()

    

    try:
        max_peak = sorted(kmer_peak_table, key=lambda x:x[6])[-1]
    except: 
        max_peak = [None]*7
    try:
        total_volume = sum(row[6] for row in kmer_peak_table)
    except:
        print(kmer_peak_table, file=sys.stderr)
        total_volume = 0
    data = (kmer_peaks, max_peak[6], max_peak[6] / total_volume if (max_peak[6] is not None and total_volume) else None,  gc_peaks, gsize, unit, kmer_freq, mean_gc)

    if data[0] is None:
        status, errmsg = "PASS", "MISSING"
    elif data[0] < 1:
        status, errmsg = "PASS", "NOPEAK"
    elif data[0] > 1:
        if max_peak[6] / total_volume < max_peak_volume_threshold:
            status, errmsg = "PASS", "MULTI_MODAL"

    return TestResult(test, status, errmsg, data)

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
    quast_report = join(QADIR, "quast", sample, "report.tsv")

    status, errmsg, assembly_size = "PASS", "", 0
    if exists(quast_report):
        assembly_size = extractAssemblySize(open(quast_report))
        if assembly_size < min_size:
            status, errmsg = "FAIL", "TOO_SMALL"
    else:
        status, errmsg = "FAIL", "NOT_ASSEMBLED"
            
    return TestResult(test, status, errmsg, str(assembly_size))

TESTS = [("FASTQC:READCOUNT", test_fastqc_readcount),
         ("KAT:PEAKS", test_kat_peaks),
         # ("KAT:HIST", test_kat_hist),
         # ("KAT:GCP", test_kat_gcp),
         ("TADPOLE:SIZE", test_tadpole_size)]


def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=str)
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("--min_readcount", type=int, default=1000)
    ap.add_argument("--min_tadpolesize", type=int, default=1e6)
    ap.add_argument("--report-dir", type=str, default="")
    args = ap.parse_args(args_in)
    
    print("Running qc:evaluation...", end="", flush=True)

    if not args.report_dir:
        report_dir = os.path.join(args.indir, "reports")
    else:
        report_dir = os.path.join(args.report_dir)
    pathlib.Path(os.path.join(report_dir, "samplesheets")).mkdir(parents=True, exist_ok=True)
    

    INPUTFILES = dict(readSamplesheet(args.input))
    global QCDIR
    # QCDIR  = join(args.indir, "Analysis", "qc")
    QCDIR = join(args.indir, "qc")
    global QADIR
    QADIR = join(args.indir, "qa", "survey")

    #Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
    
    qc_eval_outf, asm_samplesheet_f = join(report_dir, "qc_eval.tsv"), join(report_dir, "samplesheets", "samplesheet.qc_pass.tsv")
    with open(qc_eval_outf, "w") as qc_eval_out, open(asm_samplesheet_f, "w") as asm_samplesheet:
        print("SAMPLE", *tuple(test for test, testf in TESTS), sep="\t", file=qc_eval_out)
        for sample in sorted(INPUTFILES):
            results = list(testf(sample) for test, testf in TESTS)
            if all(result[1] == "PASS" for result in results):
                print(*INPUTFILES[sample], sep=",", file=asm_samplesheet)
            elif results[0][1] == "FAIL" and results[0][2] == "MISSING" and results[2][1] == "PASS":
                print("WARNING: sample {} has missing fastqc report(s), but passes survey assembly.".format(sample), file=sys.stderr)

            # p_results = tuple(",".join(result) for result in results)
            # print(p_results)

            p_result = tuple(",".join(map(str, result[1:])) for result in results)
          
            print(sample, *p_result, sep="\t", file=qc_eval_out)


    print(" Done.\n Generated qc_eval report in {} and asm-samplesheet in {}.".format(qc_eval_outf, asm_samplesheet_f))

    return True

if __name__ == "__main__":
    main()
