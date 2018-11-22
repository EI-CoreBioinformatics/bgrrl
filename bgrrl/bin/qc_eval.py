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

# from bgrrl.samplesheet import readSamplesheet, Sample, ASM_Sample, sample2asmsample
from bgrrl.samplesheet import *

TestResult = namedtuple("TestResult", "test status errmsg data".split(" "))

def test_fastqc_readcount(sample, min_reads=1000, readtype="bbnorm", **kwargs):
    def extractReadCount(fn):
        try:
            with open(fn) as fi:
                for line in fi:
                    if line.startswith("Total Sequences"):
                        return int(line.strip().split("\t")[-1])
        except FileNotFoundError:
            return 0
        return 0
    def extractReadLength(fn):
        try:
            with open(fn) as fi:
                for line in fi:
                    if line.startswith("Sequence length "):
                        return int(line.strip("Sequence length ").split("-")[-1])
        except FileNotFoundError:
            return 0
        return 0
    def checkReadFile(_file):
        return os.stat(_file).st_size if exists(_file) else "NA"

    test = "FASTQC:READCOUNT"
    fastqc_dir = join(QCDIR, "fastqc", readtype, sample)
    fastqc_r1 = join(fastqc_dir, sample + "_R1.{}_fastqc".format(readtype), "fastqc_data.txt")
    fastqc_r2 = join(fastqc_dir, sample + "_R2.{}_fastqc".format(readtype), "fastqc_data.txt")
    read_dir = join(QCDIR, readtype, sample)
    if not (exists(fastqc_r1) and exists(fastqc_r2)):
        R1_size = checkReadFile(join(read_dir, sample + "_R1.{}.fastq.gz".format(readtype)))
        R2_size = checkReadFile(join(read_dir, sample + "_R2.{}.fastq.gz".format(readtype)))
        if R1_size == "NA" or R2_size == "NA":
            print("WARNING: Missing readfiles {}/{}.".format(join(read_dir, sample + "_R1.{}.fastq.gz".format(readtype)), join(read_dir, sample + "_R2.{}.fastq.gz".format(readtype))), file=sys.stderr)
            return(test, "FAIL", "MISSING!", (0, 0, 0))
        if R1_size > 20 or R2_size > 20:
            print("WARNING: One or both FASTQC-report(s) missing for sample {}, but read files seem not to be empty: R1_fastqc={},R1_read_file_size={} R2_fastqc={},R2_read_file_size={}.".format(sample, exists(fastqc_r1), R1_size, exists(fastqc_r2), R2_size), file=sys.stderr)
        return TestResult(test, "FAIL", "MISSING?", (0, 0, 0)) 
    n1, n2 = map(extractReadCount, (fastqc_r1, fastqc_r2))
    readlen = max(map(extractReadLength, (fastqc_r1, fastqc_r2)))
    if n1 != n2: 
        return TestResult(test, "FAIL", "INCONSISTENT", (n1, n2, readlen))
    if n1 < min_reads:
        return TestResult(test, "FAIL", "LOW", (n1, n2, readlen))
    return TestResult(test, "PASS", "", (n1, n2, readlen)) 

def test_kat_peaks(sample, max_peak_volume_threshold=0.9, **kwargs):
    test = "KAT:PEAKS"
    kat_log = join(QCDIR, "kat", sample, sample + ".dist_analysis.json")
    status, errmsg, data = "PASS", "", tuple([None]*8)
    kmer_peaks, gc_peaks, gsize, unit, kmer_freq, mean_gc = None, None, None, None, None, None
    kmer_peak_table = list()
    if not exists(kat_log) or os.stat(kat_log).st_size == 0:
        return TestResult(test, "PASS", "MISSING", data)
    else:
        import json
        try:
            with open(kat_log) as _in:
                kat_data = json.load(_in)
        except:
            return(test, "PASS", "JSON_ERROR", data)

        print(sample, kat_data)

        cov_data = kat_data.get("coverage", dict())
        kmer_peaks = cov_data.get("nb_peaks", None)
        kmer_freq = cov_data.get("mean_freq", None)
        gsize = cov_data.get("est_genome_size", None)
        unit = "bp"
        if gsize is not None:
            unit = "Mbp"
            gsize /= 1e6
        gc_peaks = kat_data.get("gc", dict()).get("nb_peaks", None)
        mean_gc = kat_data.get("gc", dict()).get("mean_gc%", None)

        keys = ("mean_freq", "stddev", "count", "volume")
        kmer_peak_table = list([peak[k] for k in keys] for peak in cov_data.get("peaks", dict(zip(keys, [None]*4))))


        print(sample, kmer_peak_table)

        VOLUME_INDEX = 3

        try:
            max_peak = sorted(kmer_peak_table, key=lambda x:x[VOLUME_INDEX])[-1]
        except: 
            max_peak = [None]*(VOLUME_INDEX + 1)
        try:
            total_volume = sum(row[VOLUME_INDEX] for row in kmer_peak_table)
        except:
            total_volume = 0
            print(kmer_peak_table, file=sys.stderr)

        data = (kmer_peaks, max_peak[VOLUME_INDEX], max_peak[VOLUME_INDEX] / total_volume if (max_peak[VOLUME_INDEX] is not None and total_volume) else None,  gc_peaks, gsize, unit, kmer_freq, mean_gc)
        if total_volume:
            if data[0] is None:
                status, errmsg = "PASS", "MISSING"
            elif data[0] < 1:
                status, errmsg = "PASS", "NOPEAK"
            elif data[0] > 1:
                if max_peak[VOLUME_INDEX] / total_volume < max_peak_volume_threshold:
                    status, errmsg = "PASS", "MULTI_MODAL"
        else:
            status, errmsg = "PASS", "TABLE_EMPTY"

    return TestResult(test, status, errmsg, data)


def test_kat_peaks_old(sample, max_peak_volume_threshold=0.9, **kwargs):
    test = "KAT:PEAKS"
    kat_log = join(QCDIR, "kat", sample, sample + ".kat")
    status, errmsg, data = "PASS", "", tuple([None]*8)
    kmer_peaks, gc_peaks, gsize, unit, kmer_freq, mean_gc = None, None, None, None, None, None
    if not exists(kat_log) or os.stat(kat_log).st_size == 0:
        stats, errmsg = "PASS", "MISSING"
    else:
        with open(kat_log) as _in:
            kmer_peak_table, state, read_table = list(), "", False
            for line in _in:
                # print("LINE=", line, "STATE=", state)
                if line.startswith("ERROR"):
                    stats, errmsg = "PASS", "ERROR"
                    break
                if line.startswith("K-mer frequency spectra statistics"):
                    state = "kmer_freq"
                    try:
                        [next(_in), next(_in)]
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
                            state = ""
                        elif not line.strip(): #  or line.startswith("K-value used: "):
                            read_table = False
                        elif line.startswith("Global") or line.strip().startswith("Index") or line.startswith("--"):
                            continue
                        elif read_table:
                            try:
                                kmer_peak_table.append(list(map(float, re.sub(" +", " ", line.strip()).split(" "))))
                            except:
                                kmer_peak_table.append(list(map(float, re.sub(" +", " ", line.strip()).split(" ")[:-1])))
                    elif state == "gc_dist":
                        if line.startswith("Peaks in analysis: "):
                            try:
                                gc_peaks = int(line.replace("Peaks in analysis: ", "").strip())
                            except:
                                gc_peaks = None
                        elif line.startswith("Mean GC: "):
                                mean_gc = line.replace("Mean GC: ", "").strip()
                    elif state == "genome_stats":
                        # print("PLONK", line)
                        if line.startswith("Estimated genome size: "):
                            line = line.replace("Estimated genome size: ", "").strip()
                            try:
                                gsize, unit = line.split(" ")
                            except:
                                print("XXX", file=sys.stderr)
                                gsize, unit = line, "bp"
                            gsize = float(gsize)

    

            try:
                max_peak = sorted(kmer_peak_table, key=lambda x:x[6])[-1]
            except: 
                max_peak = [None]*7
            try:
                total_volume = sum(row[6] for row in kmer_peak_table)
            except:
                total_volume = 0
                print(kmer_peak_table, file=sys.stderr)
            data = (kmer_peaks, max_peak[6], max_peak[6] / total_volume if (max_peak[6] is not None and total_volume) else None,  gc_peaks, gsize, unit, kmer_freq, mean_gc)
            if total_volume:

                if data[0] is None:
                    status, errmsg = "PASS", "MISSING"
                elif data[0] < 1:
                    status, errmsg = "PASS", "NOPEAK"
                elif data[0] > 1:
                    if max_peak[6] / total_volume < max_peak_volume_threshold:
                        status, errmsg = "PASS", "MULTI_MODAL"
            else:
                status, errmsg = "PASS", "TABLE_EMPTY"

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
 
def test_tadpole_size(sample, min_size=1e6, **kwargs):
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
            
    return TestResult(test, status, errmsg, (str(assembly_size),))
    
TESTS = [("FASTQC:READCOUNT", 
          test_fastqc_readcount,  
          ("Status", "Description", "#R1_reads", "#R2_reads", "Read_Length"), 
          {"min_reads": 1000, "readtype": "bbnorm"}),
         ("KAT:PEAKS", 
           test_kat_peaks, 
           ("Status", "Description", "#K-mer_peaks", "Max_Peak_Volume", "Max_Peak_Volume/total", "GC_peaks", "Genome_size", "Genome_size_unit", "K-mer_freq", "mean_gc"), 
           {"max_peak_volume_threshold": 0.9}),
         ("TADPOLE:SIZE", 
          test_tadpole_size, 
          ("Status", "Description", "Assembly_size"), 
          {"min_size": 1e6})]








def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=str)
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("--min_readcount", type=int, default=1000)
    ap.add_argument("--min_tadpolesize", type=int, default=1e6)
    ap.add_argument("--readtype", type=str, choices=["bbnorm", "bbduk"], default="bbnorm")
    ap.add_argument("--report-dir", type=str, default="")
    ap.add_argument("--override_survey", action="store_true")
    args = ap.parse_args(args_in)
    
    print("Running qc:evaluation...", end="", flush=True)

    if not args.report_dir:
        report_dir = os.path.join(args.indir, "reports")
    else:
        report_dir = os.path.join(args.report_dir)
    pathlib.Path(os.path.join(report_dir, "samplesheets")).mkdir(parents=True, exist_ok=True)
    

    INPUTFILES = dict(readSamplesheet(args.input))
    SHEET = Samplesheet(args.input)
    global QCDIR
    # QCDIR  = join(args.indir, "Analysis", "qc")
    QCDIR = join(args.indir, "qc")
    global QADIR
    QADIR = join(args.indir, "qa", "survey")

    #Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
    
    qc_eval_outf, asm_samplesheet_f = join(report_dir, "qc_eval.tsv"), join(report_dir, "samplesheets", "samplesheet.qc_pass.tsv")
    with open(qc_eval_outf, "w") as qc_eval_out, open(asm_samplesheet_f, "w") as asm_samplesheet:
        # print("SAMPLE", *tuple(test for test, testf, tdata in TESTS), sep="\t", file=qc_eval_out)

        header = ["Sample", "Status"]
        header2 = ["", ""]
        # print("SAMPLE", end="\t", file=qc_eval_out)
        for test, testf, tdata, _ in TESTS:
            header.extend([test] + [""]*(len(tdata)-1))   
            header2.extend(tdata)
            # print(*((test,) + tuple([""]*len(tdata)) for test, testf, tdata in TESTS), sep="\t", file=qc_eval_out)
        print(*header, sep="\t", file=qc_eval_out)
        print(*header2, sep="\t", file=qc_eval_out)

        # for sample in sorted(INPUTFILES):
        keepSamples = set()
        for sample in SHEET:
            results = list()
            for test, testf, tdata, kwargs in TESTS:
                if test == "FASTQC:READCOUNT":
                    kwargs["readtype"] = args.readtype
                results.append(testf(sample, **kwargs))
            # results = list(testf(sample, **kwargs) for test, testf, tdata, kwargs in TESTS)
            passed = all(result[1] == "PASS" for result in results)
            if passed or args.override_survey:

                # /tgac/data/reads/Jay_Hinton_GCRF_Salmonella_LITE/180524_K00287_0041_BHVMHLBBXX/FD01543206_PRO1620_plate92_H12_ACGCAGCAA-AGTCAA_L008_R1.fastq.gz
                # /tgac/workarea/group-pb/schudomc_sandbox/bgrrl_test/Analysis_20180719/qc/bbduk/FD01543280_PRO1620_plate92_F03_GAGTATAAT-AGTCAA_L008/FD01543280_PRO1620_plate92_F03_GAGTATAAT-AGTCAA_L008_R1.bbduk.fastq.gz
                trimdata = [join(QCDIR, "bbduk", SHEET[sample].sampleID, basename(SHEET[sample].R1).replace(".fastq.gz", ".bbduk.fastq.gz")),
                            join(QCDIR, "bbduk", SHEET[sample].sampleID, basename(SHEET[sample].R2).replace(".fastq.gz", ".bbduk.fastq.gz")),
                            ""]
                normdata = [join(QCDIR, "bbnorm", SHEET[sample].sampleID, basename(SHEET[sample].R1).replace(".fastq.gz", ".bbnorm.fastq.gz")),
                            join(QCDIR, "bbnorm", SHEET[sample].sampleID, basename(SHEET[sample].R2).replace(".fastq.gz", ".bbnorm.fastq.gz")),
                            ""]
                if not exists(normdata[0]) or not exists(normdata[1]):
                    normdata = ["", "", ""]
                # asmsample = sample2asmsample(INPUTFILES[sample], trimdata, normdata=normdata)
                # print(asmsample, sep=",", file=asm_samplesheet)
                SHEET[sample].upgrade(ASM_SAMPLE_FIELDS, trimdata + normdata)
                keepSamples.add(sample)
                #print(*INPUTFILES[sample], sep=",", file=asm_samplesheet)


            elif results[0][1] == "FAIL" and results[0][2] == "MISSING" and results[2][1] == "PASS":
                print("WARNING: sample {} has missing fastqc report(s), but passes survey assembly.".format(sample), file=sys.stderr)

            # p_results = tuple(",".join(result) for result in results)
            # print(p_results)
            row = [sample, "PASS" if passed else "FAIL"]
            for result in results:
                row.extend(result[1:3] + result[3])
            print(*row, sep="\t", file=qc_eval_out)
            # p_result = tuple("\t".join(map(str, result[1:])) for result in results)          

            # print(sample, *p_result, sep="\t", file=qc_eval_out)
        SHEET.write(asm_samplesheet, keepSamples)


    print(" Done.\n Generated qc_eval report in {} and asm-samplesheet in {}.".format(qc_eval_outf, asm_samplesheet_f))

    return True

if __name__ == "__main__":
    main()
