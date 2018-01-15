#!/usr/bin/env python
import sys
import os
import glob
import re

from collections import Counter, namedtuple

def checkReadCount(r1, r2, threshold=1000):
    def getReadCount(fn):
        with open(fn) as fi:
            for line in fi:
                if line.startswith("Total Sequences"):
                    return int(line.strip().split("\t")[-1])
    n1, n2 = map(getReadCount, (r1, r2))
    return n1 == n2 and n1 >= threshold, n1, n2

READCOUNT_THRESHOLD = 1000
TADPOLE_SIZE_THRESHOLD = 1000000

ROOT = sys.argv[1]
INPUTDIR = os.path.join(ROOT, "Reads")
QCDIR = os.path.join(ROOT, "Analysis", "qc")
QCLOGDIR = os.path.join(QCDIR, "log")

INPUTFILES = dict((os.path.basename(_file).replace("_R1.fastq.gz", ""),
                   (_file, _file.replace("_R1.fastq.gz", "_R2.fastq.gz")))
                  for _file in glob.glob(os.path.join(INPUTDIR, "*_R1.fastq.gz")))

sample_counter = Counter()

TestResult = namedtuple("TestResult", "test status errmsg data".split(" "))

def test_fastqc_readcount(sample):
    def extractReadCount(fn):
        with open(fn) as fi:
            for line in fi:
                if line.startswith("Total Sequences"):
                    return int(line.strip().split("\t")[-1])
        return 0

    test = "FASTQC:READCOUNT"
    fastqc_dir = os.path.join(QCDIR, "fastqc", sample)
    fastqc_r1 = os.path.join(fastqc_dir, sample + "_R1.bbduk_fastqc", "fastqc_data.txt")
    fastqc_r2 = os.path.join(fastqc_dir, sample + "_R2.bbduk_fastqc", "fastqc_data.txt")   
    if not os.path.exists(fastqc_r1) and os.path.exists(fastqc_r2):
        return TestResult(test, "FAIL", "MISSING", (0, 0)) 
    n1, n2 = map(extractReadCount, (fastqc_r1, fastqc_r2))
    if n1 != n2: 
        return TestResult(test, "FAIL", "INCONSISTENT", (n1, n2))
    if n1 < READCOUNT_THRESHOLD:
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
    
def test_tadpole_size(sample):
    test = "TADPOLE:SIZE"
    tadpole_contigs = os.path.join(QCDIR, "tadpole", sample, "tadpole_contigs.fasta")   
    if os.path.exists(tadpole_contigs) and os.stat(tadpole_contigs).st_size < TADPOLE_SIZE_THRESHOLD:
        return TestResult(test, "FAIL", "TOO_SMALL", (os.stat(tadpole_contigs).st_size, TADPOLE_SIZE_THRESHOLD))
    return TestResult(test, "PASS", "", (os.stat(tadpole_contigs).st_size, TADPOLE_SIZE_THRESHOLD))

TESTS = [("FASTQC:READCOUNT", test_fastqc_readcount),
         ("KAT:HIST", test_kat_hist),
         ("KAT:GCP", test_kat_gcp),
         ("TADPOLE:SIZE", test_tadpole_size)]

print("\t".join(("SAMPLE",) + tuple(test for test, testf in TESTS)))
for sample in sorted(INPUTFILES):
    results = tuple(",".join(map(str, testf(sample)[1:])) for test, testf in TESTS)
    print("\t".join((sample, ) + results))


'''
print(sample_counter["PASS"], "samples of", sum(sample_counter.values()), "({:.3f}%)".format(sample_counter["PASS"]/sum(sample_counter.values())*100), "samples ready for assembly.")
    
'''

    
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
