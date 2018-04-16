from collections import namedtuple
import csv

Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))

def readSamplesheet(fn, delimiter=","):
    with open(fn) as fi:
        for row in csv.reader(fi, delimiter=delimiter):
            sample = Sample(*row)
            assert sample.sampleID
            if not sample.customerSampleID:
                row[1] = sample.sampleID
            yield (sample.sampleID, Sample(*row))

def verifySamplesheet(fn, delimiter=","):
    for sample_id, sample in readSamplesheet(fn, delimiter=delimiter):
        r1_exists, r2_exists, s_exists = map(os.path.exists, (sample.R1, sample.R2, sample.S))
        if not (r1_exists and r2_exists):
            raise ValueError("Cannot find R1/R2 data at R1={}, R2={}.".format(sample.R1, sample.R2))
        return True
