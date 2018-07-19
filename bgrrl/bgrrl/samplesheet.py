import os
from collections import namedtuple, OrderedDict
import csv

RAW_SAMPLE_FIELDS = [
    "sampleID",
    "customerSampleID",
    "R1",
    "R2",
    "S",
    "taxonomyID",
    "taxonomyTxt",
    "fastqcR1",
    "fastqcR2",
    "fastqcS"
]

ASM_SAMPLE_FIELDS = [
    "R1trim",
    "R2trim",
    "Strim",
    "R1norm",
    "R2norm",
    "Snorm"
]

Sample = namedtuple("Sample", RAW_SAMPLE_FIELDS)
ASM_Sample = namedtuple("ASM_Sample", RAW_SAMPLE_FIELDS + ASM_SAMPLE_FIELDS)

# Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
# ASM_Sample = namedtuple("ASM_Sample", "sampleID customerSampleID R1 R2 S R1trim R2trim Strim R1norm R2norm Snorm taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))




def sample2asmsample(sample, trimdata, normdata=["", "", ""]):
    return ASM_Sample(*sample, *trimdata, *normdata)


class BaseSample(object):
    def __init__(self, *args, sampleFields=RAW_SAMPLE_FIELDS):
        self._sample_field_order = list(sampleFields)
        if not len(args) == len(sampleFields):
             raise ValueError("SAMPLE ERROR: Number of arguments ({}) does not match number of fields ({}).".format(len(args), len(sampleFields)))
        for argid, arg in zip(sampleFields, args):
            setattr(self, argid, arg)
    def verifyDatasets(self, fields=["R1", "R2"]):
        for f in fields:
            path = getattr(self, f)
            if path and not os.path.exists(path):
                raise ValueError("SAMPLE ERROR: Sample {}. Cannot find {} data at {}.".format(s, f, path))       
    def upgrade(self, sampleFields, sampleData):
        if not len(sampleFields) == len(sampleData):
            raise ValueError("SAMPLE UPGRADE ERROR: Number of fields ({}) does not match expected number of fields ({}).".format(len(sampleData), len(sampleFields)))
        for argid, arg in zip(sampleFields, sampleData):
            setattr(self, argid, arg)
        self._sample_field_order.extend(sampleFields)
    def __iter__(self):
        for fname in self._sample_field_order:
            yield getattr(self, fname)

class ASM_Sample(BaseSample):
    def __init__(self, *args, sampleFields=(RAW_SAMPLE_FIELDS,ASM_SAMPLE_FIELDS)):
        super(ASM_Sample, self).__init__(*args[:len(sampleFields[0])], sampleFields=RAW_SAMPLE_FIELDS)
        self._sample_field_order = sampleFields[0] + sampleFields[1]
        if not len(args[len(sampleFields[0]):]) == len(sampleFields[1]):
            print(args[len(sampleFields[0]):], sampleFields[1], file=sys.stderr)
            raise ValueError("ASM_SAMPLE ERROR: Number of arguments ({}) does not match number of fields ({}).".format(len(args[len(sampleFields[0]):]), len(sampleFields[1])))
        for argid, arg in zip(sampleFields[1], args[len(sampleFields[0]):]):
            setattr(self, argid, arg)


class Samplesheet(OrderedDict):
    def __init__(self, _input, sampletype=BaseSample):
        super(Samplesheet, self).__init__()
        if type(_input) is str:
            if not os.path.exists(_input):
                raise ValueError("SAMPLESHEET ERORR: input file {} does not exist.".format(_input))
            with open(_input) as _in:
                for r in csv.reader(_in, delimiter=","):
                    print(r)
                    self[r[0]] = sampletype(*r)
    def write(self, stream, _filter=set()):
        for s in self:
            if not _filter or s in _filter:
                print(*tuple(self[s]), sep=",", file=stream)
                print("")
    def verifySampleData(self, fields=["R1", "R2", "R1trim", "R2trim", "R1norm", "R2norm"]):
        for s in self:
            self[s].verifyDatasets(fields=fields)
        return True



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
