#!/usr/bin/env python
import sys
import os
from os.path import exists, join, basename, dirname
import re
import csv
import argparse
import pathlib 
import yaml

from collections import Counter, namedtuple

from bgrrl.samplesheet import *

TestResult = namedtuple("TestResult", "test status errmsg data".split(" "))

def test_fastqc_readcount(sample, workdir, min_reads=1000, readtype="bbnorm", **kwargs):
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
		return os.stat(_file).st_size if exists(_file) else None

	test = "FASTQC:READCOUNT"
	fastqc_dir = join(workdir, "qc", "fastqc", readtype, sample)
	read_dir = join(workdir, "qc", readtype, sample)

	fastqc_r1 = join(fastqc_dir, sample + "_R1.{}_fastqc".format(readtype), "fastqc_data.txt")

	def check_reads(read_dir, sample, readtype, mate, test="FASTQC:READCOUNT"):
		read_file = join(read_dir, "{}_R{}.{}.fastq.gz".format(sample, mate, readtype))
		read_size = checkReadFile(read_file)
		if read_size is None:
			print("WARNING: Missing R{} read file at {}".format(mate, read_file))
			return TestResult(test, "FAIL", "R{}_MISSING".format(mate), (0, 0, 0))
		if read_size > 20:
			print("WARNING: One or both FASTQC-report(s) missing for sample {0}, but read file {3} does seem not to be empty: R{1}_read_file_size={2}".format(sample, mate, read_size, read_file), file=sys.stderr)
		return TestResult(test, "FAIL", "R{}_MISSING".format(mate), (0, 0, 0))

	if not exists(fastqc_r1):
		return check_reads(read_dir, sample, readtype, 1)

	n1, n2 = extractReadCount(fastqc_r1), "NA"
	readlen = extractReadLength(fastqc_r1)

	if not kwargs.get("single_cell", False):
		fastqc_r2 = join(fastqc_dir, sample + "_R2.{}_fastqc".format(readtype), "fastqc_data.txt")
		if not exists(fastqc_r2):
			return check_reads(read_dir, sample, readtype, 2)

		n2 = extractReadCount(fastqc_r2)
		readlen = max(readlen, extractReadLength(fastqc_r2))

		if n1 != n2: 
			return TestResult(test, "FAIL", "INCONSISTENT", (n1, n2, readlen))

	if n1 < min_reads:
		return TestResult(test, "FAIL", "LOW", (n1, n2, readlen))

	return TestResult(test, "PASS", "", (n1, n2, readlen)) 

def test_kat_peaks(sample, workdir, max_peak_volume_threshold=0.9, **kwargs):
	test = "KAT:PEAKS"
	kat_log = join(workdir, "qc", "kat", sample, sample + ".dist_analysis.json")
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

		#print(sample, kat_data)

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


		#print(sample, kmer_peak_table)

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
 
def test_tadpole_size(sample, workdir, min_size=1e6, **kwargs):
	def extractAssemblySize(qreport_in):
		for row in csv.reader(qreport_in, delimiter="\t"):
			if row[0].startswith("Total length"):
				return int(row[1])
		return 0
	test = "TADPOLE:SIZE"
	# quast_report = join(workdir, "qaa", "survey", "quast", sample, "report.tsv")
	quast_report = join(workdir, "qc", "tadpole", sample, "quast", "report.tsv")

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

class SurveyStageEvaluation:
	def __init__(self, workdir):
		print("Running survey stage evaluation...")
		self.report_dir = join(workdir, "reports")
		self.samplesheet_dir = join(self.report_dir, "samplesheets")
		pathlib.Path(self.samplesheet_dir).mkdir(exist_ok=True, parents=True)

		self.workdir = workdir
		self.survey_evaluation_results_file = join(self.report_dir, "survey_stage_evaluation.tsv")
		self.samplesheet_file = join(self.samplesheet_dir, "samplesheet.survey_pass.yaml")
	
	def __validate_sample_data(self):	
		survey_data_path = join(self.workdir, "qc", "bbduk")
		walker = os.walk(survey_data_path)
		try:
			_, samples, _ = next(walker)
		except StopIteration:
			raise ValueError("Could not find survey data at {}.".format(survey_data_path))

		self.sample_data = dict()
		for sample_path, dirs, files in walker:
			sample_id = os.path.basename(sample_path)
			self.sample_data[sample_id] = {
				"bbduk_r1": join(sample_path, sample_id + "_R1.bbduk.fastq.gz"),
				"bbduk_r2": join(sample_path, sample_id + "_R2.bbduk.fastq.gz"),
				"bbduk_singles": join(sample_path, sample_id + "_S.bbduk.fastq.gz"),
				"tadpole_singles": join(sample_path, sample_id + "_S.bbduk.corr.fastq.gz"),
				"bbmerge": join(sample_path, sample_id + ".merged.fastq.gz"),
				"bbmerge_ur1": join(sample_path, sample_id + "_R1.bbduk.corr.unmerged.fastq.gz"),
				"bbmerge_ur2": join(sample_path, sample_id + "_R2.bbduk.corr.unmerged.fastq.gz"),
				"uc_singles": join(sample_path, sample_id + ".uc_singles.fastq.gz")
			}
			for k, fn in self.sample_data[sample_id].items():
				if not os.path.exists(fn):
					print("Missing file {}".format(fn), file=sys.stderr)
					self.sample_data[sample_id][k] = None

	def __evaluate_survey_results(self, all_tests=TESTS, min_tadpole_size=1e6, override_survey=False):
		with open(self.survey_evaluation_results_file, "w") as survey_out:
			header, header2 = ["Sample", "Status"], ["", ""]
			for test, testf, tdata, _ in TESTS:
				header.extend([test] + [""]*(len(tdata)-1))
				header2.extend(tdata)
			print(*header, sep="\t", file=survey_out)
			print(*header2, sep="\t", file=survey_out)

			discard_samples = set()
			for sample in self.sample_data:
				results = list()
				for test, testf, tdata, kwargs in TESTS:
					if test == "FASTQC:READCOUNT":
						kwargs["readtype"] = "bbduk"
					if test == "TADPOLE:SIZE":
						kwargs["min_size"] = min_tadpole_size
					results.append(testf(sample, self.workdir, **kwargs))

				passed = all(result[1] == "PASS" for result in results)
				if passed or override_survey:
					pass
				elif results[0][1] == "FAIL" and results[0][2] == "MISSING" and results[2][1] == "PASS":
					print("WARNING: sample {} has missing fastqc report(s), but passes survey assembly.".format(sample), file=sys.stderr)
				else:
					discard_samples.add(sample)

				row = [sample, "PASS" if passed else "FAIL"]
				for result in results:
					row.extend(result[1:3] + result[3])
				print(*row, sep="\t", file=survey_out)

				self.sample_data = dict((k, v) for k, v in self.sample_data.items() if k not in discard_samples)


	def __generate_samplesheet(self):
		with open(self.samplesheet_file, "w") as samplesheet:
			yaml.dump(self.sample_data, samplesheet, default_flow_style=False)
		

	def run(self, min_tadpole_size=1e6, override_survey=False):
		self.__validate_sample_data()
		self.__evaluate_survey_results(min_tadpole_size=min_tadpole_size, override_survey=override_survey)
		self.__generate_samplesheet()
		print(" Done.\n Generated survey stage evaluation report in {} and asm-samplesheet in {}.".format(
			self.survey_evaluation_results_file,
			self.samplesheet_file
		))



def main(args_in=sys.argv[1:]):
	ap = argparse.ArgumentParser()
	ap.add_argument("indir", type=str, default=".")
	ap.add_argument("--min_readcount", type=int, default=1000)
	ap.add_argument("--min_tadpole_size", type=int, default=1e6)
	ap.add_argument("--readtype", type=str, choices=["bbnorm", "bbduk"], default="bbnorm")
	ap.add_argument("--report-dir", type=str, default="")
	ap.add_argument("--override_survey", action="store_true")
	ap.add_argument("--single-cell", action="store_true")
	args = ap.parse_args(args_in)

	SSE = SurveyStageEvaluation(args.indir)
	SSE.run(min_tadpole_size=args.min_tadpole_size, override_survey=args.override_survey)
	return True


if __name__ == "__main__":
	main()
