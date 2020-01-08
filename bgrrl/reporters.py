import os
import sys
import csv
import glob
from collections import Counter


def collate_quast_reports(reportfile, *quast_reports):
	header = ""
	with open(reportfile, "w") as report_out:
		for f in quast_reports:
			for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
				if i == 0 and not header:
					header = row
					print(*header, file=report_out, sep="\t")
				elif i:
					print(*row, file=report_out, sep="\t")

def collate_assembly_information(assembly_dir, stage_report_file, stage_summary_file, samplesheet_file, verbose=True, stage_prefix="asm_", supported_stages=dict()):
	assembly_stages = list()
	with open(stage_report_file, "w") as stage_out, open(samplesheet_file, "w") as samplesheet_out:
		for cdir, _, _ in os.walk(assembly_dir):
			sample = os.path.basename(cdir)
			if sample != "log" and os.path.dirname(cdir) == assembly_dir:
				try:
					stage = glob.glob(os.path.join(cdir, stage_prefix))[0]
				except:
					stage = "NA"
				assembly_stages.append(stage)
				if stage != "NA":
					print(sample, sample, os.path.join(cdir, sample + ".assembly.fasta"), sep=",", file=samplesheet_out)
				print(sample, stage, sep="\t", file=stage_out)

	assembly_stages = Counter(assembly_stages)
	if verbose:
		print(assembly_stages)
	with open(stage_summary_file, "w") as summary_out:
		N = sum(assembly_stages.values())
		for stage_id, stage in supported_stages.items():
			if assembly_stages[stage_id]:
				print(stage, assembly_stages[stage_id], assembly_stages[stage_id]/N, sep="\t", file=summary_out)
		print("Total", "", "", N, sep="\t", file=summary_out)
