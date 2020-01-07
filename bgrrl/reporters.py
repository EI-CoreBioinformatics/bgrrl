import os
import sys
import csv

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
