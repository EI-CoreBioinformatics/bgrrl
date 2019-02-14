import csv
import argparse
import os

def main():
	
	ap = argparse.ArgumentParser()
	ap.add_argument("qc_sheet", type=str)
	ap.add_argument("asm_path", type=str)
	
	args = ap.parse_args()

	if not os.path.exists(args.qc_sheet):
		raise ValueError("Samplesheet {} does not seem to exist. Please check.".format(args.qc_sheet))

	if not os.path.exists(args.asm_path):
		raise ValueError("Path to assembly data does not seem to exist. Please check.".format(args.asm_path))

	with open(args.qc_sheet.replace(".qc_pass.tsv", ".survey_pass.tsv").replace(".survey_pass.tsv", ".asm_pass.tsv"), "wt") as asm_pass_out:
		for row in csv.reader(open(args.qc_sheet), delimiter=","):
			out_row = row[:2] + [os.path.join(args.asm_path, row[0], row[0] + ".assembly.fasta")]
			print(*out_row, sep=",", file=asm_pass_out)

	


	pass



if __name__ == "__main__":
	main()
