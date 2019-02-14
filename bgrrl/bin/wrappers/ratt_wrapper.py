import sys
import os
from os.path import basename, exists, dirname, join
import shutil
import argparse
import subprocess
import pathlib
import glob
import re
import csv

# done = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann, wildcards.sample + "_" + wildcards.ref_ann + ".done")

#export RATT_HOME=/tgac/software/testing/ratt/dev_no_defined/x86_64/bin/
			#export RATT_CONFIG=$RATT_HOME/RATT.config_bac


class RATTRunner(object):
	def __init__(self, contigs, reference, path_to_ratt, outdir, container=""):
		self.contigs = contigs
		self.reference = reference
		self.path_to_ratt = path_to_ratt
		self.outdir = outdir
		self.sample = basename(dirname(self.contigs))
		self.prefix = "{}_{}".format(self.sample, basename(self.reference))
		self.logfile = open(join(self.outdir, "ratt_runner.log"), "wt")
		self.container = container

		pathlib.Path(self.outdir).mkdir(parents=True, exist_ok=True)

	def run(self):
		print("Running RATT ... ", end="", flush=True, file=self.logfile)

		if self.container:
			cmd = "cd {0}; /usr/bin/time -v " + self.container + " start.ratt.sh {1} {2} {3} Strain"
			cmd = cmd.format(self.outdir, self.reference, self.contigs, self.prefix)
		else:
			cmd = "export RATT_HOME={0}; " + \
				"export RATT_CONFIG={0}/RATT.config_bac; " + \
				"cd {1}; " + \
				"/usr/bin/time -v {0}/start.ratt.sh {2} {3} {4} Strain"
			cmd = cmd.format(self.path_to_ratt, self.outdir, self.reference, self.contigs, self.prefix)

		try:
			ratt_msg = subprocess.check_output(
				cmd,
				stderr=subprocess.STDOUT,
				shell=True
			).decode()
			self.do_postprocess = True
			print("done", file=self.logfile)
		except subprocess.CalledProcessError as e:
			ratt_msg = e.output.decode()
			self.do_postprocess = False
			print("FAILED", file=self.logfile)

		print("RATT_MESSAGE")
		print(ratt_msg)
		print("RATT_MESSAGE_END")
	
		return self.do_postprocess

	def postprocess(self):
		def __setup_directories(outdir):
			dirs = ["nucmer", "artemis/not_transfered"]
			for d in dirs:
				d = join(outdir, d)
				pathlib.Path(d).mkdir(exist_ok=True, parents=True)

		def __clean_outdir(outdir):
			file_patterns = ["DONUCMER.log", "query.*", "Reference.*.fasta", "*.tmp2.embl"]
			for fp in file_patterns:
				for f in glob.glob(join(outdir, fp)):
					try:
						os.remove(f)
					except:
						pass

			for f in glob.glob(join(outdir, "*.embl")):
				if not f.endswith(".final.embl") and not f.endswith(".not_transfered.embl"):
					os.remove(f)
			for f in glob.glob(join(outdir, "*.done*")):
				os.remove(f)

			dir_patterns = ["Sequences", "Query", "Reference"]
			for d in dir_patterns:
				try:
					shutil.rmtree(join(outdir, d))
				except:
					pass

		def __handle_transfers(merged_file, destdir, fpattern):
			tmpfile = merged_file + ".tmp"
			with open(tmpfile, "wt") as tmp_out:
				for f in glob.glob(fpattern):
					print(open(f, "rt").read().strip(), file=tmp_out)
					shutil.move(f, join(destdir, basename(f)))
			shutil.move(tmpfile, merged_file)

		def __handle_txt_reports(outdir, reference, sample, prefix):
			with open(join(outdir, prefix + ".report.tsv"), "wt") as out:
				got_header = False
				for f in glob.glob(join(outdir, "*.Report.txt")):
					if os.stat(f).st_size > 166:
						with open(f) as _in:
							for j, r in enumerate(csv.reader(_in, delimiter="\t")):
								col1, col2 = "", ""
								if j == 0:
									if not got_header:
										col1, col2 = "Sample", "Reference"
									got_header = True
								else:
									col1, col2 = sample, reference

								if col1 and col2:
									print(col1, col2, *r, sep="\t", file=out)
					os.remove(f)

		def __handle_gff_reports(outdir, prefix):
			with open(join(outdir, prefix + ".report.gff"), "wt") as out:
				for f in glob.glob(join(outdir, "*.Report.gff")):
					if os.stat(f).st_size > 0:
						print("#", basename(f).replace(".Report.gff", ""), file=out)
						print(open(f, "rt").read().strip(), file=out)
					os.remove(f)

		def __run_seqret_conversion(outdir, prefix, container):
			cmd = "cd {0}; " + \
				"/usr/bin/time -v " + container + " seqret -sequence {1}.final.embl " + \
				"-outseq {1}.final.gff -offormat gff -feature"

			print("Running seqret:")
			print(cmd.format(outdir, prefix))

			try:
				seqret_msg = subprocess.check_output(
					cmd.format(outdir, prefix),
					stderr=subprocess.STDOUT,
					shell=True
				).decode()
				do_final = True
			except subprocess.CalledProcessError as e:
				seqret_msg = e.output.decode()
				do_final = False
															
			print("SEQRET MESSAGE")
			print(seqret_msg)

		def __correct_offsets(outdir, prefix):

			with open(join(outdir, prefix + ".final.gff.raw"), "wt") as out_raw, \
					open(join(outdir, prefix + ".final.gff"), "wt") as out, \
					open(join(outdir, prefix + ".final.gff.delta"), "wt") as out_delta:

				print("##gff-version 3", file=out_raw)
				print("##gff-version 3", file=out)

				reader = csv.reader(open(join(outdir, "final.gff")), delimiter="\t")
				for row in reader:
					if row[0].startswith("#"):
						if row[0].startswith("##gff-version"):
							seq_region = next(reader)
							date = next(reader)
							ftype = next(reader)
						else:
							continue		
					else:
						if row[2] != "contig":
							row[0] = re.sub(".final$", "", row[0])
							row[0] = re.sub("^" + prefix + ".", "", row[0])
						
						if seq_region:
							print(*seq_region, sep="\t", file=out_raw)
							print(*seq_region, sep="\t", file=out)
							print(*date, sep="\t", file=out_raw)
							print(*date, sep="\t", file=out)
							print(*ftype, sep="\t", file=out_raw)
							print(*ftype, sep="\t", file=out)
							seq_region = data = ftype = ""
					
			
						print(*row, sep="\t", file=out_raw)
						if row[3] == "0":
							print(*row, sep="\t", file=out_delta)
							row[3] = "1"
							print(*row, sep="\t", file=out_delta)
							print("---", file=out_delta)	
						print(*row, sep="\t", file=out)


		assert self.do_postprocess

		print("Setting up directories ... ", end="", flush=True, file=self.logfile)
		__setup_directories(self.outdir)
		print("done.", file=self.logfile)

		print("Moving nucmer output to nucmer/ subdirectory ... ", end="", flush=True, file=self.logfile)
		for f in glob.glob(join(self.outdir, "nucmer.*")):				
			shutil.move(f, join(self.outdir, "nucmer", basename(f)))
		print("done.", file=self.logfile)

		print("Compiling transfer reports ... ", end="", flush=True, file=self.logfile)
		__handle_transfers(
			join(self.outdir, self.prefix + ".final.embl"), 
			join(self.outdir, "artemis"), 
			join(self.outdir, "*.final.embl")
		)		

		__handle_transfers(
			join(self.outdir, self.prefix + ".not_transfered.embl"),
			join(self.outdir, "artemis/not_transfered"),
			join(self.outdir, "*.NOTTransfered.embl")
		)
		
		__handle_txt_reports(self.outdir, self.reference, self.sample, self.prefix)

		__handle_gff_reports(self.outdir, self.prefix)
		print("done.", file=self.logfile)

		print("Cleaning output directory ... ", end="", flush=True, file=self.logfile)
		__clean_outdir(self.outdir)
		print("done.", file=self.logfile)

		print("Running seqret to convert transfers to gff ... ", end="", flush=True, file=self.logfile)
		__run_seqret_conversion(self.outdir, self.prefix, self.container)
		print("done.", file=self.logfile)

		print("Correcting zero-offsets ...", end="", flush=True, file=self.logfile)
		__correct_offsets(self.outdir, self.prefix)
		print("done.", file=self.logfile)

	
	

def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("--outdir", "-o", type=str, default="ratt_out")
	ap.add_argument("contigs", type=str)
	ap.add_argument("reference", type=str)
	ap.add_argument("path_to_ratt", type=str)
	ap.add_argument("--singularity-container", type=str, default="")

	args = ap.parse_args()

	singularity_prefix = ("singularity exec " + args.singularity_container + " ") if args.singularity_container else ""


	ratt_runner = RATTRunner(	
		args.contigs, args.reference, args.path_to_ratt, args.outdir, container=singularity_prefix
	)
	
	if not ratt_runner.run():
		raise ValueError("RATT Run failed ({})".format(str(args)))

	ratt_runner.postprocess()


if __name__ == "__main__":
	main()
