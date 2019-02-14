import sys
import os
from os.path import join
import argparse
import subprocess
import shutil
import pathlib

ACTIVATE_PROKKA_ENV = "source activate bgasm_env"


"""
old:
PROKKAPATH=/tgac/software/testing/prokka/git/x86_64
export PATH=$PATH:$PROKKAPATH/binaries/linux

# source hmmer-3.1b2
# source gnuparallel-20131122
source perl-5.20.1_gk
source jre-8u92

/usr/bin/time -v $PROKKAPATH/bin/prokka --cpus $5 --outdir $1 --prefix $2 --centre EI $3 --force &> $4
"""

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("contigs", type=str)
	ap.add_argument("--outdir", type=str, default=".")
	ap.add_argument("--prefix", type=str, default="sample")
	ap.add_argument("--seq-centre", type=str, default="EI")
	ap.add_argument("--force", action="store_true")
	ap.add_argument("--threads", type=int, default=8)
	ap.add_argument("--singularity-container", type=str, default="")
	ap.add_argument("--proteins", type=str, default="")

	args = ap.parse_args()

	singularity_prefix = ("singularity exec " + args.singularity_container + " ") if args.singularity_container else ""


	cmd = "mkdir -p {1}" + \
		" && cp -v {4} {1}/prokka_contigs.fasta" + \
		" && " + singularity_prefix + \
		"prokka" + \
		" --cpus {0}" + \
		" --outdir {1}" + \
		" --prefix {2}" + \
		" --centre {3}" + \
		(" --force" if args.force else "") + \
		((" --proteins " + args.proteins) if args.proteins else "") + \
		" {1}/prokka_contigs.fasta" + \
		" && rm {1}/prokka_contigs.fasta"

	try:
		host = subprocess.check_output("hostname", shell=True, stderr=subprocess.STDOUT).decode()
	except:
		host = "n/a"

	prokka_msg = "UNINITIATED PROKKA MESSAGE"
	
	print(
		cmd.format(
			args.threads,
			args.outdir,
			args.prefix,                
			args.seq_centre,
			args.contigs
		)
	)






	try:
		prokka_msg = subprocess.check_output(
			cmd.format(
				args.threads,
				args.outdir,
				args.prefix,
				args.seq_centre,
				args.contigs
			),
			stderr=subprocess.STDOUT,
			shell=True
		).decode()
	except subprocess.CalledProcessError as e:
		prokka_msg = e.output.decode()
		files = ["log", "faa", "ffn", "gff"]
		for sfx in files:
			open(join(args.outdir, args.prefix + "." + sfx), "wt").close()

		with open(join(args.outdir, "PROKKA_FAILED"), "wt") as prokka_fail_out:
			prokka_fail_out.write(host)


	print(prokka_msg)


		
				
