import sys
import subprocess
import argparse
import pathlib
import os
from os.path import join

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("--singularity-container", type=str)
	ap.add_argument("--input", "-i", type=str)
	ap.add_argument("--outdir", "-o", type=str)
	ap.add_argument("--busco-mode", choices=("geno", "genome", "g", "tran", "transcriptome", "t", "prot", "proteome", "p"), default="genome")
	ap.add_argument("--threads", type=int, default=1)
	ap.add_argument("--lineage", "-l", type=str)

	args = ap.parse_args()

	
	outdir = join(os.path.dirname(args.outdir), "run_" + os.path.basename(args.outdir))
	tmpdir = join(outdir, "tmp")
	augustus_dir = join(tmpdir, "config")
	pathlib.Path(tmpdir).mkdir(exist_ok=True, parents=True)
	
	singularity_call = "singularity exec " + args.singularity_container  

	try:
		subprocess.check_call(
			"{} cp -rv /opt/miniconda/config {}/".format(singularity_call, os.path.abspath(tmpdir)),
			shell=True
		)
	except Exception as e:
		print("EXCEPTION1:", e)
		sys.exit(1)

	cmd = "cd {} && " + \
		"export AUGUSTUS_CONFIG_PATH={} && " + \
		"{} run_BUSCO.py -i {} -c {} -m {} --force -t {} -l {} -o {}"
	cmd = cmd.format(
		os.path.abspath(os.path.join(outdir, "..")),
		os.path.abspath(augustus_dir),
		singularity_call,
		os.path.abspath(args.input),
		args.threads,
		{"g": "geno", "t": "tran", "p": "prot"}.get(args.busco_mode[0], "geno"),
		os.path.abspath(tmpdir),
		args.lineage,
		os.path.basename(args.outdir)
	)

	try:
		subprocess.check_call(cmd, shell=True)
	except Exception as e:
		print(e)
		sys.exit(1)
		

	os.rename(outdir, args.outdir)
	open(os.path.join(args.outdir, "short_summary_" + os.path.basename(args.outdir) + ".txt"), "at").close()
	
	os.rename(
		os.path.join(args.outdir, "short_summary_" + os.path.basename(args.outdir) + ".txt"),
		os.path.join(args.outdir, os.path.basename(args.outdir) + "_short_summary.txt")
	)
	
