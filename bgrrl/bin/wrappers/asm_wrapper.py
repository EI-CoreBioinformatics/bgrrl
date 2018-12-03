import sys
import os
import argparse
import subprocess
import shutil
import pathlib

ACTIVATE_ASM_ENV = "source activate bgasm_env" 

class ASM_Wrapper(object):
	def __init__(self, args):
		# self.assembler = args.assembler
		r1_reads, r2_reads = args.r1.split(","), args.r2.split(",")
	
		self.reads_r1 = r1_reads.pop(0)
		self.fb_reads_r1 = r1_reads.pop(0) if r1_reads else ""

		self.reads_r2 = r2_reads.pop(0)
		self.fb_reads_r2 = r2_reads.pop(0) if r2_reads else ""
		
		self.outdir = args.outdir
		self.threads = args.threads

	def run(self):
		self.fallback_queue = [item for item in self.fallback_queue if not item is None]

		assembly_complete = False

		while self.fallback_queue:
			step, cmd = self.fallback_queue.pop(0)
			print("BGRR|_ASM_WRAPPER::" + step)

			try:
				out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
			except CalledProcessError:
				print(out)
				shutil.rmtree(self.outdir, ignore_errors=True)
				continue

			print(out)

			open(os.path.join(self.outdir, step), "wt").close()
			assembly_complete = True
			break

		if not assembly_complete:
			pathlib.Path(self.outdir).mkdir(parents=True, exists_ok=True)
			open(os.path.join(self.outdir, "assembly.fasta"), "wt").close()		









class VelvetWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)
		
		tmpdir = os.path.join(self.outdir, "tmp")
		pathlib.Path(tmpdir).mkdir(parents=True, exists_ok=True)

		cwd = os.getcwd()

		cmd = "cd {0}" + \
			" && VelvetOptimiser.pl -v -s 19 -e 191 -x 2 -t {1} -d ../velvet_assembly -f -shortPaired -fastq.gz -separate {2} {3}" + \
			" && cd .." + \
			" && mv -v velvet_assembly/* ." + \
			" && rm -rf velvet_assembly/ tmp/" + \
			" && ln -s contigs.fa assembly.fasta" + \
			" && cd {4}"

		self.fallback_queue = [
			("asm_main_ven", cmd.format(tmpdir, self.threads, os.path.join(cwd, self.reads_r1), os.path.join(cwd, self.reads_r2), cwd)),
			("asm_fb1_vet", cmd.format(tmpdir, self.threads, os.path.join(cwd, self.fb_reads_r1), os.path.join(cwd, self.fb_reads_r2), cwd)) if (self.fb_reads_r1 and self.fb_reads_r2) else None
		]
		
class UnicyclerWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)

		spades_params = "--careful --cov-cutoff auto"
		unicycler_params = "--min_polish_size 1000"

		self.fallback_queue = [
			("asm_main_ucn", "unicycler -1 {0} -2 {1} -t {2} -o {3} {4}".format(self.reads_r1, self.reads_r2, self.threads, self.outdir, unicycler_params)),
			("asm_fb1_uct", "unicycler -1 {0} -2 {1} -t {2} -o {3} {4}".format(self.fb_reads_r1, self.fb_reads_r2, self.threads, self.outdir, unicycler_params)) if (self.fb_reads_r1 and self.fb_reads_r2) else None,
			("asm_fb2_spn", "spades.py -1 {0} -2 {1} -t {2} -o {3} -m 64 {4} && ln -s scaffolds.fasta {3}/assembly.fasta".format(self.reads_r1, self.reads_r2, self.threads, self.outdir, spades_params)),
			("asm_fb3_spt", "spades.py -1 {0} -2 {1} -t {2} -o {3} -m 64 {4} && ln -s scaffolds.fasta {3}/assembly.fasta".format(self.fb_reads_r1, self.fb_reads_r2, self.threads, self.outdir, spades_params)) if (self.fb_reads_r1 and self.fb_reads_r2) else None,
		]



def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("assembler", type=str, choices=["unicycler", "spades", "velvet"], default="unicycler")
	ap.add_argument("r1", type=str)
	ap.add_argument("r2", type=str)
	ap.add_argument("outdir", type=str)
	ap.add_argument("--threads", type=int, default=8)

	args = ap.parse_args()

	asm = None
	if args.assembler == "unicycler":
		asm = UnicyclerWrapper(args)
	
	if not asm is None:
		asm.run()		


	pass

if __name__ == "__main__": 
	main()

