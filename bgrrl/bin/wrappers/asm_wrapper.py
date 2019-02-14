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
		def check_reads(reads):
			return not reads or os.path.exists(reads)

		r1_reads, r2_reads = args.r1.split(","), args.r2.split(",")
	
		self.reads_r1 = os.path.abspath(r1_reads.pop(0))
		self.fb_reads_r1 = os.path.abspath(r1_reads.pop(0)) if r1_reads else ""
		
		if not check_reads(self.reads_r1) and not check_reads(self.fb_reads_r1):
			raise ValueError("Cannot locate r1 reads.")

		self.reads_r2 = os.path.abspath(r2_reads.pop(0))
		self.fb_reads_r2 = os.path.abspath(r2_reads.pop(0)) if r2_reads else ""

		if not check_reads(self.reads_r2) and not check_reads(self.fb_reads_r2):
			raise ValueError("Cannot locate r2 reads.")
		
		self.outdir = args.outdir
		self.threads = args.threads


	def run(self, env="bgasm_env"):
		self.fallback_queue = [item for item in self.fallback_queue if not item is None]

		assembly_complete = False
		# conda_cmd = "source deactivate && source activate {} && ".format(env)
		conda_cmd = ""
		
		while self.fallback_queue:
			step, cmd = self.fallback_queue.pop(0)
			print("BGRR|_ASM_WRAPPER::" + step)
			
			try:
				out = subprocess.check_output(conda_cmd + cmd, stderr=subprocess.STDOUT, shell=True)
			except subprocess.CalledProcessError as err:
				print(err, flush=True)
				print("Stage: {} failed.".format(step), flush=True)
				shutil.rmtree(self.outdir, ignore_errors=True)
				continue

			print(out)

			open(os.path.join(self.outdir, step), "wt").close()
			assembly_complete = True
			break

		if not assembly_complete:
			pathlib.Path(self.outdir).mkdir(parents=True, exist_ok=True)
			open(os.path.join(self.outdir, "assembly.fasta"), "wt").close()		









class VelvetWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)
		
		tmpdir = os.path.join(self.outdir, "tmp")
		pathlib.Path(tmpdir).mkdir(parents=True, exist_ok=True)

		cwd = os.getcwd()

		singularity_prefix = ("singularity exec " + args.singularity_container + " ") if args.singularity_container else ""

		cmd = "cd {0}" + \
			" && " + singularity_prefix + \
			"VelvetOptimiser.pl -v -s 19 -e 191 -x 2 -t {1} -d ../velvet_assembly -f '-shortPaired -fastq.gz -separate {2} {3}'" + \
			" && cd .." + \
			" && mv -v velvet_assembly/* ." + \
			" && rm -rf velvet_assembly/ tmp/" + \
			" && ln -s contigs.fa assembly.fasta" + \
			" && cd {4}"
	
		if not (self.fb_reads_r1 and self.fb_reads_r2):
			self.fallback_queue = [
				(
					"asm_main_vet", 
					cmd.format(
						tmpdir, 
						self.threads, 
						os.path.join(cwd, self.reads_r1), 
						os.path.join(cwd, self.reads_r2), 
						cwd
					)
				)
			]
		else:
			self.fallback_queue = [
				(
					"asm_main_ven", 
					cmd.format(
						tmpdir, 
						self.threads, 
						os.path.join(cwd, self.reads_r1), 
						os.path.join(cwd, self.reads_r2), 
						cwd
					)
				),
				(
					"asm_fb1_vet", 
					cmd.format(
						tmpdir, 
						self.threads, 
						os.path.join(cwd, self.fb_reads_r1), 
						os.path.join(cwd, self.fb_reads_r2), 
						cwd
					)
				)
			]
		
class UnicyclerWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)

		spades_params = "--careful --cov-cutoff auto"
		unicycler_params = "--min_polish_size 1000"

		singularity_prefix = ("singularity exec " + args.singularity_container + " ") if args.singularity_container else ""

		unicycler_cmd = singularity_prefix + "unicycler -1 {0} -2 {1} -t {2} -o {3} {4}"
		spades_cmd = singularity_prefix + "spades.py -1 {0} -2 {1} -t {2} -o {3} -m 64 {4} && ln -s scaffolds.fasta {3}/assembly.fasta" 


		self.fallback_queue = list()
		if not (self.fb_reads_r1 and self.fb_reads_r2):
			self.fallback_queue = [
				(
					"asm_main_uct", 
					unicycler_cmd.format(
						self.reads_r1, 
						self.reads_r2, 
						self.threads, 
						self.outdir, 
						unicycler_params
					)
				),
				(
					"asm_fb1_spt", 
					spades_cmd.format(
						self.reads_r1, 
						self.reads_r2, 
						self.threads, 
						self.outdir, 
						spades_params
					)
				)
			]
		else:
			self.fallback_queue = [
				(
					"asm_main_ucn", 
					unicycler_cmd.format(
						self.reads_r1, 
						self.reads_r2, 
						self.threads, 
						self.outdir, 
						unicycler_params
					)
				),
				(
					"asm_fb1_uct", 
					unicycler_cmd.format(
						self.fb_reads_r1, 
						self.fb_reads_r2, 
						self.threads, 
						self.outdir, 
						unicycler_params
					)
				),
				(
					"asm_fb2_spn", 
					spades_cmd.format(
						self.reads_r1, 
						self.reads_r2, 
						self.threads, 
						self.outdir, 
						spades_params
					)
				),
				(
					"asm_fb3_spt", 
					spades_cmd.format(
						self.fb_reads_r1, 
						self.fb_reads_r2, 
						self.threads, 
						self.outdir, 
						spades_params
					)
				)
			]



def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("assembler", type=str, choices=["unicycler", "spades", "velvet"], default="unicycler")
	ap.add_argument("r1", type=str)
	ap.add_argument("r2", type=str)
	ap.add_argument("outdir", type=str)
	ap.add_argument("--threads", type=int, default=8)
	ap.add_argument("--singularity-container", type=str, default="")

	args = ap.parse_args()

	asm = None
	if args.assembler == "unicycler":
		asm = UnicyclerWrapper(args)
	elif args.assembler == "velvet":
		asm = VelvetWrapper(args)
	
	if not asm is None:
		asm.run()		


	pass

if __name__ == "__main__": 
	main()

