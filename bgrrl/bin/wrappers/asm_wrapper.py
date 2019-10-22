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
		
		self.ur1, self.ur2, self.merged, self.singles, self.uc_singles = args.ur1, args.ur2, args.merged, args.singles, args.uc_singles

		#r1_reads, r2_reads = args.r1.split(","), args.r2.split(",")
			
		#self.reads_r1 = os.path.abspath(r1_reads.pop(0))
		#self.fb_reads_r1 = os.path.abspath(r1_reads.pop(0)) if r1_reads else ""
		
		#if not check_reads(self.reads_r1) and not check_reads(self.fb_reads_r1):
		#	raise ValueError("Cannot locate r1 reads.")

		#if r2_reads:
		#	self.reads_r2 = os.path.abspath(r2_reads.pop(0))
		#	self.fb_reads_r2 = os.path.abspath(r2_reads.pop(0)) if r2_reads else ""

		#	if not check_reads(self.reads_r2) and not check_reads(self.fb_reads_r2):
		#		raise ValueError("Cannot locate r2 reads.")
		#else:
		#	self.reads_r2, self.fb_reads_r2 = None, None
		
		self.outdir = args.outdir
		self.threads = args.threads

		self.singularity_prefix = ("singularity exec " + args.singularity_container + " ") if args.singularity_container else ""
		self.fallback_queue = list()


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

		cmd = "cd {0}" + \
			" && " + self.singularity_prefix + \
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


class SpadesSCWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)

		spades_params = "--careful --cov-cutoff auto --sc"
		spades_cmd = self.singularity_prefix + "spades.py -s {0} -t {1} -o {2} -m 64 {3} && ln -s scaffolds.fasta {2}/assembly.fasta"

		if not self.fb_reads_r1:
			self.fallback.queue.extend(
				[
					("asm_main_scspt", spades_cmd.format(self.reads_r1, self.threads, self.outdir, spades_params))
				]
			)
		else:
			self.fallback_queue.extend(
				[
					("asm_main_scspn", spades_cmd.format(self.reads_r1, self.threads, self.outdir, spades_params)),
					("asm_fb1_scspt", spades_cmd.format(self.fb_reads_r1, self.threads, self.outdir, spades_params))
				]
			)
		
						
		
		

		
class UnicyclerWrapper(ASM_Wrapper):
	def __init__(self, args):
		super().__init__(args)

		spades_params = "--careful --cov-cutoff auto"
		unicycler_params = "--min_polish_size 1000 --no_correct"

		unicycler_cmd = self.singularity_prefix + "unicycler -1 {0} -2 {1} -s {2} -t {3} -o {4} {5}"
		spades_cmd = self.singularity_prefix + "spades.py -1 {0} -2 {1} --merged {2}  -s {3} -t {4} -o {5} -m 64 {6} && ln -s scaffolds.fasta {5}/assembly.fasta" 

		self.fallback_queue = [(
			"asm_main_uc",
			unicycler_cmd.format(
				self.ur1,
				self.ur2,
				self.uc_singles,
				self.threads,
				self.outdir,
				unicycler_params
			)), 
			(
				"asm_fb1_sp",
				spades_cmd.format(
					self.ur1,
					self.ur2,
					self.merged,
					self.singles,
					self.threads,
					self.outdir,
					spades_params
			))		
		]



def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("assembler", type=str, choices=["unicycler", "spades", "velvet", "spades_sc"], default="unicycler")
	ap.add_argument("outdir", type=str)
	ap.add_argument("ur1", type=str)
	ap.add_argument("ur2", type=str)
	ap.add_argument("merged", type=str)
	ap.add_argument("singles", type=str)
	ap.add_argument("uc_singles", type=str)
	ap.add_argument("--threads", type=int, default=8)
	ap.add_argument("--singularity-container", type=str, default="")

	args = ap.parse_args()
	#if not args.r1:
	#	raise ValueError("Missing reads input")

	asm = None
	if args.assembler == "unicycler":
		asm = UnicyclerWrapper(args)
	elif args.assembler == "velvet":
		asm = VelvetWrapper(args)
	elif args.assembler == "spades_sc":
		asm = SpadesSCWrapper(args)
	elif args.assembler == "spades":
		raise ValueError("spades not supported yet")
	
	if not asm is None:
		asm.run()		


	pass

if __name__ == "__main__": 
	main()

