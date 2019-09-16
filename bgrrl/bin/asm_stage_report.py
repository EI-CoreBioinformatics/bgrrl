import sys
import os
from os.path import join, basename, dirname
import argparse

from collections import Counter, OrderedDict


ASSEMBLY_STAGES = OrderedDict([
	("asm_main_ucn", "Main,Unicycler,normalized"),
	("asm_fb1_uct", "Fallback1,Unicycler,trimmed"),
	("asm_fb2_spn", "Fallback2,Spades,normalized"),
	("asm_fb3_spt", "Fallback3,Spades,trimmed"),
	("asm_main_uct", "Main,Unicycler,trimmed"),
	("asm_fb1_spt", "Main,Spades,trimmed"),
	("asm_main_ven", "Main,Velvet,normalized"),
	("asm_fb1_vet", "Fallback1,Velvet,trimmed"),
	("asm_main_vet", "Main,Velvet,trimmed"),
	("asm_main_scspn", "Main,SpadesSC,normalized"),
	("asm_fb1_scspt", "Fallback1,SpadesSC,trimmed"),
	("asm_main_scspt", "Main,SpadesSC,trimmed"),
	("NA", "not_assembled")])


def compileASMInfo(asm_dir, out=sys.stdout, asm_stat_out=sys.stdout, asm_samplesheet=sys.stdout):
	asm_tag_ctr = Counter()
	for cdir, dirs, files in os.walk(asm_dir):
		sample = basename(cdir)
		if sample != "log" and os.path.dirname(cdir) == asm_dir:
			asm_tag = ([f for f in files if f.startswith("asm_")] + ["NA"])[0]
			asm_tag_ctr[asm_tag] += 1 # .setdefault(asm_tag, list()).append(cdir)
			print(sample, asm_tag, sep="\t", file=out)
			print(sample, sample, join(cdir, sample + ".assembly.fasta"), sep=",")
			print(sample, sample, join(cdir, sample + ".assembly.fasta"), sep=",", file=asm_samplesheet)
	print(asm_tag_ctr)
	for asm_tag in ASSEMBLY_STAGES:
		if asm_tag_ctr[asm_tag] > 0:
			print(ASSEMBLY_STAGES[asm_tag], asm_tag_ctr[asm_tag], asm_tag_ctr[asm_tag]/sum(asm_tag_ctr.values()), sep="\t", file=asm_stat_out)
	print("Total", "", "", sum(asm_tag_ctr.values()), sep="\t", file=asm_stat_out)


def main(args):

	ap = argparse.ArgumentParser()
	ap.add_argument("indir", type=str, default=".")
	ap.add_argument("report_dir", type=str, default=".")
	args = ap.parse_args(args)

	print("Gathering assembly stage information...")
	with open(join(args.report_dir, "assembly_stage_summary.tsv"), "wt") as asm_stage_summary, open(join(args.report_dir, "assembly_stages.tsv"), "wt") as asm_stages, open(join(args.report_dir, "samplesheets", "samplesheet.asm_pass.tsv"), "wt") as asm_samplesheet:
		compileASMInfo(join(args.indir, "assembly"), out=asm_stages, asm_stat_out=asm_stage_summary, asm_samplesheet=asm_samplesheet)
	return True

	

