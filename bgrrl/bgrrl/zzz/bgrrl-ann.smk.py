import sys
import csv
import os
from os.path import join, basename, dirname
from collections import Counter

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
PROKKA_DIR = join(ANNOTATION_DIR, "prokka")
RATT_DIR = join(ANNOTATION_DIR, "ratt")

PROKKA_WRAPPER = join(config["etc"], "wrappers", "prokka_wrapper")
RATT_WRAPPER = join(config["etc"], "wrappers", "ratt_wrapper")

# INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
# with open("ann-inputfiles.txt", "w") as input_out:
#     print(*INPUTFILES.values(), sep="\n", file=input_out)
INPUTFILES = set(row[0] for row in csv.reader(open(config["samplesheet"]), delimiter=","))

TARGETS = list()
REF_PREFIXES = dict()
if config["run_prokka"]:
	# TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".log"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".prokka.gff"), INPUTFILES))
if config["run_ratt"]:
	for d in next(os.walk(config["ratt_reference"]))[1]:     
		# TARGETS.extend(map(lambda s:join(RATT_DIR, s, d, "{}_{}.done".format(s, d)), INPUTFILES))
		TARGETS.extend(map(lambda s:join(RATT_DIR, s, d, "{}_{}.final.gff".format(s, d)), INPUTFILES))
		REF_PREFIXES[d] = list(f.strip(".embl").split(".")[-1] for f in next(os.walk(join(config["ratt_reference"], d)))[2] if f.endswith(".embl"))
		#for f in next(os.walk(join(config["ratt_reference"], d)))[2]:
		#	if f.endswith(".embl"):
		#		REF_PREFIXES.append(basename(f))
		#		TARGETS.append(join(config["ratt_reference"], d, "gff", f.replace(".embl", ".gff")))
		TARGETS.append(join(config["ratt_reference"], d, "gff", d + ".gff"))

CWD = os.getcwd()

print("CONFIG")
print(config)

with open("ann-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

def get_ref_index(wc):
	return REF_PREFIXES.get(wc.ref_ann)


def get_ref(wc):
	if len(REF_PREFIXES[wc.ref_ann]) > 1:
		return [join(config["ratt_reference"], wc.ref_ann, "gff", "{}.{}.parts_gff".format(wc.ref_ann, index)) for index in REF_PREFIXES[wc.ref_ann]]
	else:
		return [join(config["ratt_reference"], wc.ref_ann, "gff", "{}.parts_gff".format(wc.ref_ann))]


localrules: all, ann_prokka_gffconvert

rule all:
	input: TARGETS

if config["run_prokka"]:
	rule ann_prokka:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		output:
			log = join(PROKKA_DIR, "{sample}", "{sample}.log"),
			faa = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.faa"),
			ffn = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.ffn"),
			gff = join(PROKKA_DIR, "{sample}", "{sample}.gff")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
		params:
			outdir = lambda wildcards: join(PROKKA_DIR, wildcards.sample),
			prefix = lambda wildcards: wildcards.sample,
			load = loadPreCmd(config["load"]["prokka"]),
			centre = config["seqcentre"]
		threads:
			8
		shell:
			"cd {params.outdir}" + \
			" && echo $(pwd)" + \
			" && cp -v " + join(config["cwd"], "{input.contigs}") + " prokka_contigs.fasta" + \
			" && ls -l" + \
			" && {params.load} (" + TIME_CMD + \				
			" prokka --cpus {threads} --outdir . --prefix {params.prefix} --centre {params.centre} prokka_contigs.fasta --force" + \
			" && rm prokka_contigs.fasta" + \
			" && cd " + CWD + ") &> {log}"

	rule ann_prokka_gffconvert:
		input:
			gff = join(PROKKA_DIR, "{sample}", "{sample}.gff")
		output:
			gff = join(PROKKA_DIR, "{sample}", "{sample}.prokka.gff")
		shell:
			"awk -v OFS=\"\t\" -v FS=\"\t\"" + \
			" 'BEGIN {{ print \"##gff-version 3\"; }}" + \
			" /^[^#]/ {{ $1=gensub(\"[^>_]+_\", \"\", \"g\", $1); }}" + \
			" {{ if (NF > 1) print $0; }}' {input.gff} > {output.gff}"



if config["run_ratt"]:

	if True:

		rule ann_ratt_prepref:
			input:
				embl_ref = join(config["ratt_reference"], "{ref_ann}", "{prefix}.embl")
			output:
				join(config["ratt_reference"], "{ref_ann}", "gff", "{prefix}.parts_gff")
			log:
				join(config["cwd"], ANNOTATION_DIR, "log", "{ref_ann}.ann_ratt_prepref.log")
			params:
				load = loadPreCmd(config["load"]["ratt"]),
				refdir = join(config["cwd"], config["ratt_reference"], "{ref_ann}")				
			threads:
				1
			shell:
				"{params.load}" + \
				" cd {params.refdir}" + \
				" && seqret -sequence {input.embl_ref} -outseq $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff -offormat gff -feature" + \
				" && rm $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff" + \
				" && mkdir -p gff && mv $(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff gff/$(basename {input.embl_ref} .embl | sed 's/_x_/./g').parts_gff && cd " + CWD 


		rule ann_ratt_mergegff:
			input:
				# gff = expand(join(config["ratt_reference"], "{ref_ann}", "gff", "{ref_ann}.{index}.parts_gff"), ref_ann=next(os.walk(config["ratt_reference"]))[1], index=get_ref_index)
				gff = get_ref
			output:
				join(config["ratt_reference"], "{ref_ann}", "gff", "{ref_ann}.gff")
			log:
				join(config["cwd"], ANNOTATION_DIR, "log", "{ref_ann}.ann_ratt_mergegff.log")
			params:
				gffdir = join(config["cwd"], config["ratt_reference"], "{ref_ann}", "gff")
			threads:
				1
			shell:
				# "cd {params.gffdir}" + \
				"cat {input.gff} > {output}"

	rule ann_ratt:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta"),
			reference = join(config["ratt_reference"], "{ref_ann}")
		output:
			# done = join(RATT_DIR, "{sample}", "{ref_ann}", "{sample}_{ref_ann}.done")
			done = join(RATT_DIR, "{sample}", "{ref_ann}", "{sample}_{ref_ann}.final.gff")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_{ref_ann}.ann_ratt.log")
		params:	
			outdir = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann),
			load = loadPreCmd(config["load"]["ratt"]),
			done = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann, wildcards.sample + "_" + wildcards.ref_ann + ".done")
		threads:
			1
		shell:
			RATT_WRAPPER + " {params.outdir} " + join(config["cwd"], "{input.contigs}") + \
			" {input.reference} {params.done} &> {log}"
			# " touch " + join(config["cwd"], "{output.done}") + ")" 
			# " && extractfeat -sequence *.final.embl -outseq {params.prefix}.fasta -type CDS" + \
			# " && transeq -sequence {params.prefix}.fasta {params.prefix}.pep.fasta" + \
			# " && mv {params.prefix}.fasta {params.prefix}.ffn" + \
			# " && mv {params.prefix}.pep.fasta {params.prefix}.faa" + \
