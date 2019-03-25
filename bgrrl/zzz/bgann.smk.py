import sys
import csv
import os
from os.path import join, basename, dirname
from collections import Counter

from bgrrl.samplesheet import readSamplesheet
from bgrrl import TIME_CMD
from bgrrl.snakemake_helper import loadPreCmd

DEBUG = config.get("debugmode", False)
CWD = os.getcwd()

# tools
RATT_WRAPPER = join(config["etc"], "wrappers", "ratt_wrapper")

# setup i/o
OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
PROKKA_DIR = join(ANNOTATION_DIR, "prokka")
RATT_DIR = join(ANNOTATION_DIR, "ratt")

#INPUTFILES = set(row[0] for row in csv.reader(open(config["samplesheet"]), delimiter=","))
INPUTFILES = {row[0]: row[2] for row in csv.reader(open(config["samplesheet"]), delimiter=",")}

# generate target list
TARGETS = list()
if True:
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".prokka.gff"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".ffn"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".ffn.16S"), INPUTFILES))

	
REF_PREFIXES = dict()
if config["run_ratt"]:
	for d in next(os.walk(config["ratt_reference"]))[1]:     
		TARGETS.extend(map(lambda s:join(RATT_DIR, s, d, "{}_{}.final.gff".format(s, d)), INPUTFILES))
		TARGETS.append(join(config["ratt_reference"], d, "gff", d + ".gff"))
		REF_PREFIXES[d] = list(f.strip(".embl").split(".")[-1] for f in next(os.walk(join(config["ratt_reference"], d)))[2] if f.endswith(".embl"))


if DEBUG:
	print("CONFIG")
	print(config)

	with open("ann-targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

# helpers
def get_ref_index(wc):
	return REF_PREFIXES.get(wc.ref_ann)

def get_ref(wc):
	if len(REF_PREFIXES[wc.ref_ann]) > 1:
		return [join(config["ratt_reference"], wc.ref_ann, "gff", "{}.{}.parts_gff".format(wc.ref_ann, index)) for index in REF_PREFIXES[wc.ref_ann]]
	else:
		return [join(config["ratt_reference"], wc.ref_ann, "gff", "{}.parts_gff".format(wc.ref_ann))]

def get_assembly(wc):
    return INPUTFILES[wc.sample]

def get_singularity_call(cfg, cmd):
    assert "singularity_env" in cfg
    return "singularity exec {0} {1}".format(cfg.get("singularity_env"), cmd)


### RULES ###

localrules: all, ann_prokka_gffconvert, ann_prokka_16S

rule all:
	input: TARGETS

if True:
	rule ann_prokka:
		message:
			"Running de novo gene/functional annotation with prokka (incl. barrnap)..."
		input:
			# contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
			contigs = get_assembly
		output:
			log = join(PROKKA_DIR, "{sample}", "{sample}.log"),
			faa = join(PROKKA_DIR, "{sample}", "{sample}.faa"),
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn"),
			gff = join(PROKKA_DIR, "{sample}", "{sample}.gff")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
		params:
			outdir = lambda wildcards: join(PROKKA_DIR, wildcards.sample),
			prefix = lambda wildcards: wildcards.sample,
			load = loadPreCmd(config["load"]["prokka"]),
			centre = config["misc"]["seqcentre"],
			container = ("--singularity-container " + config["singularity_env"]) if "singularity_env" in config else "",
			custom_proteins = ("--proteins " + config["custom_prokka_proteins"]) if config.get("custom_prokka_proteins", "") else ""
		threads:
			8
		shell:
			# "set +u && source deactivate && source activate prokka_env" + \
			# " && " + \
			"prokka_wrapper {input.contigs} {params.container} --prefix {params.prefix} --outdir {params.outdir} --seq-centre {params.centre} --force --threads {threads}" + \
			" {params.custom_proteins}" + \
			" &> {log}"

	rule ann_prokka_gffconvert:
		message:
			"Renaming contig/scaffold names in prokka-generated gff..." # !TODO might need to do that with .gbk, too?
		input:
			gff = join(PROKKA_DIR, "{sample}", "{sample}.gff")
		output:
			gff = join(PROKKA_DIR, "{sample}", "{sample}.prokka.gff")
		shell:
			"awk -v OFS=\"\\t\" -v FS=\"\\t\"" + \
			" 'BEGIN {{ print \"##gff-version 3\"; }}" + \
			" /^[^#]/ {{ $1=gensub(\"[^>_]+_\", \"\", \"g\", $1); }}" + \
			" {{ if (NF > 1) print $0; }}' {input.gff} > {output.gff}"

	rule ann_prokka_16S:
		message:
			"Extracting 16S rRNA sequences from barrnap-annotation..."
		input:
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn")
		output:
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn.16S")
		params:
			# load = loadPreCmd(config["load"]["seqtk"])
			seqtk = get_singularity_call(config, "seqtk")
		shell:
			# "{params.load} " + \
			"{params.seqtk} subseq {input.ffn} <(grep -o \">[^ ]\+ 16S ribosomal RNA\" {input.ffn} | cut -f 1 -d \" \" | cut -f 2 -d \>) > {output.ffn}"



if config["run_ratt"]:
	rule ann_ratt_prepref:
		message:
			"Preprocessing reference annotations for ratt-annotation transfer..."
		input:
			embl_ref = join(config["ratt_reference"], "{ref_ann}", "{prefix}.embl")
		output:
			join(config["ratt_reference"], "{ref_ann}", "gff", "{prefix}.parts_gff")
		# log:
		#	join(config["cwd"], ANNOTATION_DIR, "log", "{ref_ann}.ann_ratt_prepref.log")
		params:
			load = loadPreCmd(config["load"]["ratt"]),
			refdir = join(config["cwd"], config["ratt_reference"], "{ref_ann}")				
		threads:
			1
		shell:
			#"{params.load}" + \
			"set +u && " + \
			"source deactivate && source activate ratt_env && " + \
			"cd {params.refdir} && " + \
			"seqret -sequence {input.embl_ref} -outseq $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff -offormat gff -feature && " + \
			"rm $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff && " + \
			"mkdir -p gff && mv $(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff gff/$(basename {input.embl_ref} .embl | sed 's/_x_/./g').parts_gff && cd " + CWD 
			



			#" set +u && source ratt-dev_no_defined && source deactivate && source activate ratt_env &&" + \
			#" cd {params.refdir}" + \
			#" && seqret -sequence {input.embl_ref} -outseq $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff -offormat gff -feature" + \
			#" && rm $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff" + \
			#" && mkdir -p gff && mv $(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff gff/$(basename {input.embl_ref} .embl | sed 's/_x_/./g').parts_gff && cd " + CWD 

	rule ann_ratt_mergegff:
		message:
			"Merging single-gene .gffs from ratt-annotation transfer..."
		input:
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
			"cat {input.gff} > {output}"

	rule ann_ratt:
		message:
			"Running annotation-transfer with ratt..."
		input:
			#contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta"),
			contigs = get_assembly,
			reference = join(config["ratt_reference"], "{ref_ann}")
		output:
			done = join(RATT_DIR, "{sample}", "{ref_ann}", "{sample}_{ref_ann}.final.gff")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_{ref_ann}.ann_ratt.log")
		params:	
			outdir = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann),
			load = loadPreCmd(config["load"]["ratt"]),
			done = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann, wildcards.sample + "_" + wildcards.ref_ann + ".done"),
			path_to_ratt = "/tgac/software/testing/ratt/dev_no_defined/x86_64/bin"
		threads:
			1
		shell:
			#export RATT_HOME=/tgac/software/testing/ratt/dev_no_defined/x86_64/bin/
			#export RATT_CONFIG=$RATT_HOME/RATT.config_bac
			"set +u && source deactivate && source activate ratt_env && source ratt-dev_no_defined && " + \
			"ratt_wrapper -o {params.outdir} " + \
			join(config["cwd"], "{input.contigs}") + " " + \
			"{input.reference} {params.path_to_ratt} &> {log}"

			# ratt_wrapper [-h] [--outdir OUTDIR] contigs reference path_to_ratt
			#RATT_WRAPPER + " {params.outdir} " + join(config["cwd"], "{input.contigs}") + \
			#" {input.reference} {params.done} &> {log}"
