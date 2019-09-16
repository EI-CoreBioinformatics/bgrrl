import sys
import csv
import os
from os.path import join, basename, dirname
from collections import Counter

from bgrrl.samplesheet import readSamplesheet, Samplesheet, ASM_Sample, ANN_Sample
from bgrrl.snakemake_helper import get_cmd_call

DEBUG = config.get("debugmode", False)
SKIP_NORMALIZATION = config.get("no_normalization", False)
SINGLE_CELL_MODE = config.get("single_cell_mode", False)

TIME_V = config.get("tools", dict()).get("time", "time")

# set up i/o
OUTPUTDIR = os.path.abspath(config["out_dir"])
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
QC_OUTDIR = join(OUTPUTDIR, "qc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
PROKKA_DIR = join(ANNOTATION_DIR, "prokka")
RATT_DIR = join(ANNOTATION_DIR, "ratt")

# tools
RATT_WRAPPER = join(config["etc"], "wrappers", "ratt_wrapper")

if config["module"] == "bgasm":
	sampletype = ASM_Sample
	
elif config["module"] == "bgann":
	sampletype = ANN_Sample
	# INPUTFILES = {row[0]: row[2] for row in csv.reader(open(config["samplesheet"]), delimiter=",")}
else:
	raise ValueError("Module not recognized as bgasm/bgann: " + config["module"])

INPUTFILES = Samplesheet(config["samplesheet"], sampletype=sampletype)

# generate target list
TARGETS = list()

if config["run_prokka"]:
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".gff"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".fna"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".ffn.16S"), INPUTFILES))
else:
	TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))
	
REF_PREFIXES = dict()
if config["run_ratt"]:
	RATT_REF_PATH = os.path.abspath(config.get("ratt_reference", "."))
	
	for d in next(os.walk(RATT_REF_PATH))[1]:     
		TARGETS.extend(map(lambda s:join(RATT_DIR, s, d, "{}_{}.final.gff".format(s, d)), INPUTFILES))
		TARGETS.append(join(RATT_REF_PATH, d, "gff", d + ".gff"))
		REF_PREFIXES[d] = list(f.strip(".embl").split(".")[-1] for f in next(os.walk(join(RATT_REF_PATH, d)))[2] if f.endswith(".embl"))

# define normalization/no-normalization behaviour
if SKIP_NORMALIZATION:
	PRIMARY_READDIR = BBDUK_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbduk"
	SECONDARY_READ_ID = "bbduk"
else:
	PRIMARY_READDIR = BBNORM_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbnorm"
	SECONDARY_READ_ID = "bbduk"

if DEBUG:
	print("CONFIG")
	print(config)
	with open("asm-inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("asm-targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)


# helpers
def get_sample_files(wc):
	s = INPUTFILES[wc.sample]
	# default case is with normalization, i.e. if --no-normalization option isn't present, it should be "False"
	return (get_r1(wc).split(",") + get_r2(wc).split(",")) if not SINGLE_CELL_MODE else get_r1(wc).split(",")
	if SKIP_NORMALIZATION:
		return s.R1trim, s.R2trim
	else:
		return s.R1norm, s.R1trim, s.R2norm, s.R2trim

def get_r1(wc):
	s = INPUTFILES[wc.sample]
	return ",".join((s.R1norm, s.R1trim)) if not SKIP_NORMALIZATION else s.R1trim

def get_r2(wc):
	s = INPUTFILES[wc.sample]
	return ",".join((s.R2norm, s.R2trim)) if not SKIP_NORMALIZATION else s.R2trim

def get_ref_index(wc):
	return REF_PREFIXES.get(wc.ref_ann)

def get_ref(wc):
	if len(REF_PREFIXES[wc.ref_ann]) > 1:
		return [join(RATT_REF_PATH, wc.ref_ann, "gff", "{}.{}.parts_gff".format(wc.ref_ann, index)) for index in REF_PREFIXES[wc.ref_ann]]
	else:
		return [join(RATT_REF_PATH, wc.ref_ann, "gff", "{}.parts_gff".format(wc.ref_ann))]

def get_assembly(wc):
    return INPUTFILES[wc.sample]

CMD_CALL = get_cmd_call(config, "bgrrl_container")
CONTAINER_PARAM = ("--singularity-container " + " ".join(CMD_CALL.split(" ")[2:])) if CMD_CALL else ""


### RULES ###

if config["run_prokka"]:
	localrules: all, ann_prokka_gffconvert, ann_prokka_16S
else:
	localrules: all

rule all:
	input: TARGETS

if config["module"] == "bgasm":
	rule asm_assembly:
		message:
			"Generating assembly with " + config["assembler"] + "..."
		input:
			get_sample_files
		output:
			join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
		params:
			outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
			assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "assembly.fasta"),
			assembler = config["assembler"] if not SINGLE_CELL_MODE else "spades_sc",
			reads = lambda wildcards: (get_r1(wildcards) + " " + get_r2(wildcards)) if not SINGLE_CELL_MODE else get_r1(wildcards),
			r1 = lambda wildcards: get_sample_files(wildcards)[0] if SKIP_NORMALIZATION else ",".join(get_sample_files(wildcards)[:2]),
			r2 = lambda wildcards: get_sample_files(wildcards)[1] if SKIP_NORMALIZATION else ",".join(get_sample_files(wildcards)[2:]),
			container = CONTAINER_PARAM
		threads:
			8
		shell:
			"asm_wrapper --threads {threads} {params.container} {params.assembler} {params.outdir} {params.reads} &> {log}"

	rule asm_postprocess:
		message:
			"Postprocessing assembly (minimum contig length={})...".format(config["asm_lengthfilter_contig_minlen"]) 
		input:
			assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
		output:
			filtered_assembly = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_postprocess.log")
		params:
			minlen = int(config["asm_lengthfilter_contig_minlen"]),
			cmd = CMD_CALL + "reformat.sh"
		threads:
			1
		shell:
			TIME_V + " {params.cmd}" + \
			" in={input.assembly} out={output[0]} minlength={params.minlen}"


if config["run_prokka"]:
	rule ann_prokka:
		message:
			"Running de novo gene/functional annotation with prokka (incl. barrnap)..."
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		output:
			log = join(PROKKA_DIR, "{sample}", "{sample}.log"),
			faa = join(PROKKA_DIR, "{sample}", "{sample}.faa"),
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn"),
			gff = join(PROKKA_DIR, "{sample}", "{sample}.gff"),
			fna = join(PROKKA_DIR, "{sample}", "{sample}.fna")
		log:
			join(ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
		params:
			outdir = lambda wildcards: join(PROKKA_DIR, wildcards.sample),
			prefix = lambda wildcards: wildcards.sample,
			centre = config["misc"]["seqcentre"],
			container = CONTAINER_PARAM,
			custom_proteins = ("--proteins " + config["custom_prokka_proteins"]) if config.get("custom_prokka_proteins", "") else ""
		threads:
			8
		shell:
			"prokka_wrapper {input.contigs} {params.container} --prefix {params.prefix} --outdir {params.outdir} --seq-centre {params.centre} --force --threads {threads}" + \
            " {params.custom_proteins} && " + \
			" touch {PROKKA_DIR}/{wildcards.sample}/{wildcards.sample}.txt && " + \
			" sed -i \"s/strain/{wildcards.sample}/\" {PROKKA_DIR}/{wildcards.sample}/{wildcards.sample}.txt" + \
			" &> {log}"

	rule ann_prokka_16S:
		message:
			"Extracting 16S rRNA sequences from barrnap-annotation..."
		input:
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn")
		output:
			ffn = join(PROKKA_DIR, "{sample}", "{sample}.ffn.16S")
		params:
			cmd = CMD_CALL + "seqtk"
		shell:
			"{params.cmd} subseq {input.ffn} <(grep -o \">[^ ]\+ 16S ribosomal RNA\" {input.ffn} | cut -f 1 -d \" \" | cut -f 2 -d \>) > {output.ffn}"


if config["run_ratt"]:
	rule ann_ratt_prepref:
		message:
			"Preprocessing reference annotations for ratt-annotation transfer..."
		input:
			embl_ref = join(RATT_REF_PATH, "{ref_ann}", "{prefix}.embl")
		output:
			join(RATT_REF_PATH, "{ref_ann}", "gff", "{prefix}.parts_gff")
		params:
			refdir = join(RATT_REF_PATH, "{ref_ann}"),
			cmd = CMD_CALL + "seqret"			
		threads:
			1
		shell:
			"cd {params.refdir} && " + \
			"{params.cmd} -sequence $(basename {input.embl_ref}) -outseq $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff -offormat gff -feature && " + \
			"rm $(basename {input.embl_ref} .embl).$(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff && " + \
			"mkdir -p gff && mv $(basename {input.embl_ref} .embl | sed 's/\./_x_/g').gff gff/$(basename {input.embl_ref} .embl | sed 's/_x_/./g').parts_gff && cd -" 

	rule ann_ratt_mergegff:
		message:
			"Merging single-gene .gffs of ratt-annotation transfer reference..."
		input:
			gff = get_ref
		output:
			join(config["ratt_reference"], "{ref_ann}", "gff", "{ref_ann}.gff")
		log:
			join(ANNOTATION_DIR, "log", "{ref_ann}.ann_ratt_mergegff.log")
		params:
			gffdir = join(RATT_REF_PATH, "{ref_ann}", "gff")
		threads:
			1
		shell:
			"cat {input.gff} > {output}"

	rule ann_ratt:
		message:
			"Running annotation-transfer with ratt..."
		input:
			contigs = rules.ann_prokka.output.fna,# join(PROKKA_DIR, "{sample}", "{sample}.ffn"),
			reference = join(RATT_REF_PATH, "{ref_ann}")
		output:
			done = join(RATT_DIR, "{sample}", "{ref_ann}", "{sample}_{ref_ann}.final.gff")
		log:
			join(ANNOTATION_DIR, "log", "{sample}_{ref_ann}.ann_ratt.log")
		params:	
			outdir = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann),
			container = CONTAINER_PARAM, 
		threads:
			1
		shell:
			"ratt_wrapper {params.container} -o {params.outdir} {input.contigs} {input.reference} some_dummy_path &> {log}"
