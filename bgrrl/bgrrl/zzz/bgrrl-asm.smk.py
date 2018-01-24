import sys
import csv
import os
import glob
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

SOFTWAREPATH = "/tgac/software/testing"
BBSUITE_DIR = os.path.join(SOFTWAREPATH, "bbmap", "37.24", "bbmap")
ADAPTERS = os.path.join(BBSUITE_DIR, "resources", "adapters.fa")

# tools
UNICYCLER_WRAPPER = join(config["etc"], "wrappers", "unicycler_wrapper")
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"
# needed for BUSCO
CWD = os.getcwd()

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
QA_DIR = join(OUTPUTDIR, "qa")

QC_OUTDIR = os.path.join(OUTPUTDIR, "qc")
FASTQC_DIR = os.path.join(QC_OUTDIR, "fastqc")
BBDUK_DIR = os.path.join(QC_OUTDIR, "bbduk")
KAT_DIR = os.path.join(QC_OUTDIR, "kat")
BBNORM_DIR = os.path.join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = os.path.join(QC_OUTDIR, "tadpole")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("asm-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(QA_DIR, "quast", s, "quast.log"), INPUTFILES))
# TARGETS.extend(map(lambda s:join(config["cwd"], QA_DIR, "busco", s, "run_geno", "short_summary_geno.txt"), INPUTFILES))
TARGETS.extend(map(lambda s:join(config["cwd"], QA_DIR, "busco", s, "short_summary_{}.txt".format(s)), INPUTFILES))
TARGETS.extend(map(lambda s:join(QA_DIR, "blobtools", "blob", s, s + ".blobDB.table.txt"), INPUTFILES))

# TARGETS.extend(map(lambda s:join(QA_DIR, "blobtools", "bwa", s, s + ".blob_bwa.bam"), INPUTFILES))
# TARGETS.extend(map(lambda s:join(QA_DIR, "blobtools", "blast", s, s + ".blob_blast.tsv"), INPUTFILES))
#if config.get("use_asm_lengthfilter", ""):
# TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))

print("CONFIG")
print(config)

print(config.get("use_asm_length_filter", "NA"))

with open("asm-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all, asm_link_assembly

rule all:
	input: TARGETS

rule asm_assembly:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
		assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "assembly.fasta"),
		final_assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, wildcards.sample + ".assembly.fasta")
		# load = loadPreCmd(config["load"]["unicycler"])
	threads:
		8
	shell:
		UNICYCLER_WRAPPER + \
		" -1 {input.r1} -2 {input.r2} -t {threads} -o {params.outdir} &> {log}" 
		# + \
		# " && ln -s {params.assembly} {params.final_assembly} &> {log}"

"""
if 0:
	rule asm_link_assembly:
		input:
			assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
		output:
			link = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		# params:
		#	linktarget = lambda wildcards: os.path.basename(input.assembly) 
		shell:
			"ln -s assembly.fasta {output.link}"

else: 
	#config.get("use_asm_lengthfilter", ""):
"""
rule asm_lengthfilter:
	input:
		assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	output:
		filtered_assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_lengthfilter.log")
	params:
		minlen = int(config["asm_lengthfilter_contig_minlen"])
	threads:
		1
	run:
		import shutil
		from ktio.ktio import readFasta
		with open(output[0], "w") as seqout:
			for _id, _seq in readFasta(input.assembly):
				if len(_seq) >= params.minlen:
					print(_id, _seq, sep="\n", file=seqout)
		# original_assembly = join(os.path.dirname(input.assembly), "assembly.fasta")
		# shutil.copy2(original_assembly, original_assembly + ".full")
		# shutil.copy2(output[0], wildcards.sample + "." + original_assembly)


rule qa_busco_geno:
	input:
		assembly = join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		# join(config["cwd"], QA_DIR, "busco", "{sample}", "run_geno", "short_summary_geno.txt")
		join(config["cwd"], QA_DIR, "busco", "{sample}", "short_summary_{sample}.txt")
	log:	
		join(config["cwd"], QA_DIR, "log", "{sample}_busco_geno.log")
	params:
		#outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", wildcards.sample, "run_geno"),
		outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "run_" + wildcards.sample),
		final_outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", wildcards.sample),
		tmp = lambda wildcards: join(config["cwd"], QA_DIR, "busco",  "tmp", wildcards.sample),
		load = loadPreCmd(config["load"]["busco"])
	threads:
		8
	shell:
		BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
		" {params.load}" + TIME_CMD + \
		" run_BUSCO.py -i {input.assembly} -c {threads} -m geno" + \
		" --force -t {params.tmp} -l " + BUSCO_DATA + " -o {wildcards.sample} &> {log} && cd " + CWD + \
		" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
		" && rm -rf {params.outdir}" + \
		" &> {log}"

rule qa_quast:
        input:
                assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
        output:
                os.path.join(QA_DIR, "quast", "{sample}", "quast.log")
        params:
                load = loadPreCmd(config["load"]["quast"]),
                outdir = lambda wildcards: os.path.join(QA_DIR, "quast", wildcards.sample),
                cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False)
        log:
                join(QA_DIR, "log", "{sample}.asm_quast_assembly.log")
        threads:
                2
        shell:
                "{params.load}" + TIME_CMD + \
                " {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig 1000 &> {log}"

rule qa_blob_bwa_mem:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
                assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		bam = join(QA_DIR, "blobtools", "bwa", "{sample}", "{sample}.blob_bwa.bam")
	log:
		join(QA_DIR, "log", "{sample}.blob_bwa.log")
	params:
		outdir = lambda wildcards: join(QA_DIR, "blobtools", "bwa", wildcards.sample),
		index = lambda wildcards: join(QA_DIR, "blobtools", "bwa", wildcards.sample, wildcards.sample + ".assembly.fasta"),
		outbam = lambda wildcards: join(QA_DIR, "blobtools", "bwa", wildcards.sample, wildcards.sample + ".blob_bwa.bam"),
		load = loadPreCmd(config["load"]["blob_bwa"])
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + \
		" bwa index -p {params.index} {input.assembly} &&" + \
		" bwa mem -t {threads} {params.index} {input.r1} {input.r2} | samtools view -buSh - > {output.bam} 2> {log}"

rule qa_blob_blast:
	input:
                assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		tsv = join(QA_DIR, "blobtools", "blast", "{sample}", "{sample}.blob_blast.tsv")
	log:
		join(QA_DIR, "log", "{sample}.blob_blast.log")
	params:
		outdir = lambda wildcards: join(QA_DIR, "blobtools", "blast", wildcards.sample),
		load = loadPreCmd(config["load"]["blob_blast"])
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + \
		" blastn -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 " + \
		" -num_threads {threads} -query {input.assembly} -db " + config["blob_blastdb"] + " -out {output.tsv} &> {log}"

rule qa_blobtools:
	input:
		bwa = join(QA_DIR, "blobtools", "bwa", "{sample}", "{sample}.blob_bwa.bam"),
		blast = join(QA_DIR, "blobtools", "blast", "{sample}", "{sample}.blob_blast.tsv"),
		assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		blobtable = join(QA_DIR, "blobtools", "blob", "{sample}", "{sample}.blobDB.table.txt")
	log:
		join(QA_DIR, "log", "{sample}.blobtools.log")
	params:
		prefix = lambda wildcards: join(QA_DIR, "blobtools", "blob", wildcards.sample, wildcards.sample),
		load = loadPreCmd(config["load"]["blobtools"])
	threads:
		1
	shell:
		"{params.load}" + TIME_CMD + \
		" blobtools create -t {input.blast} -b {input.bwa} -i {input.assembly} -o {params.prefix} -x bestsumorder &&" + \
		" blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ -x bestsumorder -r species &&" + \
		" blobtools blobplot -r species -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ &> {log}"























