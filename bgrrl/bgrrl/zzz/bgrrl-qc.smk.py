import sys
import csv
import os
import glob
from os.path import join, basename

from bgrrl.bgrrl import Sample
from bgrrl import loadPreCmd

TIME_CMD = " /usr/bin/time -v"

BGRRL = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl"
READCOUNT = os.path.join(BGRRL, "scripts", "readcount")

SOFTWAREPATH = "/tgac/software/testing"
BBSUITE_DIR = os.path.join(SOFTWAREPATH, "bbmap", "37.24", "bbmap")
ADAPTERS = os.path.join(BBSUITE_DIR, "resources", "adapters.fa")

# tools
BBDUK = os.path.join(BBSUITE_DIR, "bbduk.sh")
BBNORM = os.path.join(BBSUITE_DIR, "bbnorm.sh")
TADPOLE = os.path.join(BBSUITE_DIR, "tadpole.sh")
FASTQC = os.path.join(SOFTWAREPATH, "fastqc", "0.11.5", "x86_64", "bin", "fastqc")

# wrappers
WRAPPERS = os.path.join(BGRRL, "scripts", "wrappers")
KAT_WRAPPER = os.path.join(WRAPPERS, "kat_wrapper")

# OUTPUTDIR = os.path.join(os.getcwd(), "Analysis")
OUTPUTDIR = config["out_dir"]

QC_OUTDIR = os.path.join(OUTPUTDIR, "qc")
READCOUNT_DIR = os.path.join(QC_OUTDIR, "readcount")
FASTQC_DIR = os.path.join(QC_OUTDIR, "fastqc")
BBDUK_DIR = os.path.join(QC_OUTDIR, "bbduk")
KATHIST_DIR = os.path.join(QC_OUTDIR, "kat", "hist")
KATGCP_DIR = os.path.join(QC_OUTDIR, "kat", "gcp")
BBNORM_DIR = os.path.join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = os.path.join(QC_OUTDIR, "tadpole")

'''
INPUTDIR = os.path.join(os.getcwd(), "Reads")
INPUTFILES = dict((os.path.basename(_file).replace("_R1.fastq.gz", ""),
                   (_file, _file.replace("_R1.fastq.gz", "_R2.fastq.gz")))
                  for _file in glob.glob(os.path.join(INPUTDIR, "*_R1.fastq.gz")))

FASTQS = [os.path.basename(_file).strip(".fastq.gz")
          for _file in glob.glob(os.path.join(INPUTDIR, "*.fastq.gz"))]
'''

# Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
INPUTFILES = dict((row[0], Sample(*row)) for row in csv.reader(open(config["samplesheet"]), delimiter=","))
with open('inputfiles.txt', 'w') as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(FASTQC_DIR, s, s + '_R1.bbduk_fastqc.html'), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, s, s + '_R2.bbduk_fastqc.html'), INPUTFILES))



# TARGETS = list()
# TARGETS.extend(map(lambda s:join(FASTQC_DIR, s.replace('_R2', '_R1').replace('_R1', ''), s + '.bbduk_fastqc.html'), FASTQS))
# TARGETS.extend(map(lambda s:join(KATHIST_DIR, s.replace('_R1', ''), s.replace('_R1', '') + '.kat.hist'), filter(lambda x:x.endswith('_R1'), FASTQS)))
# TARGETS.extend(map(lambda s:join(KATGCP_DIR, s.replace('_R1', ''), s.replace('_R1', '') + '.kat-gcp.mx'), filter(lambda x:x.endswith('_R1'), FASTQS)))
# TARGETS.extend(map(lambda s:join(BBNORM_DIR, s.replace('_R2', '_R1').replace('_R1', ''), s + '.bbnorm.fastq.gz'), FASTQS))
# TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s.replace('_R1', ''), 'tadpole_contigs.fasta'), filter(lambda x:x.endswith('_R1'), FASTQS)))

with open('targets.txt', 'w') as targets_out:
	print(*TARGETS, sep='\n', file=targets_out)

def get_sample_files(wc):
	return INPUTFILES[wc.sample].R1, INPUTFILES[wc.sample].R2


rule all:
	input: TARGETS

rule qc_bbduk:
	input:
		get_sample_files
	output:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_bbduk.log")
	threads:
		8
	params:
		load = loadPreCmd(config["load"]["bbduk"])
	shell:
		"{params.load}" + TIME_CMD + " " + BBDUK + \
		" -Xmx30g t={threads} in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2}" + \
		" ref=" + ADAPTERS + \
		" ktrim=r k=21 mink=11 hdist=2 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"

rule qc_fastqc:
	input:
		os.path.join(BBDUK_DIR, "{sample}", "{fastq}.bbduk.fastq.gz")
	output:
		fqc = os.path.join(FASTQC_DIR, "{sample}", "{fastq}.bbduk_fastqc.html")
	params:
		outdir = os.path.join(FASTQC_DIR, "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"])
	log:
		os.path.join(QC_OUTDIR, "log", "{fastq}.qc_fastqc.log")
	threads:
		2
	shell:
		"{params.load}" + TIME_CMD + " " + FASTQC + \
		" --extract --threads={threads} --outdir={params.outdir} {input} &> {log}"

'''
rule qc_kat_hist:
	input:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	output:
		kat_hist = os.path.join(KATHIST_DIR, "{sample}", "{sample}.kat.hist")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_kathist.log")
	params:
		prefix = lambda wildcards: os.path.join(KATHIST_DIR, wildcards.sample, wildcards.sample + ".kat.hist")
	threads:
		8
	shell:
		KAT_WRAPPER + \
		#" hist -o {params.prefix} -t {threads} -v \<\(zcat {input.r1} {input.r2}\) &> {log}"
		" {input.r1} {input.r2} hist -o {params.prefix} -t {threads} -v &> {log}"

rule qc_kat_gcp:
	input:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	output:
		kat_gcp = os.path.join(KATGCP_DIR, "{sample}", "{sample}.kat-gcp.mx")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_katgcp.log")
	params:
		prefix = lambda wildcards: os.path.join(KATGCP_DIR, wildcards.sample, wildcards.sample + ".kat-gcp")
	threads:
		8
	shell:
		KAT_WRAPPER + \
		" {input.r1} {input.r2} gcp -o {params.prefix} -t {threads} -v &> {log}"

rule qc_bbnorm:
    input:
        r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
        r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
    output:
        r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
        r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
        prehist = os.path.join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
        posthist = os.path.join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
    log:
        os.path.join(QC_OUTDIR, "log", "{sample}.qc_bbnorm.log")
    threads:
        8
    shell:
        TIME_CMD + \
        " " + BBNORM + \
        " -Xmx30g t={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2}" + \
        " target=100 min=2 prefilter" + \
        # ecc" + \
        " khist={output.prehist} khistout={output.posthist} &> {log}"


rule qc_tadpole:
    input:
        r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
        r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
    output:
        contigs = os.path.join(TADPOLE_DIR, "{sample}", "tadpole_contigs.fasta")
    log:
        os.path.join(QC_OUTDIR, "log", "{sample}.qc_tadpole.log")
    threads:
        8
    shell:
        TIME_CMD + \
        " " + TADPOLE + \
        " -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.contigs} &> {log}"
'''
