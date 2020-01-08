import sys
import csv
import os
from os.path import join

from bgrrl.samplesheet import readSamplesheet, Samplesheet 
from eicore.snakemake_helper import get_cmd_call

TIME_V = config.get("tools", dict()).get("time", "time")

DEBUG = config.get("debugmode", False)

SINGLE_CELL_MODE = config.get("single_cell_mode", False)
SKIP_NORMALIZATION = config.get("no_normalization", False)

# this version experiments without normalization
SKIP_NORMALIZATION = True 

# set up i/o
OUTPUTDIR = config["out_dir"]
QC_OUTDIR = join(OUTPUTDIR, "qc")
QC_LOGDIR = join(QC_OUTDIR, "log")
FASTQC_DIR = join(QC_OUTDIR, "fastqc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
KAT_DIR = join(QC_OUTDIR, "kat")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = join(QC_OUTDIR, "tadpole")

INPUTFILES = Samplesheet(config["samplesheet"])

# generate target list
TARGETS = list()
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R1.bbduk_fastqc.html"), INPUTFILES))
if not SINGLE_CELL_MODE:
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R2.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, s + "_tadpole_contigs.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(KAT_DIR, s, s + ".dist_analysis.json"), INPUTFILES))
TARGETS.extend(map(lambda s:join(BBDUK_DIR, s, s + ".merged.fastq.gz"), INPUTFILES))
TARGETS.extend(map(lambda s:join(BBDUK_DIR, s, s + ".uc_singles.fastq.gz"), INPUTFILES))

TARGETS.append(join(OUTPUTDIR, "reports", "quast_survey_report.tsv"))
TARGETS.append(join(OUTPUTDIR, "reports", "survey_stage_evaluation.tsv"))


print(config)


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
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R1.bbnorm_fastqc.html"), INPUTFILES))
	if not SINGLE_CELL_MODE:
		TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R2.bbnorm_fastqc.html"), INPUTFILES))

if DEBUG:
	with open("inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

# helper
def get_sample_files(wc):
	return INPUTFILES[wc.sample].R1, INPUTFILES[wc.sample].R2

def get_r1(wc):
	return INPUTFILES[wc.sample].R1
def get_r2(wc):
	return INPUTFILES[wc.sample].R2

CMD_CALL = get_cmd_call(config, "bgrrl_container")
QAA_CMD_CALL = get_cmd_call(config, "qaa_container")


### RULES ###

localrules: all, survey_evaluate, collate_quast_reports

rule all:
	input: TARGETS


bbduk_command = TIME_V + \
	" {{params.cmd}}" + \
	" -Xmx30g t={{threads}}" + \
	" in1={{input[0]}}{}" + \
	" out1={{output[0]}}{}" + \
	" ref={{params.adapters}}" + \
	" {{params.bbduk_params}}" + \
	" &> {{log}}"

tadpole_command = TIME_V + \
	" {{params.cmd}}" + \
	" -Xmx30g threads={{threads}}" + \
	" in={{input.r1}}{}" + \
	" out={{output.contigs}} &> {{log}}"

katgcp_command = " ({{params.cmd}} gcp -o {{params.prefix}} -t {{threads}} -v {{input.r1}}{} || touch {{output.katgcp}}) &> {{log}}"


if True:
	rule qc_bbduk:
		message:
			"Preprocessing read data with bbduk..."
		input:
			get_sample_files
		output:
			r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
			r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz"),
			singles = join(BBDUK_DIR, "{sample}", "{sample}_S.bbduk.fastq.gz")
		log:
			join(QC_LOGDIR, "{sample}", "{sample}.qc_bbduk.log")
		resources:
			mem_mb = 32000
		threads:
			8
		params:
			cmd = CMD_CALL + "bbduk.sh",
			adapters = config["resources"]["bb_adapters"],
			bbduk_params = config["params"]["bbduk"]
		shell:
			bbduk_command.format(" in2={input[1]}", " out2={output.r2} outs={output.singles}")


rule qc_fastqc_bbduk:
	message:
		"Generating post-preprocessing report with FastQC..."
	input:
		join(BBDUK_DIR, "{sample}", "{sample}_{mate}.bbduk.fastq.gz")
	output:
		fqc = join(FASTQC_DIR, "bbduk", "{sample}", "{sample}_{mate}.bbduk_fastqc.html")
	params:
		outdir = join(FASTQC_DIR, "bbduk", "{sample}"),
		cmd = CMD_CALL + "fastqc" 
	log:
		join(QC_LOGDIR, "{sample}", "{sample}_{mate}.qc_fastqc_bbduk.log")
	resources:
		mem_mb = 8000
	threads:
		2
	shell:
		"(" + TIME_V + " {params.cmd}" + \
		" --extract --threads={threads} --outdir={params.outdir} {input} " + \
		" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

"""
BB: 13th June 2019
https://twitter.com/BBToolsBio/status/1139272120125386752
In my experience, BBDuk for adapter removal, then Tadpole error correction, then BBMerge (with rsem flag), then Spades with no error correction, enables the fastest, lowest-memory, and most accurate assemblies.  For single cells BBNorm is suggested also, but not for isolates.
"""

rule qc_tadpole_error_correction:
	message:
		"Performing tadpole error correction on reads..."
	input:
		r1 = rules.qc_bbduk.output[0],
		r2 = rules.qc_bbduk.output[1],
		singles = rules.qc_bbduk.output[2]
	output:
		r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.corr.fastq.gz"),
		r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.corr.fastq.gz"),
		discarded_pairs = join(BBDUK_DIR, "{sample}", "{sample}.bbduk.tp_errc.discarded_pairs.fastq.gz"),
		singles = join(BBDUK_DIR, "{sample}", "{sample}_S.bbduk.corr.fastq.gz"),
		discarded_singles = join(BBDUK_DIR, "{sample}", "{sample}.bbduk.tp_errc.discarded_singles.fastq.gz"),
	log:
		join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole_errc.log")
	resources:
		mem_mb = lambda wildcards, attempt: 16000 * attempt
	threads:
		8
	params:
		cmd = CMD_CALL + "tadpole.sh"
	shell:
		TIME_V + \
		" {params.cmd}" + \
		" -Xmx{resources.mem_mb}m threads={threads}" + \
		"  in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outd={output.discarded_pairs} mode=correct &&" + \
		" {params.cmd}" + \
		" -Xmx{resources.mem_mb}m threads={threads}" + \
		" in={input.singles} out={output.singles} outd={output.discarded_singles} mode=correct" + \
		" &> {log}"


rule qc_bbmerge:
	message:
		"Merging reads with bbmerge..."
	input:
		r1 = rules.qc_tadpole_error_correction.output[0],
		r2 = rules.qc_tadpole_error_correction.output[1]
	output:
		merged = join(BBDUK_DIR, "{sample}", "{sample}.merged.fastq.gz"),
		ur1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.corr.unmerged.fastq.gz"),
		ur2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.corr.unmerged.fastq.gz")
	log:
		join(QC_LOGDIR, "{sample}", "{sample}.qc_bbmerge.log")
	resources:
		mem_mb = lambda wildcards, attempt: 8000 * attempt,
		mem_gb = lambda wildcards, attempt: (8000 * attempt) // 1000
	threads:
		8
	params:
		cmd = CMD_CALL + "bbmerge.sh",
		extend2 = 50,
		ext_iterations = 5
	shell:
		TIME_V + " {params.cmd}" + \
		" -Xmx{resources.mem_gb}g threads={threads}" + \
		" in={input.r1} in2={input.r2}" + \
		" out={output.merged} outu1={output.ur1} outu2={output.ur2}" + \
		" rsem=t extend2={params.extend2} iterations={params.ext_iterations} ecct vstrict"

rule qc_concat_singles:
	message:
		"Concatenating single reads for unicycler..."
	input:
		singles = rules.qc_tadpole_error_correction.output.singles,
		merged = rules.qc_bbmerge.output.merged
	output:
		uc_single_reads = join(BBDUK_DIR, "{sample}", "{sample}.uc_singles.fastq.gz")
	resources:
		mem_mb = 2000
	shell:
		"cat {input.singles} {input.merged} > {output.uc_single_reads}"


rule qc_tadpole_survey_assembly:
	message:
		"Generating survey assemblies with tadpole..."
	input:
		ur1 = rules.qc_bbmerge.output.ur1,
		ur2 = rules.qc_bbmerge.output.ur2,
		merged = rules.qc_bbmerge.output.merged,
		singles = rules.qc_tadpole_error_correction.output.singles
	output:
		contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	params:
		cmd = CMD_CALL + "tadpole.sh"
	log:
		join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole_survey_assembly.log")
	resources:
		mem_mb = 32000
	threads:
		8
	shell:
		TIME_V + \
		" {params.cmd}" + \
		" -Xmx30g threads={threads}" + \
		" in={input.ur1} in2={input.ur2} extra={input.merged},{input.singles}" + \
		" out={output.contigs} &> {log}"


rule qaa_quast:
	input:
		assembly = rules.qc_tadpole_survey_assembly.output.contigs
	output:
		join(TADPOLE_DIR, "{sample}", "quast", "transposed_report.tsv")
	params:
		outdir = lambda wildcards: join(TADPOLE_DIR, wildcards.sample, "quast"),
		cmd = QAA_CMD_CALL + "quast.py",
		contiglen = 0 #config["quast_mincontiglen"]
	log:
		join(QC_LOGDIR, "{sample}.quast_tadpole_assembly.log")
	resources:
		mem_mb = 8000
	threads:
		2
	shell:
		" ({params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig {params.contiglen}" + \
		" || touch {params.outdir}/transposed_report.tsv {params.outdir}/report.tsv) 2> {log}" + \
		" && cut -f 1,2 {params.outdir}/report.tsv > {params.outdir}/report.tsv.12" + \
		" && mv {params.outdir}/report.tsv {params.outdir}/report.tsv.full" + \
		" && mv {params.outdir}/report.tsv.12 {params.outdir}/report.tsv"


rule collate_quast_reports:
	input:
		quast_reports = expand(join(TADPOLE_DIR, "{sample}", "quast", "transposed_report.tsv"), sample=INPUTFILES)
	output:
		report = join(OUTPUTDIR, "reports", "quast_survey_report.tsv")
	run:
		from bgrrl.reporters.reporters import collate_quast_reports as report
		report(output.report, *input.quast_reports)


rule survey_evaluate:
	input:
		assembly_stats = expand(join(TADPOLE_DIR, "{sample}", "quast", "transposed_report.tsv"), sample=INPUTFILES),
		read_stats = expand(join(FASTQC_DIR, "bbduk", "{sample}", "{sample}_{mate}.bbduk_fastqc.html"), sample=INPUTFILES, mate=["R1","R2"]),
		seq_stats = expand(join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json"), sample=INPUTFILES)
	output:
		join(OUTPUTDIR, "reports", "survey_stage_evaluation.tsv"),
		join(OUTPUTDIR, "reports", "samplesheets", "samplesheet.survey_pass.yaml")
	run:
		from bgrrl.reporters.survey_stage_evaluation import main as survey_eval_main
		readtype = "bbduk"
		min_tadpole_size = config["minimum_survey_assembly_size"]
		qc_eval_args = list(map(str, ["--readtype", readtype, "--min_tadpole_size", min_tadpole_size, OUTPUTDIR]))
		survey_eval_main(qc_eval_args)


rule qc_katgcp:
	message:
		"Analyzing k-mer distribution, GC content and estimated genome size with kat..."
	input:
		ur1 = rules.qc_bbmerge.output.ur1,                        		
		ur2 = rules.qc_bbmerge.output.ur2,
		merged = rules.qc_bbmerge.output.merged,
		singles = rules.qc_tadpole_error_correction.output.singles
	output:
		katgcp = join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json")
	log:
		join(KAT_DIR, "{sample}", "{sample}.kat")
	params:
		prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample),
		cmd = CMD_CALL + "kat gcp"
	resources:
		mem_mb = 8000
	threads:
		2
	shell:
		"({params.cmd} -o {params.prefix} -t {threads} -v {input.ur1} {input.ur2} {input.merged} {input.singles}" + \
		" || touch {output.katgcp}) &> {log}"



if not SKIP_NORMALIZATION:

	bbnorm_command = TIME_V + " {{params.cmd}}" + \
		" -Xmx30g t={{threads}} in={{input.r1}}{}" + \
		" out={{output.r1}}{}" + \
		" {{params.bbnorm_params}}" + \
		" khist={{output.prehist}} khistout={{output.posthist}} &> {{log}}"

	if True:
		rule qc_bbnorm:
			message:
				"Normalizing read data with bbnorm..."
			input:
				r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
				r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
			output:
				r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
				r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
				prehist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
				posthist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
			params:
				cmd = CMD_CALL + "bbnorm.sh",
				bbnorm_params = config["params"]["bbnorm"]
			log:
				join(QC_LOGDIR, "{sample}", "{sample}.qc_bbnorm.log")
			resources:
				mem_mb = lambda wildcards, attempt: (attempt + 1) * 16000
			threads:
				8
			shell:
				bbnorm_command.format(" in2={input.r2}", " out2={output.r2}")


	rule qc_fastqc_bbnorm:
		message:
			"Generating post-normalization report with FastQC..."
		input:
			join(BBNORM_DIR, "{sample}", "{sample}_{mate}.bbnorm.fastq.gz")
		output:
			fqc = join(FASTQC_DIR, "bbnorm", "{sample}", "{sample}_{mate}.bbnorm_fastqc.html")
		params:
			outdir = join(FASTQC_DIR, "bbnorm", "{sample}"),
			cmd = CMD_CALL + "fastqc"
		log:
			join(QC_LOGDIR, "{sample}", "{sample}_{mate}.qc_fastqc_bbnorm.log")
		resources:
			mem_mb = 8000
		threads:
			2
		shell:
			"({params.cmd}" + \
			" --extract --threads={threads} --outdir={params.outdir} {input}" + \
			" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

