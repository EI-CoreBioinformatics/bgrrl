import sys
import csv
import os
import glob
from os.path import join, basename

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

SOFTWAREPATH = "/tgac/software/testing"
INPUTDIR = config["out_dir"]
OUTPUTDIR = config["package_dir"]


INPUTFILES = dict(readSamplesheet(config["samplesheet"]))

TARGETS = list()
TARGETS.append(join(OUTPUTDIR, "ASSEMBLY_PKG_DONE"))

localrules: all, fin_compile_qa

EB_ORGANISMS = ["Salmonella"]

rule all:
	input: 
		expand(join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE"), organism=EB_ORGANISMS),
		expand(join(OUTPUTDIR, "{organism}_READ_PKG_DONE"), organism=EB_ORGANISMS),
		expand(join(OUTPUTDIR, "{organism}_RESULTS_PKG_DONE"), organism=EB_ORGANISMS)
	# TARGETS

rule fin_compile_assembly:
	input:
		eb_samples = join(INPUTDIR, "reports", "{organism}_EB_samples.txt")		
	output:
		done = join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE")
	params:
		outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_assemblies"),
		prefix = config["project_prefix"]
	shell:
		"mkdir -p {params.outdir} &&" + \
		" (for s in $(cat {input.eb_samples});" + \
		" do ln -s ../../Analysis/assembly/$s/$s.assembly.fasta {params.outdir}/$s.assembly.1k.fasta;" + \
		" done) && tar chvzf {params.outdir}.tar.gz {params.outdir} &&" + \
		" md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5 && touch {output.done}"

rule fin_compile_reads:
	input:
		eb_samples = join(INPUTDIR, "reports", "{organism}_EB_samples.txt"),
		samplesheet = config["samplesheet"]
	output:
		done = join(OUTPUTDIR, "{organism}_READ_PKG_DONE")
	params:
		outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_rawreads"),
		prefix = config["project_prefix"]
	shell:
		"mkdir -p {params.outdir} &&" + \
		" (for r in $(grep -F -f {input.eb_samples} {input.samplesheet} | cut -f 3 -d ,); do" + \
		" ln -s $r {params.outdir}/$(basename $r);" + \
		" r=$(dirname $r)/$(basename $r _R1.fastq.gz)_R2.fastq.gz;" + \
		" ln -s $r {params.outdir}/$(basename $r); done)" + \
		" && tar chvzf {params.outdir}.tar.gz {params.outdir}" + \
		" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
		" && touch {output.done}"

rule fin_compile_qa:
	input:
		eb_samples = join(INPUTDIR, "reports", "{organism}_EB_samples.txt")
	output:
		done = join(OUTPUTDIR, "{organism}_RESULTS_PKG_DONE")
	params:
		outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_results"),
		prefix = config["project_prefix"],
		indir = INPUTDIR
	log:
		join(OUTPUTDIR, config["project_prefix"] + "_{organism}.finalize.log")
	shell:
		"mkdir -p {params.outdir}/{{blobtools,busco,fastqc,quast}}" + \
		" && (for s in $(cat {input.eb_samples}); do" + \
		" echo \"$s -> blobtools\";" + \
		" ln -s ../../../Analysis/qa/blobtools/blob/$s/$s.blobDB.table.txt {params.outdir}/blobtools/$s.blobDB.table.txt;" + \
		" echo \"$s -> busco\";" + \
		" ln -s ../../../Analysis/qa/busco/$s/full_table_$s.tsv {params.outdir}/busco/full_table_$s.tsv;" + \
		" ln -s ../../../Analysis/qa/busco/$s/missing_busco_list_$s.tsv {params.outdir}/busco/missing_busco_list_$s.tsv;" + \
		" ln -s ../../../Analysis/qa/busco/$s/short_summary_$s.txt {params.outdir}/busco/short_summary_$s.txt;" + \
		" echo \"$s -> fastqc\";" + \
		" ln -s ../../../Analysis/qc/fastqc/bbnorm/$s/$s\"_R1.bbnorm_fastqc.html\" {params.outdir}/fastqc/$s\"_R1.bbnorm_fastqc.html\";" + \
		" ln -s ../../../Analysis/qc/fastqc/bbnorm/$s/$s\"_R2.bbnorm_fastqc.html\" {params.outdir}/fastqc/$s\"_R2.bbnorm_fastqc.html\";" + \
		" ln -s ../../../Analysis/qc/fastqc/bbnorm/$s/$s\"_R1.bbnorm_fastqc.zip\" {params.outdir}/fastqc/$s\"_R1.bbnorm_fastqc.zip\";" + \
		" ln -s ../../../Analysis/qc/fastqc/bbnorm/$s/$s\"_R2.bbnorm_fastqc.zip\" {params.outdir}/fastqc/$s\"_R2.bbnorm_fastqc.zip\";" + \
		" echo \"$s -> quast\";" + \
		" ln -s ../../../Analysis/qa/quast/$s {params.outdir}/quast/$s;" + \
		" done)" + \
		" && (tar chvzf {params.outdir}.tar.gz {params.outdir}" + \
		" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
		" && touch {output.done}) &> {log}"
