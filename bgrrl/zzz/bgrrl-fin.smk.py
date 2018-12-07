import sys
import csv
import os
from os.path import join, basename, dirname
import glob

from bgrrl import TIME_CMD
from bgrrl.snakemake_helper import loadPreCmd

DEBUG = config.get("debugmode", False)
# needed for tar-ball generation
EB_ORGANISMS = config.get("enterobase_groups", list())  # ["Salmonella"]
CWD = os.getcwd()

# setup i/o
INPUTDIR = config["out_dir"]
OUTPUTDIR = config["package_dir"]

if DEBUG:
	print("IN_SNAKE:")
	print(config)
	print(EB_ORGANISMS)


### RULES ###

localrules: all

if EB_ORGANISMS:
	rule all:
		input: 
			expand(join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE"), organism=EB_ORGANISMS),
			expand(join(OUTPUTDIR, "{organism}_READ_PKG_DONE"), organism=EB_ORGANISMS)

	rule fin_compile_assembly:
		message:
			"Packaging assemblies by organism..."
		input:
			eb_samples = join(INPUTDIR, "reports", "eb_{organism}_samples.txt")		
		output:
			done = join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE")
		params:
			package_dir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_assemblies"),
			outdir = basename(INPUTDIR),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.package_dir} &&" + \
			" (for s in $(cat {input.eb_samples}); do" + \
			" ln -sf ../../{params.outdir}/assembly/$s/$s.assembly.fasta {params.package_dir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
			" && cd " + CWD + \
			" && md5sum {params.package_dir}.tar.gz > {params.package_dir}.tar.gz.md5" + \
			" && touch {output.done}"

	rule fin_compile_reads:
		message:
			"Packaging reads by organism..."
		input:
			eb_samples = join(INPUTDIR, "reports", "eb_{organism}_samples.txt"),
			samplesheet = config["samplesheet"]
		output:
			done = join(OUTPUTDIR, "{organism}_READ_PKG_DONE")
		params:
			outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_rawreads"),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.outdir} &&" + \
			" (for r in $(grep -F -f {input.eb_samples} {input.samplesheet} | cut -f 3 -d ,); do" + \
			" ln -sf $r {params.outdir}/$(basename $r);" + \
			" r=$(dirname $r)/$(basename $r _R1.fastq.gz)_R2.fastq.gz;" + \
			" ln -sf $r {params.outdir}/$(basename $r); done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
			" && cd " + CWD + \
			" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
			" && touch {output.done}"

elif config["package_mode"] == "processed_reads" or "processed_reads" in config["package_mode"]:
	sentinel = join(OUTPUTDIR, "READ_PKG_DONE")

	rule all:
		input:
			sentinel	
	
	rule fin_package_processed_reads:
		message:
			"Packaging processed_reads..."
		input:
			samplesheet = config["samplesheet"]
		output:
			done = sentinel
		params:
			outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_processed_reads"),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.outdir}" + \
			" && (for r in $(cut -f 11 -d , {input.samplesheet}); do" + \
			" ln -sf $r {params.outdir}/$(basename $r); done)" + \
			" r=$(dirname $r)/$(basename $r _R1.fastq.gz)_R2.fastq.gz;" + \
			" ln -sf $r {params.outdir}/$(basename $r); done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.outdir}.tar.gz $(basename {params.outdir})" + \
			" && cd " + CWD + \
			" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
			" && touch {output.done}"

elif config["package_mode"] == "analysis" or "analysis" in config["package_mode"]:
	sentinel = join(OUTPUTDIR, "ANALYSIS_PKG_DONE")
	
	rule all:
		input:
			sentinel

	rule fin_package_analysis:
		message:
			"Packaging analysis files..."
		input:
			samplesheet = config["samplesheet"]
		output:
			done = sentinel
		params:
			outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_analysis"),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.outdir}" + \
			" && cd {params.outdir}" + \
			" && ln -sf ../../Analysis/reports/16S_report.tsv" + \
			" && ln -sf ../../Analysis/qaa/asm/blobtools" + \
			" && ln -sf ../../Analysis/reports/blobtools_report.tsv" + \
			" && ln -sf ../../Analysis/qa/asm/busco" + \
			" && ln -sf ../../Analysis/reports/multiqc/asm/{params.prefix}_asm_multiqc_report.html" + \
			" && ln -sf ../../Analysis/reports/qc_eval.tsv" + \
			" && ln -sf ../../Analysis/qa/asm/quast" + \
			" && ln -sf ../../Analysis/reports/quast_report.tsv" +\
			# annotation_report.xlsx???
			" && cd .." + \
			" && tar chvzf $(basename {params.outdir}.tar.gz $(basename {params.outdir})" + \
            " && cd " + CWD + \
            " && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
            " && touch {output.done}"


elif config["package_mode"] == "asm" or "asm" in config["package_mode"]:
	sentinel = join(OUTPUTDIR, "ASSEMBLY_PKG_DONE") 

	rule all:
		input: 
			sentinel

	rule fin_package_assembly:
		message:
			"Packaging assemblies..."
		input:
			samples = join(INPUTDIR, "reports", "quast_report.tsv")
		output:
			done = sentinel
		params:
			package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_assemblies"),
			outdir = basename(INPUTDIR),
			prefix = config["misc"]["project"]
		shell:
			"mkdir -p {params.package_dir} &&" + \
			" (for s in $(tail -n +2 {input.samples} | cut -f 1 | grep -v _broken); do" + \
			" ln -sf ../../{params.outdir}/assembly/$s/$s.assembly.fasta {params.package_dir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd $(dirname {params.package_dir})" + \
			" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
			" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
			" && cd -" + \
			" && touch {output.done}"

elif config["package_mode"] == "ann" or "ann" in config["package_mode"]:
	sentinel = join(OUTPUTDIR, "ANNOTATION_PKG_DONE") 

	if config["run_prokka"] and not config["run_ratt"]:
		rule all:
			input: 
				sentinel				

		shell_str = "" + \
			"mkdir -p {params.package_dir}" + \
			" && cwd=$(pwd)" + \
			" && cd {params.package_dir}" + \
			" && {0}" + \
			" && cd .." + \
			" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
			" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
			" && cd $cwd" + \
			" && touch {output.done}"

		if config.get("prokka_package_style", "by_sample") == "by_sample":
			link_command = "ln -sf ../../{0}/annotation/prokka".format(basename(INPUTDIR))
		else:			
			link_command = "find ../../{0}/annotation/prokka -type f -print -exec ln -sf {{}} \;".format(basename(INPUTDIR))

		rule fin_package_annotation:
			message:
				"Packaging prokka annotations (style={})...".format(config.get("prokka_package_style", "by_sample"))
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				done = sentinel
			params:
				package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_prokka_denovo_annotation"),
				outdir = basename(INPUTDIR),
				prefix = config["misc"]["project"],
				lcmd = link_command
			shell:
				"mkdir -p {params.package_dir}" + \
				" && cwd=$(pwd)" + \
				" && cd {params.package_dir}" + \
				" && {params.lcmd}" + \
				" && cd .." + \
				" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
				" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
				" && cd $cwd" + \
				" && touch {output.done}"

	elif config["run_ratt"]:
		sentinel = join(OUTPUTDIR, "ANNOTATION_PKG_DONE")
	
		from math import ceil
		annotation_report = join(INPUTDIR, "reports", "annotation_report.tsv")
		NUM_BATCHES = ceil((len(list(row for row in csv.reader(open(annotation_report), delimiter="\t"))) - 1)/30)

		rule all:
			input:
				sentinel,
				expand(join(OUTPUTDIR, config["misc"]["project"] + "_ratt_annotation_batch.{batch_id}.tar.gz"), batch_id=list(range(1, NUM_BATCHES + 1)))

		rule fin_package_annotation:
			message:
				"Packaging prokka/delta annotations..."
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				done = sentinel
			params:
				package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
				outdir = basename(INPUTDIR),
				prefix = config["misc"]["project"]
			shell:
				"mkdir -p {params.package_dir}/ratt/{{reports,gff}}" + \
				" && ln -sf ../../{params.outdir}/annotation/prokka {params.package_dir}/prokka" + \
				" && ln -sf ../../{params.outdir}/annotation/delta {params.package_dir}/delta" + \
				" && cd {params.package_dir}/ratt/reports" + \
				" && find ../../../../{params.outdir}/annotation/ratt -name '*.ratt_report.tsv' -exec ln -sf {{}} \;" + \
				" && cd - && cd {params.package_dir}/ratt/gff" + \
				" && find ../../../../{params.outdir}/annotation/ratt -name '*.final.gff' -exec ln -sf {{}} \;" + \
				" && cd - && cd $(dirname {params.package_dir})" + \
				" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
				" && echo TARBALL_DONE" + \
				" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
				" && cd -" + \
				" && touch {output.done}"

		rule fin_package_annotation_make_ratt_batches:
			message:
				"Preparing ratt annotation batches..."
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				batches = expand(join(OUTPUTDIR, config["misc"]["project"] + "_ratt_annotation_batch.{batch_id}", "good_to_go"), batch_id=list(range(1, NUM_BATCHES + 1)))
			params:
				package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_ratt_annotation"),
				outdir = basename(INPUTDIR),
				prefix = config["misc"]["project"]
			run:
				import pathlib
				batchid = 0
				for i, row in enumerate(csv.reader(open(input.samples), delimiter="\t"), start=-1):
					if i > -1:
						bid = i % 30
						if bid == 0:
							batchid += 1
							bdir = params.package_dir + "_batch.{}".format(batchid)
							pathlib.Path(bdir).mkdir(parents=True, exist_ok=True)
							open(join(bdir, "good_to_go"), "wt").write("")
						os.symlink(join("..", "..", params.outdir, "annotation", "ratt", row[0]), join(bdir, row[0]))

		rule fin_package_annotation_ratt_tarballs:
			message:
				"Packaging ratt annotations..."
			input:
				goflag = join(OUTPUTDIR, config["misc"]["project"] + "_ratt_annotation_batch.{batch_id}", "good_to_go")
			output:
				tarball = join(OUTPUTDIR, config["misc"]["project"] + "_ratt_annotation_batch.{batch_id}.tar.gz")
			params:
				outdir = OUTPUTDIR, 
				prefix = config["misc"]["project"],
				indir = dirname("{input.goflag}")
			shell:
				"""
				cd {params.outdir} &&
				indir=$(dirname {input.goflag}) &&
				rm -vf {input.goflag} &&
				tar chvzf $(basename {output.tarball}) $(basename $indir) && 
				md5sum $(basename {output.tarball}) > $(basename {output.tarball}).md5
				"""
