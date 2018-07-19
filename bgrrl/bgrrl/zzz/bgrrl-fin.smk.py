import sys
import csv
import os
import glob
from os.path import join, basename

from bgrrl.samplesheet import readSamplesheet
from bgrrl import TIME_CMD
from eicore.external_process.snakemake_helper import loadPreCmd

INPUTDIR = config["out_dir"]
OUTPUTDIR = config["package_dir"]

# INPUTFILES = dict(readSamplesheet(config["samplesheet"]))

TARGETS = list()
if config["finalize_mode"] == "asm":
	TARGETS.append(join(OUTPUTDIR, "ASSEMBLY_PKG_DONE"))
if config["finalize_mode"] == "ann":
	TARGETS.append(join(OUTPUTDIR, "ANNOTATION_PKG_DONE"))

localrules: all

EB_ORGANISMS = config.get("enterobase_groups", list())  # ["Salmonella"]
# needed for tar-ball generation
CWD = os.getcwd()


print("IN_SNAKE:")
print(config)
print(EB_ORGANISMS)

if EB_ORGANISMS:
	rule all:
		input: 
			expand(join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE"), organism=EB_ORGANISMS),
			expand(join(OUTPUTDIR, "{organism}_READ_PKG_DONE"), organism=EB_ORGANISMS),
		#	expand(join(OUTPUTDIR, "{organism}_RESULTS_PKG_DONE"), organism=EB_ORGANISMS),
		#	expand(join(OUTPUTDIR, config["project_prefix"] + "_{organism}_multiqc_report.html"), organism=EB_ORGANISMS)
		# TARGETS

	rule fin_compile_assembly:
		input:
			eb_samples = join(INPUTDIR, "reports", "eb_{organism}_samples.txt")		
		output:
			done = join(OUTPUTDIR, "{organism}_ASSEMBLY_PKG_DONE")
		params:
			outdir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_assemblies"),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.outdir} &&" + \
			" (for s in $(cat {input.eb_samples}); do" + \
			" ln -s ../../Analysis/assembly/$s/$s.assembly.fasta {params.outdir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
			" && cd " + CWD + \
			" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
			" && touch {output.done}"

	rule fin_compile_reads:
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
			" ln -s $r {params.outdir}/$(basename $r);" + \
			" r=$(dirname $r)/$(basename $r _R1.fastq.gz)_R2.fastq.gz;" + \
			" ln -s $r {params.outdir}/$(basename $r); done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
			" && cd " + CWD + \
			" && md5sum {params.outdir}.tar.gz > {params.outdir}.tar.gz.md5" + \
			" && touch {output.done}"


elif config["finalize_mode"] == "asm":
	rule all:
		input: 
			join(OUTPUTDIR, "ASSEMBLY_PKG_DONE")

	rule fin_package_assembly:
		input:
			samples = join(INPUTDIR, "reports", "quast_report.tsv")
		output:
			done = join(OUTPUTDIR, "ASSEMBLY_PKG_DONE")
		params:
			outdir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_assemblies"),
			prefix = config["misc"]["project"]
		shell:
			"mkdir -p {params.outdir} &&" + \
			" (for s in $(tail -n +2 {input.samples} | cut -f 1 | grep -v _broken); do" + \
			" ln -s ../../Analysis/assembly/$s/$s.assembly.fasta {params.outdir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd $(dirname {params.outdir})" + \
			" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
			" && md5sum $(basename {params.outdir}).tar.gz > $(basename {params.outdir}).tar.gz.md5" + \
			" && cd -" + \
			" && touch {output.done}"

elif config["finalize_mode"] == "ann":
	if not config["run_ratt"]:
		rule all:
			input: 
				join(OUTPUTDIR, "ANNOTATION_PKG_DONE")

		rule fin_package_annotation:
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				done = join(OUTPUTDIR, "ANNOTATION_PKG_DONE")
			params:
				outdir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
				prefix = config["misc"]["project"]
			shell:
				"mkdir -p {params.outdir}" + \
				"&& ln -s ../../Analysis/annotation/prokka/ {params.outdir}/prokka" + \
				" && cd $(dirname {params.outdir})" + \
				" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
				" && echo TARBALL_DONE" + \
				" && md5sum $(basename {params.outdir}).tar.gz > $(basename {params.outdir}).tar.gz.md5" + \
				" && cd -" + \
				" && touch {output.done}"
				
	else:
		from math import ceil
		annotation_report = join(INPUTDIR, "reports", "annotation_report.tsv")
                NUM_BATCHES = ceil((len(list(row for row in csv.reader(open(annotation_report), delimiter="\t"))) - 1)/30)

		rule all:
			input:
				join(OUTPUTDIR, "ANNOTATION_PKG_DONE"),
				expand(join(OUTPUTDIR, config["misc"]["project"] + "_annotation_batch.{batch_id}.tar.gz"), batch_id=list(range(1, NUM_BATCHES + 1)))

		rule fin_package_annotation:
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				done = join(OUTPUTDIR, "ANNOTATION_PKG_DONE"),
			params:
				outdir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
				prefix = config["misc"]["project"]
			shell:
				"mkdir -p {params.outdir}/ratt/{{reports,gff}}" + \
				" && ln -s ../../Analysis/annotation/prokka/ {params.outdir}/prokka" + \
				" && ln -s ../../Analysis/annotation/delta/ {params.outdir}/delta" + \
				" && cd {params.outdir}/ratt/reports" + \
				" && find ../../../../Analysis/annotation/ratt -name '*.ratt_report.tsv' -exec ln -s {{}} \;" + \
				" && cd - && cd {params.outdir}/ratt/gff" + \
				" && find ../../../../Analysis/annotation/ratt -name '*.final.gff' -exec ln -s {{}} \;" + \
				" && cd - && cd $(dirname {params.outdir})" + \
				" && tar chvzf $(basename {params.outdir}).tar.gz $(basename {params.outdir})" + \
				" && echo TARBALL_DONE" + \
				" && md5sum $(basename {params.outdir}).tar.gz > $(basename {params.outdir}).tar.gz.md5" + \
				" && cd -" + \
				" && touch {output.done}"

		rule fin_package_annotation_make_ratt_batches:
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				batches = expand(join(OUTPUTDIR, config["misc"]["project"] + "_annotation_batch.{batch_id}"), batch_id=list(range(1, NUM_BATCHES + 1)))
			params:
				outdir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
				prefix = config["misc"]["project"]
			run:
				import pathlib
				batchid = 0
				for i, row in enumerate(csv.reader(open(input.samples), delimiter="\t"), start=-1):
					if i > -1:
						bid = i % 30
						if bid == 0:
							batchid += 1
							bdir = params.outdir + "_batch.{}".format(batchid)
							pathlib.Path(bdir).mkdir(parents=True, exist_ok=True)
						os.symlink(join("..", "..", "Analysis", "annotation", "ratt", row[0]), join(bdir, row[0]))

		rule fin_package_annotation_ratt_tarballs:
			input:
				indir = join(OUTPUTDIR, config["misc"]["project"] + "_annotation_batch.{batch_id}")
			output:
				tarball = join(OUTPUTDIR, config["misc"]["project"] + "_annotation_batch.{batch_id}.tar.gz")
			params:
                        	outdir = OUTPUTDIR, 
	                        prefix = config["misc"]["project"]
			shell:
				"""
				cd {params.outdir} &&
				tar chvzf $(basename {output.tarball}) $(basename {input.indir}) &&
        	                md5sum $(basename {output.tarball}) > $(basename {output.tarball}).md5
				"""
