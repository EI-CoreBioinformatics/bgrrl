import sys
import csv
import os
import glob
from os.path import join, basename

from bgrrl.samplesheet import readSamplesheet
from bgrrl import TIME_CMD
from bgrrl.snakemake_helper import loadPreCmd

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
			package_dir = lambda wildcards: join(OUTPUTDIR, config["project_prefix"] + "_" + wildcards.organism + "_assemblies"),
			outdir = basename(INPUTDIR),
			prefix = config["project_prefix"]
		shell:
			"mkdir -p {params.package_dir} &&" + \
			" (for s in $(cat {input.eb_samples}); do" + \
			" ln -s ../../{params.outdir}/assembly/$s/$s.assembly.fasta {params.package_dir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd " + OUTPUTDIR + \
			" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
			" && cd " + CWD + \
			" && md5sum {params.package_dir}.tar.gz > {params.package_dir}.tar.gz.md5" + \
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
			package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_assemblies"),
			outdir = basename(INPUTDIR),
			prefix = config["misc"]["project"]
		shell:
			"mkdir -p {params.package_dir} &&" + \
			" (for s in $(tail -n +2 {input.samples} | cut -f 1 | grep -v _broken); do" + \
			" ln -s ../../{params.outdir}/assembly/$s/$s.assembly.fasta {params.package_dir}/$s.assembly.fasta;" + \
			" done)" + \
			" && cd $(dirname {params.package_dir})" + \
			" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
			" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
			" && cd -" + \
			" && touch {output.done}"

elif config["finalize_mode"] == "ann":
	if config["run_prokka"] and not config["run_ratt"]:

		rule all:
			input: 
				join(OUTPUTDIR, "ANNOTATION_PKG_DONE")

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

			link_command = "ln -s ../../{0}/annotation/prokka".format(basename(INPUTDIR))

		else:
			
			link_command = "find ../../{0}/annotation/prokka -type f -print -exec ln -sf {{}} \;".format(basename(INPUTDIR))

		rule fin_package_annotation:
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				done = join(OUTPUTDIR, "ANNOTATION_PKG_DONE")
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
				

				# "mkdir -p {params.package_dir}" + \
				# "&& ln -s ../../{params.outdir}/annotation/prokka/ {params.package_dir}/prokka" + \
				# " && cd $(dirname {params.package_dir})" + \
				# " && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
				# " && echo TARBALL_DONE" + \
				# " && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
				# " && cd -" + \
				# " && touch {output.done}"

	elif config["run_ratt"]:		
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
				package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
				outdir = basename(INPUTDIR),
				prefix = config["misc"]["project"]
			shell:
				"mkdir -p {params.package_dir}/ratt/{{reports,gff}}" + \
				" && ln -s ../../{params.outdir}/annotation/prokka/ {params.package_dir}/prokka" + \
				" && ln -s ../../{params.outdir}/annotation/delta/ {params.package_dir}/delta" + \
				" && cd {params.package_dir}/ratt/reports" + \
				" && find ../../../../{params.outdir}/annotation/ratt -name '*.ratt_report.tsv' -exec ln -s {{}} \;" + \
				" && cd - && cd {params.package_dir}/ratt/gff" + \
				" && find ../../../../{params.outdir}/annotation/ratt -name '*.final.gff' -exec ln -s {{}} \;" + \
				" && cd - && cd $(dirname {params.package_dir})" + \
				" && tar chvzf $(basename {params.package_dir}).tar.gz $(basename {params.package_dir})" + \
				" && echo TARBALL_DONE" + \
				" && md5sum $(basename {params.package_dir}).tar.gz > $(basename {params.package_dir}).tar.gz.md5" + \
				" && cd -" + \
				" && touch {output.done}"

		rule fin_package_annotation_make_ratt_batches:
			input:
				samples = join(INPUTDIR, "reports", "annotation_report.tsv")
			output:
				batches = expand(join(OUTPUTDIR, config["misc"]["project"] + "_annotation_batch.{batch_id}"), batch_id=list(range(1, NUM_BATCHES + 1)))
			params:
				package_dir = lambda wildcards: join(OUTPUTDIR, config["misc"]["project"] + "_annotation"),
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
						os.symlink(join("..", "..", params.outdir, "annotation", "ratt", row[0]), join(bdir, row[0]))

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
