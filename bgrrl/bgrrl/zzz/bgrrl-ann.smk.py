import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
PROKKA_DIR = join(ANNOTATION_DIR, "prokka")
RATT_DIR = join(ANNOTATION_DIR, "ratt")

PROKKA_WRAPPER = join(config["etc"], "wrappers", "prokka_wrapper")
RATT_WRAPPER = join(config["etc"], "wrappers", "ratt_wrapper")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("ann-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
if config["run_prokka"]:
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".log"), INPUTFILES))
if config["run_ratt"]:
	for d in next(os.walk(config["ratt_reference"]))[1]:     
	    TARGETS.extend(map(lambda s:join(RATT_DIR, s, d, "{}_{}.done".format(s, d)), INPUTFILES))

CWD = os.getcwd()

print("CONFIG")
print(config)

with open("ann-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all

rule all:
	input: TARGETS

if config["run_prokka"]:
	rule ann_prokka:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		output:
			log = join(PROKKA_DIR, "{sample}", "{sample}.log"),
			faa = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.faa"),
			ffn = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.ffn")
		log:
			join(ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
		params:
			outdir = lambda wildcards: join(PROKKA_DIR, wildcards.sample),
			prefix = lambda wildcards: wildcards.sample
		threads:
			8
		shell:
			PROKKA_WRAPPER + \
			" {params.outdir} {params.prefix} {input.contigs} {log} {threads}"

if config["run_ratt"]:
	rule ann_ratt:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta"),
			reference = join(config["ratt_reference"], "{ref_ann}")
		output:
			done = join(RATT_DIR, "{sample}", "{ref_ann}", "{sample}_{ref_ann}.done")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_{ref_ann}.ann_ratt.log")
		params:	
			outdir = lambda wildcards: join(RATT_DIR, wildcards.sample, wildcards.ref_ann),
			load = loadPreCmd(config["load"]["ratt"]),
			# prefix1 = lambda wildcards: wildcards.sample,
			# prefix2 = lambda wildcards: wildcards.ref_ann
		threads:
			1
		shell:
			RATT_WRAPPER + " {params.outdir} " + join(config["cwd"], "{input.contigs}") + \
			" {input.reference} {output.done} &> {log}"
			# " touch " + join(config["cwd"], "{output.done}") + ")" 
			# " && extractfeat -sequence *.final.embl -outseq {params.prefix}.fasta -type CDS" + \
			# " && transeq -sequence {params.prefix}.fasta {params.prefix}.pep.fasta" + \
			# " && mv {params.prefix}.fasta {params.prefix}.ffn" + \
			# " && mv {params.prefix}.pep.fasta {params.prefix}.faa" + \
			#"{params.load}" + \
			#" (cd {params.outdir} && " + TIME_CMD + \
			#" $RATT_HOME/start.ratt.sh {input.reference} " + \
			#join(config["cwd"], "{input.contigs}") + " {params.prefix1}_{params.prefix2} Strain" + \
			#" && mkdir -p nucmer && mv nucmer.* nucmer/" + \
			#" && rm -f DONUCMER.log Sequences query.* Reference.*.fasta" + \
			#" && find . -maxdepth 1 -name '*.tmp2.embl' -exec rm {} \;" + \
			#" && mkdir -p artemis" + \
			#"  && mv *.final.embl artemis/ && cat artemis/*.final.embl > {params.prefix1}_{params.prefix2}.final.embl" + \
			#" && find . -maxdepth 1 -name '*.Report.txt' -size -167c -size +165c -type f -exec rm {} \;" + \
			#" && find . -maxdepth 1 -name '*.Report.gff' -size 0 -exec rm {} \;" + \
			#" && header=0 && (for f in $(find . -maxdepth 1 -name '*.Report.txt'); do" + \
			#" if [[ $header == 0 ]]; then head -n 1 $f | awk -v FS='\t' -v OFS='\t' '{{ print \"Sample\",\"Reference\",$0; }}'; header=1; fi;" + \
			#" tail -n +2 $f | awk -v FS='\t' -v OFS='\t' -v col1='{params.prefix1}' -v col2='{params.prefix2}' '{{ print col1,col2,$0; }}';" + \
			#" rm $f;" + \
			#" done > {params.prefix1}_{params.prefix2}.report.tsv)" + \
			#"" + \
			#" && (for f in $(find . -maxdepth 1 -name '*.Report.gff'); do" + \
			#" echo \"# \"$(basename $f .Report.gff);" + \
			#" cat $f;" + \
			#" rm $f;" + \
			#" done > {params.prefix1}_{params.prefix2}.report.gff)" + \
			#" && find . -maxdepth 1 -name '*.embl' -and -not -name '*.final.embl' -and -not -name '*.NOTTransfered.embl' -exec rm {} \;" + \
			#" && touch " + join(config["cwd"], "{output.done}") + ")" + \
			#" && cd " + CWD + " &> {log}"
