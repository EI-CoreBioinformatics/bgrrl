import sys
import os
import glob

TIME_CMD = " /usr/bin/time -v "
BGRRL = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl"

INPUTDIR = os.path.join(os.getcwd(), "Reads")
OUTPUTDIR = os.path.join(os.getcwd(), "Analysis")

print(sys.argv)

SOFTWAREPATH = "/tgac/software/testing"
BUSCO_DATA = os.path.join(BGRRL, "data", "busco", "bacteria_odb9")

BBSUITE_DIR = os.path.join(SOFTWAREPATH, "bbmap", "37.24", "bbmap")

ADAPTERS = os.path.join(BBSUITE_DIR, "resources", "adapters.fa")

# tools
BBDUK = os.path.join(BBSUITE_DIR, "bbduk.sh")
BBNORM = os.path.join(BBSUITE_DIR, "bbnorm.sh")
FASTQC = os.path.join(SOFTWAREPATH, "fastqc", "0.11.5", "x86_64", "bin", "fastqc")

# wrappers - only used until dependencies are stable
BGRRL_WRAPPERS = os.path.join(BGRRL, "scripts", "wrappers")
EST_GSIZE = os.path.join(BGRRL_WRAPPERS, "estimate_genomesize")
SPADES_WRAPPER = os.path.join(BGRRL_WRAPPERS, "spades_wrapper")
UNICYCLER_WRAPPER = os.path.join(BGRRL_WRAPPERS, "unicycler_wrapper")
PROKKA_WRAPPER = os.path.join(BGRRL_WRAPPERS, "prokka_wrapper")
DFAST_WRAPPER = os.path.join(BGRRL_WRAPPERS, "dfast_wrapper")
QUAST_WRAPPER = os.path.join(BGRRL_WRAPPERS, "quast_wrapper")
BUSCO_WRAPPER = os.path.join(BGRRL_WRAPPERS, "busco_wrapper")
REAPR_WRAPPER = os.path.join(BGRRL_WRAPPERS, "reapr_wrapper")

# directories
QC_OUTDIR = os.path.join(OUTPUTDIR, 'qc')
FASTQC_DIR = os.path.join(QC_OUTDIR, 'fastqc')
GSIZE_DIR = os.path.join(QC_OUTDIR, 'gsize')
BBNORM_DIR = os.path.join(QC_OUTDIR, 'bbnorm')


ASSEMBLY_OUTDIR = os.path.join(OUTPUTDIR, 'assembly')
SPADES_OUTDIR = os.path.join(ASSEMBLY_OUTDIR, 'spades')
UNICYCLER_OUTDIR = os.path.join(ASSEMBLY_OUTDIR, 'unicycler')
REAPR_OUTDIR = os.path.join(ASSEMBLY_OUTDIR, 'reapr')

ANNOTATION_OUTDIR = os.path.join(OUTPUTDIR, 'annotation')
PROKKA_OUTDIR = os.path.join(ANNOTATION_OUTDIR, 'prokka')
DFAST_OUTDIR = os.path.join(ANNOTATION_OUTDIR, 'dfast')

QA_OUTDIR = os.path.join(OUTPUTDIR, 'qa')
QUAST_OUTDIR = os.path.join(QA_OUTDIR, 'quast')
BUSCO_OUTDIR = os.path.join(QA_OUTDIR, 'busco')


quarantined = set()
if os.path.exists('quarantine.txt'):
	quarantined = set(line.strip().replace('.fastq.gz', '') for line in open('quarantine.txt'))

INPUTFILES = dict((os.path.basename(_file).replace('_R1.fastq.gz', ''), (_file, _file.replace('_R1.fastq.gz', '_R2.fastq.gz')))
				  for _file in glob.glob(os.path.join(INPUTDIR, '*_R1.fastq.gz')))
INPUTFILES = dict((k, v) for k, v in INPUTFILES.items() if k not in quarantined)
# print(INPUTFILES)

def adequateFileSize(fn):
    return not os.stat(fn).st_size < 1000000



FASTQS = [os.path.basename(_file).replace('.fastq.gz', '')
          for _file in glob.glob(os.path.join(INPUTDIR, '*.fastq.gz'))
	  if adequateFileSize(os.path.join("Analysis", "qc", "tadpole", os.path.basename(_file).replace('_R2', '_R1').replace('_R1.fastq.gz', ''), "tadpole_contigs.fasta"))]
FASTQS = list(filter(lambda s:s.replace('_R1', '') not in quarantined, FASTQS))
print(FASTQS)

TARGETS = list()
TARGETS.extend(map(lambda s:os.path.join(UNICYCLER_OUTDIR, s.replace('_R1', ''), 'assembly.fasta'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(QUAST_OUTDIR, s.replace('_R1', ''), 'quast.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(BUSCO_OUTDIR, s.replace('_R1', ''), 'run_geno', 'short_summary_geno.txt'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))

'''
TARGETS.extend(map(lambda s:os.path.join(PROKKA_OUTDIR, s.replace('_R1', ''), s.replace('_R1', '') + '.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(DFAST_OUTDIR, s.replace('_R1', ''), 'application.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
'''
# stolen from Dan's eipap
def get_sample_files(wc):
	return INPUTFILES[wc.sample]

rule all:
	input: TARGETS

rule unicycler_assembly:
	input:
		r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
	log:
		os.path.join(UNICYCLER_OUTDIR, "log", "{sample}_assembly_unicycler.log")
	params:
		outdir = lambda wildcards: os.path.join(UNICYCLER_OUTDIR, wildcards.sample)
	threads:
		8
	shell:
		UNICYCLER_WRAPPER + \
		" -1 {input.r1} -2 {input.r2} -t {threads} -o {params.outdir} &> {log}"
'''
rule asm_reapr:
	input:
		scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta"),
		r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		os.path.join(REAPR_OUTDIR, "{sample}", "04.break.broken_assembly.fa")
	log:
		os.path.join(REAPR_OUTDIR, "log", "{sample}_reapr.log")
	params:
		outdir = os.path.join(REAPR_OUTDIR, "{sample}")
	threads:
		8
	shell:
		REAPR_WRAPPER + \
		" {input.scaffolds} {params.outdir} {input.r1} {input.r2} {threads}"
'''

rule qa_busco_geno:
	input:
		scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
	output:
		os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno", "short_summary_geno.txt")
	log:
		os.path.join(BUSCO_OUTDIR, "log", "run_geno", "{sample}_busco_geno.log")
	params:
		outdir = os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno"),
                tmp = os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno", "tmp")
	threads:
		8
	shell:
		BUSCO_WRAPPER + \
		" {params.outdir} -i {input.scaffolds} -c {threads} -m geno --force -t {params.tmp} -l " + BUSCO_DATA + " -o geno &> {log}"



rule qa_quast:
	input:
		# scaffolds = os.path.join(SPADES_OUTDIR, "{sample}", "scaffolds.fasta")
		scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
	output:
		os.path.join(QUAST_OUTDIR, "{sample}", "quast.log")
	log:
		os.path.join(QUAST_OUTDIR, "log", "{sample}_qc_quast.log")
	params:
		outdir = lambda wildcards: os.path.join(QUAST_OUTDIR, wildcards.sample)
	threads:
		8
	shell:
		QUAST_WRAPPER + \
		" -o {params.outdir} -t {threads} -L -s {input.scaffolds} &> {log}"

'''
# FOR POST-ASSEMBLY MODULES CHANGE INPUT TO CONSOLIDATED_ASSEMBLY_OUTDIR LATER ON?
rule prokka_annotation:
	input:
		# scaffolds = os.path.join(SPADES_OUTDIR, "{sample}", "scaffolds.fasta")
		scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
	output:
		os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.log")
	log:
		os.path.join(PROKKA_OUTDIR, "log", "{sample}_annotation_prokka.log")
	params:
		outdir = lambda wildcards: os.path.join(PROKKA_OUTDIR, wildcards.sample),
                prefix = lambda wildcards: wildcards.sample
	threads:
		8
	shell:
		PROKKA_WRAPPER + \
		" {params.outdir} {params.prefix} {input.scaffolds} {log} {threads}"

# /usr/bin/time -v python "${DFASTPATH}/dfast" --cpu "${NTHREADS}" --center_name "EI" -g "${INPUT}" -o "${OUTDIR}"  --force &> "${LOG}"
rule dfast_annotation:
	input:
		# scaffolds = os.path.join(SPADES_OUTDIR, "{sample}", "scaffolds.fasta")
		scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
	output:
		os.path.join(DFAST_OUTDIR, "{sample}", "application.log")
	log:
		os.path.join(DFAST_OUTDIR, "log", "{sample}_annotation_dfast.log")
	params:
		outdir = lambda wildcards: os.path.join(DFAST_OUTDIR, wildcards.sample)
	threads:
		8
	shell:
		DFAST_WRAPPER + \
		" {threads} {input.scaffolds} {params.outdir} {log}"
'''
