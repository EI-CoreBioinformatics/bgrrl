import sys
import os
import glob

from bgrrl_conf import *

INPUTDIR = os.path.join(os.getcwd(), "Reads")
OUTPUTDIR = os.path.join(os.getcwd(), "Analysis")

quarantined = set()
if os.path.exists('quarantine.txt'):
        quarantined = set(line.strip().replace('.fastq.gz', '') for line in open('quarantine.txt'))

INPUTFILES = dict((os.path.basename(_file).replace('_R1.fastq.gz', ''), (_file, _file.replace('_R1.fastq.gz', '_R2.fastq.gz')))
                                  for _file in glob.glob(os.path.join(INPUTDIR, '*_R1.fastq.gz')))
INPUTFILES = dict((k, v) for k, v in INPUTFILES.items() if k not in quarantined)

def adequateFileSize(fn):
    return True
    return os.stat(fn).st_size >= 1000000

def get_sample_files(wc):
    return INPUTFILES[wc.sample]

FASTQS = [os.path.basename(_file).replace('.fastq.gz', '')
          for _file in glob.glob(os.path.join(INPUTDIR, '*.fastq.gz'))
          if adequateFileSize(os.path.join("Analysis", "qc", "tadpole", os.path.basename(_file).replace('_R2', '_R1').replace('_R1.fastq.gz', ''), "tadpole_contigs.fasta"))]
FASTQS = list(filter(lambda s:s.replace('_R2', '_R1').replace('_R1', '') not in quarantined, FASTQS))
print(FASTQS)
print(quarantined)

TARGETS = list()
TARGETS.extend(map(lambda s:os.path.join(PROKKA_OUTDIR, s.replace('_R1', ''), s.replace('_R1', '') + '.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
# TARGETS.extend(map(lambda s:os.path.join(DFAST_OUTDIR, s.replace('_R1', ''), 'application.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(BUSCO_OUTDIR, s.replace('_R1', ''), 'run_tran', 'short_summary_tran.txt'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(BUSCO_OUTDIR, s.replace('_R1', ''), 'run_prot', 'short_summary_prot.txt'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))

rule all:
	input: TARGETS

rule ann_prokka:
    input:
        # needs to point at reapr contigs if used!
        scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
    output:
        log = os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.log"),
        faa = os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.faa"),
        ffn = os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.ffn")
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

'''
rule ann_dfast:
    input:
        # needs to point at reapr contigs if used!
        scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
    output:
        os.path.join(DFAST_OUTDIR, "{sample}", "application.log")        
    log:
        os.path.join(DFAST_OUTDIR, "log", "{sample}_annotation_dfast.log")
    params:
        outdir = lambda wildcards: os.path.join(DFAST_OUTDIR, wildcards.sample)
    shell:
        DFAST_WRAPPER + \
        " {threads} {input.scaffolds} {params.outdir} {log}"
'''


rule qa_busco_tran:
    input:
        # needs to point at ann-consolidate if used!
        transcripts = os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.ffn")
    output:
        os.path.join(BUSCO_OUTDIR, "{sample}", "run_tran", "short_summary_tran.txt")
    log:
        os.path.join(BUSCO_OUTDIR, "log", "{sample}_busco_tran.log")
    params:
        outdir = os.path.join(BUSCO_OUTDIR, "{sample}", "run_tran"),
        tmp = os.path.join(BUSCO_OUTDIR, "{sample}", "run_tran", "tmp")
    threads:
        8
    shell:
        BUSCO_WRAPPER + \
        " {params.outdir} -i {input.transcripts} -c {threads} -m tran" + \
        " --force -t {params.tmp} -l " + BUSCO_DATA + " -o tran &> {log}"


rule qa_busco_prot:
    input:
        # needs to point at ann-consolidate if used!
        proteins = os.path.join(PROKKA_OUTDIR, "{sample}", "{sample}.faa")
    output:
        os.path.join(BUSCO_OUTDIR, "{sample}", "run_prot", "short_summary_prot.txt")
    log:
        os.path.join(BUSCO_OUTDIR, "{sample}_busco_prot.log")
    params:
        outdir = os.path.join(BUSCO_OUTDIR, "{sample}", "run_prot"),
        tmp = os.path.join(BUSCO_OUTDIR, "{sample}", "run_prot", "tmp")
    threads:
        8
    shell:
        BUSCO_WRAPPER + \
        " {params.outdir} -i {input.proteins} -c {threads} -m prot" + \
        " --force -t {params.tmp} -l " + BUSCO_DATA + " -o prot &> {log}"
