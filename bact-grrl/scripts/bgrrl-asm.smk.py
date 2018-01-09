import sys
import os
import glob

from bgrrl_conf import *


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
FASTQS = list(filter(lambda s:s.replace('_R1', '') not in quarantined, FASTQS))
print(FASTQS)


TARGETS = list()
TARGETS.extend(map(lambda s:os.path.join(UNICYCLER_OUTDIR, s.replace('_R1', ''), 'assembly.fasta'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(QUAST_OUTDIR, s.replace('_R1', ''), 'quast.log'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))
TARGETS.extend(map(lambda s:os.path.join(BUSCO_OUTDIR, s.replace('_R1', ''), 'run_geno', 'short_summary_geno.txt'), (fastq for fastq in FASTQS if fastq.endswith('_R1'))))


rule all:
    input: TARGETS

rule asm_unicycler:
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
# untested!
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
        # needs to point at reapr contigs if used!
        scaffolds = os.path.join(UNICYCLER_OUTDIR, "{sample}", "assembly.fasta")
    output:
        os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno", "short_summary_geno.txt")
    log:
        os.path.join(BUSCO_OUTDIR, "{sample}_busco_geno.log")
    params:
        outdir = os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno"),
        tmp = os.path.join(BUSCO_OUTDIR, "{sample}", "run_geno", "tmp")
    threads:
        8
    shell:
        BUSCO_WRAPPER + \
        " {params.outdir} -i {input.scaffolds} -c {threads} -m geno" + \
        " --force -t {params.tmp} -l " + BUSCO_DATA + " -o geno &> {log}"

rule qa_quast:
    input:
        # needs to point at reapr contigs if used!
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
