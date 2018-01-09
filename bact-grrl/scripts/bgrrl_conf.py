import os

TIME_CMD = " /usr/bin/time -v "
BGRRL = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl"

INPUTDIR = os.path.join(os.getcwd(), "Reads")
OUTPUTDIR = os.path.join(os.getcwd(), "Analysis")



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
