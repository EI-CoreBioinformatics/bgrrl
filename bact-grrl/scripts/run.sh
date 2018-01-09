#!/bin/bash
source perl-5.20.1_gk
source gcc-5.2.0
source python_anaconda-4.4.0_py3_dm
source jre-8u92

BACT_GRRL=/tgac/workarea/group-pb/schudomc_bact/bact-grrl

SLURM_DIR=Analysis/slurm_logs
mkdir -p $SLURM_DIR

if [[ "$1" == "ann" ]]
then
 mod="${BACT_GRRL}/scripts/bgrrl-ann.smk.py"
else
 mod="${BACT_GRRL}/scripts/bgrrl-asm.smk.py"
fi


snakemake --snakefile "${mod}" --cluster-config "${BACT_GRRL}/conf/hpc_config.json" --verbose --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.memory} -J {cluster.J} -x {cluster.xnodes} -e ${SLURM_DIR}/slurm-%N.%j.err -o ${SLURM_DIR}/slurm-%N.%j.out" --cores 32 --latency-wait 120 $2 # --unlock
