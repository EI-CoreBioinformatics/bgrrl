import sys
import subprocess
import argparse
import pathlib
import os
from os.path import join

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("--singularity-container", type=str)
	ap.add_argument("--input", "-i", type=str)
	ap.add_argument("--outdir", "-o", type=str)
	ap.add_argument("--busco-mode", choices=("geno", "genome", "g", "tran", "transcriptome", "t", "prot", "proteome", "p"), default="genome")
	ap.add_argument("--threads", type=int, default=1)
	ap.add_argument("--lineage", "-l", type=str)

	args = ap.parse_args()

	
	outdir = join(os.path.dirname(args.outdir), "run_" + os.path.basename(args.outdir))
	tmpdir = join(outdir, "tmp")
	augustus_dir = join(tmpdir, "config")
	pathlib.Path(tmpdir).mkdir(exist_ok=True, parents=True)
	
	singularity_call = "singularity exec " + args.singularity_container  

	try:
		subprocess.check_call(
			"{} cp -rv /opt/miniconda/config {}/".format(singularity_call, os.path.abspath(tmpdir)),
			shell=True
		)
	except Exception as e:
		print("EXCEPTION1:", e)
		sys.exit(1)

	cmd = "cd {} && " + \
		"export AUGUSTUS_CONFIG_PATH={} && " + \
		"{} run_BUSCO.py -i {} -c {} -m {} --force -t {} -l {} -o {}"
	cmd = cmd.format(
		os.path.abspath(os.path.join(outdir, "..")),
		os.path.abspath(augustus_dir),
		singularity_call,
		os.path.abspath(args.input),
		args.threads,
		{"g": "geno", "t": "tran", "p": "prot"}.get(args.busco_mode[0], "geno"),
		os.path.abspath(tmpdir),
		args.lineage,
		os.path.basename(args.outdir)
	)

	try:
		subprocess.check_call(cmd, shell=True)
	except Exception as e:
		print(e)
		sys.exit(1)
		

	os.rename(outdir, args.outdir)
	open(os.path.join(args.outdir, "short_summary_" + os.path.basename(args.outdir) + ".txt"), "at").close()
	
	os.rename(
		os.path.join(args.outdir, "short_summary_" + os.path.basename(args.outdir) + ".txt"),
		os.path.join(args.outdir, os.path.basename(args.outdir) + "_short_summary.txt")
	)
	

"""


sbatch -p TempProject2 --mem 8GB -c 8 -J BUSCOTEST --wrap "export AUGUSTUS_CONFIG_PATH=$(pwd)/test; singularity exec /tgac/scratch/schudomc/qaa_exec.img_old run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9 -c 8 -f"

==> this will end up with a folder called run_test (-o test) in the current directory
	



 qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
            " {params.load}" + TIME_CMD + \
            " bash -c \"(run_BUSCO.py -i {params.inputpath} -c {threads} -m {params.busco_mode}" + \
            " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample}" + \
            " && touch {params.outdir}/short_summary_{wildcards.sample}.txt)\" &> {log}" + \
            " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
            " && rm -rvf {params.outdir}" + \
            " && (mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt || touch {params.final_outdir}/{wildcards.sample}_short_summary.txt)"




1009  singularity exec /tgac/scratch/schudomc/qaa_exec.img cp -rv /opt/miniconda/config/* .
 1010  singularity exec /tgac/scratch/schudomc/qaa_exec.img ls /opt/miniconda
 1011  singularity exec /tgac/scratch/schudomc/qaa_exec.img ls /opt/miniconda/config
 1012  singularity exec /tgac/scratch/schudomc/qaa_exec.img ls '/opt/miniconda/config/*'
 1013  singularity exec /tgac/scratch/schudomc/qaa_exec.img cp -rv /opt/miniconda/config .
 1014  ls ./config/
 1015  vi config/
 1016  config.ini
 1017  vi config/config.ini
 1018  singularity exec /tgac/scratch/schudomc/qaa_exec.img run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9
 1019  export AUGUSTUS_CONFIG_PATH=$(pwd)/test; singularity exec /tgac/scratch/schudomc/qaa_exec.img run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9
 1020  echo $AUGUSTUS_CONFIG_PATH
 1021  export AUGUSTUS_CONFIG_PATH=$(pwd)/test; singularity exec /tgac/scratch/schudomc/qaa_exec.img_old run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9
 1022  mkdir test
 1023  export AUGUSTUS_CONFIG_PATH=$(pwd)/test; singularity exec /tgac/scratch/schudomc/qaa_exec.img_old run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9
 1024  mv config/* test/
 1025  export AUGUSTUS_CONFIG_PATH=$(pwd)/test; singularity exec /tgac/scratch/schudomc/qaa_exec.img_old run_BUSCO.py -i Analysis/assembly/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001/PRO1946_P1_A1_LTAD-De_CGGTTGGTT-AACCAACCG_L001.assembly.fasta -o test -m geno -l /ei/workarea/group-pb/BUSCO_DATABASES/odb9/bacteria_odb9


    rule qaa_busco_prot:
        input:
            busco_input = getProteins
        output:
            join(qaa_env.busco_prot_dir, "{sample}", "{sample}_short_summary.txt")
        log:
            join(qaa_env.log_dir, "{sample}_busco_prot.log")
        params:
            outdir = lambda wildcards: join(qaa_env.busco_prot_dir, "run_" + wildcards.sample),
            final_outdir = lambda wildcards: join(qaa_env.busco_prot_dir, wildcards.sample),
            tmp = lambda wildcards: join(qaa_env.busco_prot_dir, "tmp", wildcards.sample),
            load = loadPreCmd(config["load"]["busco"]),
            busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
            busco_mode = "prot",
            inputpath = lambda wildcards: getProteins(wildcards) if getProteins(wildcards).startswith("/") else join(qaa_env.cwd, getProteins(wildcards))
        threads:
            8
        shell:
            qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
            " {params.load}" + TIME_CMD + \
            " bash -c \"(run_BUSCO.py -i {params.inputpath} -c {threads} -m {params.busco_mode}" + \
            " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample}" + \
            " && touch {params.outdir}/short_summary_{wildcards.sample}.txt)\" &> {log}" + \
            " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
            " && rm -rvf {params.outdir}" + \
            " && (mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt || touch {params.final_outdir}/{wildcards.sample}_short_summary.txt)"



if config["run_transcriptome_module"]:
    rule qaa_busco_tran:
        input:
            busco_input = getTranscripts
        output:
            join(qaa_env.busco_tran_dir, "{sample}", "{sample}_short_summary.txt")
        log:
            join(qaa_env.log_dir, "{sample}_busco_tran.log")
        params:
            outdir = lambda wildcards: join(qaa_env.busco_tran_dir, "run_" + wildcards.sample),
            final_outdir = lambda wildcards: join(qaa_env.busco_tran_dir, wildcards.sample),
            tmp = lambda wildcards: join(qaa_env.busco_tran_dir, "tmp", wildcards.sample),
            load = loadPreCmd(config["load"]["busco"]),
            busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
            busco_mode = "tran",
            inputpath = lambda wildcards: getTranscripts(wildcards) if getTranscripts(wildcards).startswith("/") else join(qaa_env.cwd, getTranscripts(wildcards))
        threads:
            8
        shell:
            qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
            " {params.load}" + TIME_CMD + \
            " bash -c \"(run_BUSCO.py -i {params.inputpath} -c {threads} -m {params.busco_mode}" + \
            " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample}" + \
            " && touch {params.outdir}/short_summary_{wildcards.sample}.txt)\" &> {log}" + \
            " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
            " && rm -rvf {params.outdir}" + \
            " && (mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt || touch {params.final_outdir}/{wildcards.sample}_short_summary.txt)"


if config["run_genome_module"]:
    if config["run_busco"]:
        rule qaa_busco_geno:
            input:
                busco_input = getAssembly
            output:
                join(qaa_env.busco_geno_dir, "{sample}", "{sample}_short_summary.txt")
            log:
                join(qaa_env.log_dir, "{sample}_busco_geno.log")
            params:
                outdir = lambda wildcards: join(qaa_env.busco_geno_dir, "run_" + wildcards.sample),
                final_outdir = lambda wildcards: join(qaa_env.busco_geno_dir, wildcards.sample),
                tmp = lambda wildcards: join(qaa_env.busco_geno_dir, "tmp", wildcards.sample),
                load = loadPreCmd(config["load"]["busco"]),
                busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
                busco_mode = "geno",
                inputpath = lambda wildcards: getAssembly(wildcards) if getAssembly(wildcards).startswith("/") else join(qaa_env.cwd, getAssembly(wildcards))
            threads:
                8
            shell:
                qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
                " {params.load}" + TIME_CMD + \
                " bash -c \"(run_BUSCO.py -i {params.inputpath} -c {threads} -m {params.busco_mode}" + \
                " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample}" + \
                " && touch {params.outdir}/short_summary_{wildcards.sample}.txt)\" &> {log}" + \
                " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
                " && rm -rvf {params.outdir}" + \
                " && (mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt || touch {params.final_outdir}/{wildcards.sample}_short_summary.txt)"


"""
