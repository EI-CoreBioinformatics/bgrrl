#!/usr/bin/env python3

import os
from os.path import join
import csv
import sys
import shutil
import yaml
# https://stackoverflow.com/a/600612
import pathlib 
from argparse import ArgumentParser
from textwrap import dedent
import copy

from snakemake.utils import min_version
min_version("4.0")

from eicore import NOW
from eicore.external_process.snakemake_helper import make_exeenv_arg_group, ExecutionEnvironment 
from . import DEFAULT_HPC_CONFIG_FILE, DEFAULT_BGRRL_CONFIG_FILE, PipelineStep, __version__, BGRRLModuleRunner, BGRRLRunner

"""
from bgrrl.bgrrl import validateEnterobaseInput, ENTERO_CRITERIA
from bgrrl.bin.qc_eval import main as qc_eval_main
from bgrrl.bin.asm_report import main as asm_report_main
from bgrrl.bin.ann_report import main as ann_report_main
from bgrrl.bin.asm_stage_report import main as asm_stage_report_main
from bgrrl.bin.annocmp import main as annocmp_main

from qaa import QAA_Runner, QAA_ArgumentsAdapter as QAA_Args, DEFAULT_CONFIG_FILE as qaa_config_file, QAA_ID
print("QAA_ID="+QAA_ID)

DEFAULT_QAA_ARGS = QAA_Args(make_input_stream=True, 
                            no_blobtools=False, 
                            blobtools_no_bwa=False, 
                            qaa_mode="genome", 
                            busco_db="bacteria_odb", 
                            no_multiqc=True)
"""


VALID_ASSEMBLERS = ["unicycler", "velvet"]
VALID_ANNOTATERS = ["prokka", "ratt", "both"]

	
def main():
    print("Starting EI BGRRL V " + __version__)
    print()

    # bgrrl_config = yaml.load(open(DEFAULT_BGRRL_CONFIG_FILE))
    bgrrl_config_file = DEFAULT_BGRRL_CONFIG_FILE

    parser = ArgumentParser("The Earlham Institute Bacterial Genome Reconstruction & Recognition Pipeline (BGRR|)",
                            description="""This program controls the various Snakemake pipelines making up the EI-BGRR| pipeline.""")

    parser.add_argument("input", 
                        help="""The (EI-PAP) samplesheet to process. This is a comma-separated file, 
                                containing location and meta-information for each sample to be processed.""")
    parser.add_argument("-o", "--output_dir", default="BGRRL_<timestamp>",
	                help="If specified BGRRL will output data to this directory.")

    parser.add_argument("--project-prefix", 
                        default=os.path.basename(os.getcwd()), 
                        help="Reports and resultfiles/-folders will be prefixed by this.")

    parser.add_argument("-f", "--force", 
                        action="store_true", 
                        help="Force overwriting existing output directory, causes pipeline to be restarted. (disabled)")

    parser.add_argument("--bgrrl_config", 
                        help="""Configuration file for BGRRL. This file specifies details for accessing services and commands 
                                to be executed prior to running each pipeline tool.  
                                The default config file is located at {}""".format(DEFAULT_BGRRL_CONFIG_FILE))

    parser.add_argument("-m", "--module", 
                        choices=[ps.name.lower() for ps in PipelineStep], 
                        default="read_qc",
                        help=dedent("""This option controls which part of the pipeline to run. Due to the task at hand 
                                       being highly dependent on input data quality,
                                       it is recommended to run the 4 major steps of EI-BGRRL manually, carrying out manual/visual 
                                       data checks, and additional data manipulation steps in between steps.
                                       \"read_qc\" - This step preprocesses the read sets and assesses the assemblability 
                                       of individual samples. It will produce an input-samplesheet for the assembly-step.
                                       \"assembly\" - This step performs genome assembly and post-assembly quality control for 
                                       each sample.
                                       \"annotation\" - This step performs genome annotation and post-annotation quality control 
                                       for each assembly.
                                       \"report_and_package\" - This step performs report generation and data packaging.
                                       \"attempt_full\" - If you know the data and what you're doing, this option attempts to run the 
                                       individual pipeline steps consecutively without breaks (not implemented yet)"""))

    parser.add_argument("--report-only", 
                        action="store_true", 
                        help="Only runs reporting modules, no snakemake pipelines [False]")

    parser.add_argument("--no-normalization", 
                        action="store_true", 
                        help="Use non-normalized reads in asm module (kills fallback-mode!) [False]")

    parser.add_argument("--finalize-mode", 
                        choices=["ann", "asm"],
                        help="""Packaging runs automatically after the assembly/annotation steps. 
                                In case packages need to be regenerated, the finalize step can be ran for each module separately.""")

    parser.add_argument("--fin-report-only", 
                        action="store_true", 
                        help="""If the finalize-module is called with this parameter, then it will only generate reports 
                               and omit generating data packages. [False]""")

    parser.add_argument("--enterobase-groups", 
                        type=str, 
                        default="", 
                        help="""Comma-separated list of Enterobase microbial organisms. 
                                The set of assemblies is tested against organism-specific criteria and assemblies are 
                                packaged according to their species. [NEEDS REWORDING!]. 
                                By default, the enterobase mode is disabled.""")

    parser.add_argument("--assembler", 
                        type=str, 
                        choices=VALID_ASSEMBLERS,
                        default="unicycler", 
                        help="""Assembly software to use for genome assembly. [unicycler]""")

    parser.add_argument("--contig-minlen", 
                        type=int, 
                        default=0, 
                        help="Minimum length [bp] of contigs retained in filtering step [0].")

    parser.add_argument("--annotation", 
                        type=str, 
                        choices=VALID_ANNOTATERS, 
                        help="""Annotation software to use for genome annotation. [prokka]""")

    parser.add_argument("--ratt-reference-dir", 
                        type=str, 
                        help="Path to reference data for ratt", 
                        default="")

    make_exeenv_arg_group(parser)	# Add in cluster and DRMAA options
    args = parser.parse_args()

    bgrrl_runner = BGRRLRunner(args)
    run_result = bgrrl_runner.run()

"""
    # Set run mode
    run_mode = PipelineStep[args.module.upper()]

    # Establish a valid cluster configuration... may throw if invalid
    print("Configuring execution environment ... ", end="", flush=True)
    exe_env = ExecutionEnvironment(args, NOW, job_suffix=args.input + "_" + args.output_dir, log_dir=join(args.output_dir, "hpc_logs"))
    print("done.")
    print(str(exe_env))
    print()

    if args.bgrrl_config:
        print("Custom BGRRL configuration file specifed, overriding defaults ... ", end="", flush=True)
        bgrrl_config = yaml.load(open(args.bgrrl_config))
        bgrrl_config_file = args.bgrrl_config
        print("done.")
        print()
    else:
        print("Loading default BGRRL configuration ...", end="", flush=True)
        bgrrl_config = yaml.load(open(bgrrl_config_file))
        print("done.")
        print()

    if os.path.exists(args.output_dir):
        if args.force:
            print("Output directory already exists and overwrite was requested (-f option).  Deleting directory contents ... ",
                  end="",
                  flush=True)
            print("DEACTIVATED DUE TO TOO MANY ACCIDENTS.")
            # shutil.rmtree(args.output_dir)
            # os.makedirs(args.output_dir)
            print("done.")
        else:
            print("Output already exists, attempting to resume.", flush=True)
    else:
        print("Output directory doesn't exist creating ... ", end="", flush=True)
        os.makedirs(args.output_dir)
        print("done.")

    logs_dir = join(args.output_dir, "hpc_logs")
    if not os.path.exists(logs_dir) and exe_env.use_scheduler:
        print("HPC log dir doesn't exist.  Creating " + logs_dir + " now ... ", end="", flush=True)
        os.makedirs(logs_dir)
        print("done.")

    print()
    print(args.module.upper())

    #qaa_args = {
    #    # "config": qaa_config_file,
    #    "make_input_stream": True,
    #    "no_blobtools": False,
    #    "blobtools_no_bwa": False,
    #    "quast_mincontiglen": 1000,
    #    "busco_db": "bacteria_odb9",
    #    "qaa_mode": "genome",
    #    "no_multiqc": True,
    #    "project_prefix": args.project_prefix,
    #    "config": bgrrl_config_file,
    #    "hpc_config": args.hpc_config}
    qaa_args = DEFAULT_QAA_ARGS
    qaa_args.update(quast_mincontiglen=1000, 
                    project_prefix=args.project_prefix, 
                    config=bgrrl_config_file, 
                    hpc_onfig=args.hpc_config)

    if run_mode == PipelineStep.READ_QC:
        readtype = "bbduk" if args.no_normalization else "bbnorm"

        if args.report_only:
            run_result = qc_eval_main(["--readtype", readtype, args.input, args.output_dir])
        else:
            run_result = BGRRLModuleRunner("bgrrl-qc", args, exe_env, config=bgrrl_config).run()
            # run_result = run_qc(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
            if run_result:
                qaa_args.update(**vars(args),
                                survey_assembly=True, 
                                runmode="survey", 
                                no_blobtools=True, 
                                no_busco=True, 
                                normalized=not args.no_normalization)
                # qaa_args["survey_assembly"] = True
                # qaa_args["runmode"] = "survey"
                # qaa_args["no_blobtools"] = True
                # qaa_args["no_busco"] = True
                # qaa_args["normalized"] = not args.no_normalization
                qaa_run = QAA_Runner(qaa_args).run()					
                if qaa_run:
                    run_result = qc_eval_main(["--readtype", readtype, args.input, args.output_dir])
                    if run_result:
                        args.input = join(args.output_dir, "reports", "samplesheets", "samplesheet.qc_pass.tsv")
                        qaa_args.update(**vars(args),
                                        no_blobtools=False, 
                                        no_busco=False, 
                                        no_multiqc=False, 
                                        multiqc_dir=join(args.output_dir, "reports", "multiqc", "qc"))
                        #qaa_args["no_blobtools"] = False
                        #qaa_args["no_busco"] = False
                        #qaa_args["no_multiqc"] = False
                        #qaa_args["multiqc_dir"] = join(args.output_dir, "reports", "multiqc", "qc") 	
                        qaa_run = QAA_Runner(qaa_args).run()

    elif run_mode == PipelineStep.ASSEMBLY:
        bgrrl_config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
        bgrrl_config["cwd"] = os.getcwd()
        # config["assembler"] = args.assembler
        bgrrl_config["reapr_correction"] = False
        # config["no_normalization"] = args.no_normalization
    
        if args.contig_minlen:
            bgrrl_config["use_asm_lengthfilter"] = True
            bgrrl_config["asm_lengthfilter_contig_minlen"] = args.contig_minlen
        else:
            bgrrl_config["use_asm_lengthfilter"] = False
            bgrrl_config["asm_lengthfilter_contig_minlen"] = 0

        if args.report_only:
            run_result = asm_stage_report_main([args.output_dir, join(args.output_dir, "reports")])
            if args.enterobase_groups: # needs validation?
                run_result = asm_report_main([args.output_dir, args.enterobase_groups])
        else:
            # run_result = run_asm(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
            run_result = BGRRLModuleRunner("bgrrl-asm", args, exe_env, config=bgrrl_config).run() 
            if run_result:
                run_result = asm_stage_report_main([args.output_dir, join(args.output_dir, "reports")])
                if run_result:
                        qaa_args.update(**vars(args),
                                        survey_assembly=False,
                                        runmode="asm",
                                        no_multiqc=False,
                                        multiqc_dir=join(args.output_dir, "reports", "multiqc", "asm"))
                        # qaa_args["survey_assembly"] = False
                        # qaa_args["runmode"] = "asm"
                        # qaa_args["no_multiqc"] = False
                        # qaa_args["multiqc_dir"] = join(args.output_dir, "reports", "multiqc", "asm")
                        qaa_run = QAA_Runner(qaa_args).run()
                        if qaa_run:
                            args.finalize_mode = "asm"
                            if args.enterobase_groups:
                            	run_result = asm_report_main([args.output_dir, args.enterobase_groups])
                            if run_result:
                                run_result = run_fin(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)

    elif run_mode == PipelineStep.ANNOTATION:
        bgrrl_config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
        bgrrl_config["cwd"] = os.getcwd()
    
        bgrrl_config["run_ratt"] = args.annotation in ("both", "ratt")
        bgrrl_config["run_prokka"] = args.annotation in ("both", "prokka")
        bgrrl_config["ratt_reference"] = args.ratt_reference_dir

        assert not bgrrl_config["run_ratt"] or os.path.exists(bgrrl_config["ratt_reference"]), "Missing reference data for ratt. Please make sure to use the --ratt-reference-dir parameter."

        if args.report_only:
            run_result = False
            if args.annotation in ("ratt", "both"):
                run_result = ann_report_main(["--ref-dir", args.ratt_reference_dir, join(args.output_dir, "annotation", "ratt")])
                annocmp_main([join(args.output_dir, "annotation", "prokka"), join(args.output_dir, "annotation", "ratt"), join(args.output_dir, "reports")])
        else:
            if args.annotation in ("prokka", "both"):
                print("WARNING: Prokka annotation selected. If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails).")
                # run_result = run_ann(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
                run_result = BGRRLModuleRunner("bgrrl-ann", args, exe_env, config=bgrrl_config).run()
                if not run_result:
                    print("ANNOTATION RUN FAILED?")
                if run_result:
                    qaa_args.update(**vars(args),
                                    survey_assembly=False, 
                                    qaa_mode="transcriptome,proteome", 
                                    runmode="ann")
                    # qaa_args["survey_assembly"] = False
                    # qaa_args["qaa_mode"] = "transcriptome,proteome"
                    # qaa_args["runmode"] = "ann"
                    qaa_run = QAA_Runner(qaa_args).run()
    
                    if qaa_run and args.annotation in ("ratt", "both"):
                        ann_report_main(["--ref-dir", args.ratt_reference_dir, join(args.output_dir, "annotation", "ratt")])
                        annocmp_main([join(args.output_dir, "annotation", "prokka"), join(args.output_dir, "annotation", "ratt"), join(args.output_dir, "reports")])
                    if qaa_run:
                        args.finalize_mode = "ann"
                        run_result = run_fin(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)

    elif run_mode == PipelineStep.FINALIZE:
        bgrrl_config["package_dir"] = os.path.join(os.path.dirname(args.output_dir), "Data_Package")
        # config["project_prefix"] = args.project_prefix
        bgrrl_config["enterobase_groups"] = validateEnterobaseInput(args.enterobase_groups, ENTERO_CRITERIA) if "enterobase_groups" in args else list()

        # if "finalize_mode" in args:
        #    config["finalize_mode"] = args.finalize_mode

        #if not args.fin_report_only:
        #run_result = run_fin(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
        # else:
        # run_result = True
        # args.finalize_mode 
        # run_result = run_fin(args.input, args.output_dir, args, exe_env, bgrrl_config=bgrrl_config)
        run_result = BGRRLModuleRunner("bgrrl-fin", args, exe_env, config=bgrrl_config).run()
    else:
        print("Wrong runmode: (ATTEMPT_FULL is not implemented yet)", run_mode)
        exit(1)

    print()
    if run_result:
        print("BGRRL completed successfully.")
    else:
        print("BGRRL failed.  Please consult logs to debug run.")
        exit(1)
"""


if __name__ == "__main__":
    main()
