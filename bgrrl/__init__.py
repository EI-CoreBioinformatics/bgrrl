import pkg_resources

__title__ = "ei:bgrr|"
__author__ = "Christian Schudoma"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2017-2018 Earlham Institute"
__version__ = pkg_resources.require("bgrrl")[0].version

import os
from os.path import join, dirname, basename, exists
import sys
from enum import Enum, unique
from collections import namedtuple, Counter
import yaml
import csv
import shutil
# https://stackoverflow.com/a/600612
import pathlib
from copy import copy

from eicore.external_process.snakemake_helper import *
from eicore import NOW

from bgrrl.enterobase_helpers import validateEnterobaseInput, loadEnterobaseCriteria 
from bgrrl.bin.qc_eval import main as qc_eval_main
from bgrrl.bin.asm_report import main as asm_report_main
from bgrrl.bin.ann_report import main as ann_report_main
from bgrrl.bin.asm_stage_report import main as asm_stage_report_main
from bgrrl.bin.annocmp import main as annocmp_main
from bgrrl.samplesheet import verifySamplesheet, Samplesheet, BaseSample, ASM_Sample

from qaa import QAA_Runner, QAA_ArgumentsAdapter as QAA_Args, DEFAULT_CONFIG_FILE as qaa_config_file, QAA_ID
print("QAA_ID="+QAA_ID)

DEFAULT_QAA_ARGS = QAA_Args(make_input_stream=True,
                            run_blobtools=True,
                            align_reads="bowtie2",
                            qaa_mode="genome",
                            busco_db="bacteria_odb9",
                            run_multiqc=False)

TIME_CMD = " /usr/bin/time -v"

@unique
class PipelineStep(Enum):
    READ_QC = 0
    ASSEMBLY = 1
    ANNOTATION = 2
    DATA_QA = 3
    FINALIZE = 4
    ATTEMPT_FULL = 5


class BGRRLModuleRunner(object):
    def __init__(self, module, args, exe_env, hpc_config, config=dict()):
        print(config)
        self.config = dict(config)
        self.module = module
        self.outdir = args.output_dir
        self.unlock = args.unlock
        self.exe_env = exe_env
        self.config["hpc_config"] = hpc_config

        try:
            sample_sheet = Samplesheet(args.input, sampletype=ASM_Sample)
        except:
            print("This samplesheet {} is not an ASM_Sample. Trying to parse it as BaseSample.")
            sample_sheet = Samplesheet(args.input, sampletype=BaseSample)

        if sample_sheet.verifySampleData(fields=["R1", "R2"]):
            self.config["samplesheet"] = args.input

        # sampletype, readcols = (BaseSample, ["R1", "R2"]) if self.module == "bgrrl-qc" else (ASM_Sample, ["R1", "R2"])
        # if Samplesheet(args.input, sampletype=sampletype).verifySampleData(fields=readcols):
        #     self.config["samplesheet"] = args.input
        self.config["out_dir"] = self.outdir

        for k, v in args._get_kwargs():
            # need to find a solution here how to not override already sanitised input, e.g. enterobase-groups
            # this should do it for the moment, but the whole args-structure needs to be changed
            if k not in self.config:                
                print("WRITING {}: {} TO CONFIG".format(k, v))
                self.config[k] = v

        self.config_file = os.path.join(self.outdir, module + ".conf.yaml")
        with open(self.config_file, "w") as conf_out:
            yaml.dump(self.config, conf_out, default_flow_style=False)

    def run(self):
 
        print("Running " + self.module)
        snake = os.path.join(os.path.dirname(__file__), "zzz", self.module + ".smk.py")
        run_result = run_snakemake(snake, self.outdir, self.config_file, self.exe_env, dryrun=False, unlock=self.unlock)

        return run_result


class BGRRLRunner(object):

    def _handle_output_dir(self, output_dir, overwrite=False):
        outdir_exists = os.path.exists(output_dir)
        if outdir_exists:
            if overwrite:
                print("Output directory already exists and overwrite was requested (-f option).  Deleting directory contents ... ",
                      end="", flush=True)
                print("DEACTIVATED DUE TO TOO MANY ACCIDENTS.")
                # shutil.rmtree(output_dir)
                # os.makedirs(output_dir)
            else:
                print("Output already exists, attempting to resume.", flush=True)
        else:
            pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
        self.logs_dir = join(output_dir, "hpc_logs")
        if not os.path.exists(self.logs_dir) and self.exe_env.use_scheduler:
            print("HPC log dir doesn't exist.  Creating " + self.logs_dir + " now ... ", end="", flush=True)
            pathlib.Path(self.logs_dir).mkdir(parents=True, exist_ok=True)

        print("done.")
        print()

        return outdir_exists

    def __init__(self, args, **kwargs):
        self.args = copy(args) # need to find better solution for this.      

        # Establish a valid cluster configuration... may throw if invalid
        print("Configuring execution environment ... ", end="", flush=True)
        self.exe_env = ExecutionEnvironment(self.args, 
                                            NOW, 
                                            job_suffix=args.input + "_" + args.output_dir, 
                                            log_dir=join(args.output_dir, "hpc_logs"))
        print("done.")
        print(str(self.exe_env))

        # make sure output-directory exists and create hpclog-directory
        outdir_exists = self._handle_output_dir(args.output_dir, overwrite=args.force)
        bginit_bgrrl_cfg = os.path.join(args.output_dir, "config", "bgrrl_config.yaml")
        bginit_hpc_cfg = os.path.join(args.output_dir, "config", "hpc_config.json")
            
        if args.bgrrl_config and os.path.exists(args.bgrrl_config):
            print("Custom BGRRL configuration file specified, overriding defaults")
            self.bgrrl_config_file = args.bgrrl_config
        elif os.path.exists(bginit_bgrrl_cfg):
            print("Found BGRRL configuration file at bginit location ({}), using this.".format(bginit_bgrrl_cfg))
            self.bgrrl_config_file = bginit_bgrrl_cfg
        else:
            raise ValueError("No valid BGRRL configuration specified ({}). Please run bginit or provide a valid configuration file with --bgrrl_config".format(args.bgrrl_config))

        print("Loading BGRRL configuration from {} ...".format(self.bgrrl_config_file), end="", flush=True)
        self.bgrrl_config = yaml.load(open(self.bgrrl_config_file))
        print("done.")
        print()

        if args.hpc_config and os.path.exists(args.hpc_config):
            print("Custom HPC configuration file specified, overriding defaults")
            self.hpc_config = args.hpc_config
            self.args.hpc_config = args.hpc_config # quick fix: this duplicate needs to be dealt with at some point.
        elif os.path.exists(bginit_hpc_cfg):
            print("Found HPC configuration at bginit location ({}), using this.".format(bginit_hpc_cfg))
            self.hpc_config = bginit_hpc_cfg
            self.args.hpc_config = bginit_hpc_cfg # quick fix: this duplicate needs to be dealt with at some point.
        else:       
            raise ValueError("No valid HPC configuration specified ({}). Please run bginit or provide a valid configuration file with --hpc_config".format(args.hpc_config))









        print()

        # Set run mode
        print(args.module.upper())
        self.run_mode = PipelineStep[args.module.upper()]


    def _run_qc(self, args): 
        readtype = "bbduk" if self.args.no_normalization else "bbnorm"
        qaa_args = copy(args)

        if self.args.report_only:
            run_result = qc_eval_main(["--readtype", readtype, self.args.input, self.args.output_dir])
        else:
            run_result = BGRRLModuleRunner("bgrrl-qc", self.args, self.exe_env, self.hpc_config, config=self.bgrrl_config).run()
            if run_result:
                qaa_args.update(**vars(self.args),
                                survey_assembly=True, 
                                runmode="survey", 
                                run_blobtools=False,
                                run_busco=False, 
                                normalized=not self.args.no_normalization)
                qaa_run = QAA_Runner(qaa_args).run()					
                if qaa_run:
                    run_result = qc_eval_main(["--readtype", readtype, self.args.input, self.args.output_dir])
                    if run_result:
                        self.args.input = join(self.args.output_dir, "reports", "samplesheets", "samplesheet.qc_pass.tsv")
                        qaa_args.update(**vars(self.args),
                                        run_blobtools=True, 
                                        run_busco=True, 
                                        run_multiqc=True, 
                                        multiqc_dir=join(self.args.output_dir, "reports", "multiqc", "qc"))
                        qaa_run = QAA_Runner(qaa_args).run()

        return run_result


    def _run_asm(self, args):
        qaa_args = copy(args)

        self.bgrrl_config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
        self.bgrrl_config["cwd"] = os.getcwd()
        self.bgrrl_config["reapr_correction"] = False
    
        if self.args.contig_minlen: # this is really, really, really bad! 
            self.bgrrl_config["use_asm_lengthfilter"] = True
            self.bgrrl_config["asm_lengthfilter_contig_minlen"] = self.args.contig_minlen
        else:
            self.bgrrl_config["use_asm_lengthfilter"] = False
            self.bgrrl_config["asm_lengthfilter_contig_minlen"] = 0

        eb_criteria = self.bgrrl_config.get("enterobase_criteria", "")

        if self.args.report_only:
            run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
            if self.args.enterobase_groups: # needs validation?
                run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
        else:
            run_result = BGRRLModuleRunner("bgrrl-asm", self.args, self.exe_env, self.hpc_config, config=self.bgrrl_config).run() 
            if run_result:
                run_result = asm_stage_report_main([self.args.output_dir, join(self.args.output_dir, "reports")])
                if run_result:
                        qaa_args.update(**vars(self.args),  # <- s. above, this is kinda stupid, redesign ASAP!
                                        survey_assembly=False,
                                        runmode="asm",
                                        run_multiqc=True,
                                        multiqc_dir=join(self.args.output_dir, "reports", "multiqc", "asm"))
                        qaa_run = QAA_Runner(qaa_args).run()
                        if qaa_run:
                            self.args.finalize_mode = "asm"
                            if self.args.enterobase_groups:
                                run_result = asm_report_main([self.args.output_dir, self.args.enterobase_groups, eb_criteria])
                            if run_result:
                                run_result = self._run_fin() 

        return run_result

    def _run_ann(self, args):
        qaa_args = copy(args)
        self.bgrrl_config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
        self.bgrrl_config["cwd"] = os.getcwd()
    
        self.bgrrl_config["run_ratt"] = self.args.annotation in ("both", "ratt")
        self.bgrrl_config["run_prokka"] = self.args.annotation in ("both", "prokka")
        self.bgrrl_config["ratt_reference"] = self.args.ratt_reference_dir

        assert not self.bgrrl_config["run_ratt"] or os.path.exists(self.bgrrl_config["ratt_reference"]), "Missing reference data for ratt. Please make sure to use the --ratt-reference-dir parameter."
        run_result = False
        print(self.args)

        if self.args.report_only:
            run_result = False
            if self.args.annotation in ("ratt", "both"):
                run_result = ann_report_main(["--ref-dir", self.args.ratt_reference_dir, join(self.args.output_dir, "annotation", "ratt")])
                annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
        else:
            if self.args.annotation in ("prokka", "both"):
                print("WARNING: Prokka annotation selected. If your jobs fail, you might have to update tbl2asn and/or exclude nodes (hmmscan/GNU parallel fails).")
                run_result = BGRRLModuleRunner("bgrrl-ann", self.args, self.exe_env, self.hpc_config, config=self.bgrrl_config).run()
                if not run_result:
                    print("ANNOTATION RUN FAILED?")
                if run_result:
                    qaa_args.update(**vars(self.args),
                                    survey_assembly=False, 
                                    qaa_mode="transcriptome,proteome", 
                                    runmode="ann")
                    qaa_run = QAA_Runner(qaa_args).run()
    
                    if qaa_run and self.args.annotation in ("ratt", "both"):
                        ann_report_main(["--ref-dir", self.args.ratt_reference_dir, join(self.args.output_dir, "annotation", "ratt")])
                        annocmp_main([join(self.args.output_dir, "annotation", "prokka"), join(self.args.output_dir, "annotation", "ratt"), join(self.args.output_dir, "reports")])
                    else:
                        open(join(self.args.output_dir, "reports", "annotation_report.tsv"), "at")
                    if qaa_run:
                        self.args.finalize_mode = "ann"
                        run_result = self._run_fin() 
        return run_result

    def _run_fin(self):
        self.bgrrl_config["package_dir"] = os.path.join(os.path.dirname(self.args.output_dir), "Data_Package")
        if "enterobase_groups" in self.args: 
            try:
                eb_criteria = loadEnterobaseCriteria(self.bgrrl_config["enterobase_criteria"])
            except:
                print("Enterobase-groups selected but missing enterobase_criteria entry in bgrrl_config, please add path to enterobase criteria.", file=sys.stderr)
                sys.exit(1)
            self.bgrrl_config["enterobase_groups"] = validateEnterobaseInput(self.args.enterobase_groups, eb_criteria)
        else:
            self.bgrrl_config["enterobase_groups"] = list()
        if "project_prefix" in self.args and self.args.project_prefix:
            self.bgrrl_config["misc"]["project"] = self.args.project_prefix
        print("_FIN_CONFIG")
        print(self.bgrrl_config)

        run_result = BGRRLModuleRunner("bgrrl-fin", self.args, self.exe_env, self.hpc_config, config=self.bgrrl_config).run()
        return run_result

    def _run_all(self, args):
        qaa_args = copy(args)
        print("Wrong runmode: (ATTEMPT_FULL is not implemented yet)", self.run_mode)
        exit(1)

    def run(self):
        qaa_args = DEFAULT_QAA_ARGS
        qaa_args.update(quast_mincontiglen=1000, 
                        project_prefix=self.args.project_prefix, 
                        config=self.bgrrl_config_file, 
                        hpc_onfig=self.hpc_config)
    
        if self.run_mode == PipelineStep.READ_QC:
            run_result = self._run_qc(qaa_args)
        elif self.run_mode == PipelineStep.ASSEMBLY:
            run_result = self._run_asm(qaa_args) 
        elif self.run_mode == PipelineStep.ANNOTATION:
            run_result = self._run_ann(qaa_args)
        elif self.run_mode == PipelineStep.FINALIZE:
            run_result = self._run_fin() 
        else:
            run_result = self._run_all(qaa_args)        
    
        print()
        if run_result:
            print("BGRRL completed successfully.")
        else:
            print("BGRRL failed.  Please consult logs to debug run.")
            exit(1)

        return True

