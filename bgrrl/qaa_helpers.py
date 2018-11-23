from qaa import QAA_ArgumentsAdapter as QAA_Args


DEFAULT_QAA_ARGS = QAA_Args(make_input_stream=True,
                            run_blobtools=True,
                            align_reads="bowtie2",
                            qaa_mode="genome",
                            busco_db="bacteria_odb9",
                            run_multiqc=False)

class QAA_ArgumentManager(object):
    @staticmethod
    def get_qaa_args(args, bgrrl_config_file, hpc_config_file, stage="init"):

        qaa_args = DEFAULT_QAA_ARGS
        qaa_args.update(quast_mincontiglen=1000,
                        project_prefix=args.project_prefix,
                        config=bgrrl_config_file,
                        hpc_config=hpc_config_file)

        if stage == "init":
            pass
        if stage == "qc_survey":
            qaa_args.update(**vars(args),
                            survey_assembly=True,
                            runmode="survey",
                            run_blobtools=False,
                            align_reads=False,
                            run_busco=False,
                            normalized=not args.no_normalization)
        elif stage == "qc_report":
            qaa_args.update(**vars(args),
                            run_blobtools=True,
                            align_reads=aligner,
                            run_busco=True,
                            run_multiqc=True,
                            multiqc_dir=join(args.output_dir, "reports", "multiqc", "qc"))
        elif stage == "asm":
            qaa_args.update(**vars(args),  # <- s. above, this is kinda stupid, redesign ASAP!
                            survey_assembly=False,
                            runmode="asm",
                            run_multiqc=True,
                            multiqc_dir=join(args.output_dir, "reports", "multiqc", "asm"))
        elif stage == "ann":
            qaa_args.update(**vars(args),
                            survey_assembly=False,
                            qaa_mode="transcriptome,proteome",
                            runmode="ann")
        else:
            raise ValueError("Invalid stage '{}' in __manage_qaa_args.".format(stage))

        return qaa_args
