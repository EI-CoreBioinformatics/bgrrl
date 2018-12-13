# from os.path import join

# from qaa.qaa_args import QAA_ArgumentsAdapter as QAA_Args

STAGE_QAA_ARGS = {
	"init": {
		"make_input_stream": True,
		"run_blobtools": True,
		"run_qualimap": True,
		"align_reads": "bowtie2",
		"qaa_mode": "genome",
		"busco_db": "bacteria_odb9",
		"run_multiqc": False
	},
	"qc_survey": {
		"survey_assembly": True,
		"runmode": "survey",
		"run_blobtools": False,
		"run_qualimap": False,
		"align_reads": False,
		"run_busco": False,
		# "normalized": not self.no_normalization
	},
	"qc_report": {
		"survey_assembly": True,
		"runmode": "survey",
		# "normalized": not self.no_normalization,
		"run_blobtools": True,
		"run_qualimap": True,
		"align_reads": "bowtie2",
		"run_busco": True,
		"run_multiqc": True,
		# "multiqc_dir": join(self.report_dir, "multiqc", "qc")
	},
	"asm": {
		"survey_assembly": False,
		"runmode": "asm",
		"run_multiqc": True,
		# "multiqc_dir": join(self.report_dir, "multiqc", "asm")
	},
	"ann": {
		"survey_assembly": False,
		"qaa_mode": "transcriptome,proteome",
		"runmode": "ann"
	}
}


#DEFAULT_QAA_ARGS = {
#	"make_input_stream": True,
#	"run_blobtools": True,
#	"run_qualimap": True,
#	"align_reads": "bowtie2",
#	"qaa_mode": "genome",
#	"busco_db": "bacteria_odb9",
#	"run_multiqc": False
#}
#
#class QAA_ArgumentManager(object):
#	@staticmethod
#	def get_qaa_args(args, bgrrl_config_file, hpc_config_file, stage="init"):
#
#		qaa_args = QAA_Args(**DEFAULT_QAA_ARGS)
#		qaa_args.update(
#			quast_mincontiglen=1000,
#			project_prefix=args.project_prefix,
#			config=bgrrl_config_file,
#			hpc_config=hpc_config_file
#		)
#
#		if stage == "init":
#			pass
#		else:
#			qaa_args.update(**vars(args))
#
#			if stage == "qc_survey":
#				qaa_args.update(
#					survey_assembly=True,
#					runmode="survey",
#					run_blobtools=False,
#					run_qualimap=False,
#					align_reads=False,
#					run_busco=False,
#					normalized=not args.no_normalization
#				)
#			elif stage == "qc_report":
#				qaa_args.update(
#					survey_assembly=True,
#					runmode="survey",
#					normalized=not args.no_normalization,
#					run_blobtools=True,
#					run_qualimap=True,
#					align_reads="bowtie2",
#					run_busco=True,
#					run_multiqc=True,
#					multiqc_dir=join(args.output_dir, "reports", "multiqc", "qc")
#				)
#			elif stage == "asm":
#				qaa_args.update(
#					survey_assembly=False,
#					runmode="asm",
#					run_multiqc=True,
#					multiqc_dir=join(args.output_dir, "reports", "multiqc", "asm")
#				)
#			elif stage == "ann":
#				qaa_args.update(
#					survey_assembly=False,
#					qaa_mode="transcriptome,proteome",
#					runmode="ann"
#				)
#			else:
#				raise ValueError("Invalid stage '{}' in __manage_qaa_args.".format(stage))
#
#		return qaa_args
