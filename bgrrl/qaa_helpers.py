
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
	},
	"qc_report": {
		"survey_assembly": True,
		"runmode": "survey",
		"run_blobtools": True,
		"run_qualimap": True,
		"align_reads": "bowtie2",
		"run_busco": True,
		"run_multiqc": True,
	},
	"asm": {
		"survey_assembly": False,
		"runmode": "asm",
		"run_multiqc": True,
	},
	"ann": {
		"survey_assembly": False,
		"qaa_mode": "transcriptome,proteome",
		"runmode": "ann"
	}
}

