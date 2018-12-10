from collections import OrderedDict, namedtuple

ExecutionEnvironmentArguments = namedtuple(
	"ExecutionEnvironmentArguments",
	[
		"scheduler",
		"partition",
		"no_drmaa",
		"max_nodes",
		"max_cores",
		"hpc_config"
	]
)


class ConfigurationManager(OrderedDict):
	def __make_exe_env_args(self, ap_args):
		self.exe_env = ExecutionEnvironmentArguments(
			ap_args.scheduler,
			ap_args.partition,
			ap_args.no_drmaa,
			ap_args.max_nodes,
			ap_args.max_cores,
			ap_args.hpc_config
		)

	def __init__(self, ap_args):		
		# hand this over to ExecutionEnvironment
		# or: TODO - lets ConfigurationManager create an ExeEnv
		self.__make_exe_env_args(ap_args)


	def __str__(self):
		return super(ConfigurationManager, self).__str__() + "\n" + str(self.exe_env)

	def setConfiguration(self):
		pass



class BGRRLConfigurationManager(ConfigurationManager):
	def __init__(self, ap_args):

		ap_args.alt_hpc_config_warning = "Please run bginit or provide a valid HPC configuration file with --hpc_config."
        ap_args.alt_config_warning = "Please run bginit or provide a valid configuration file with --bgrrl_config/--config."

		super(BGRRLConfigurationManager, self).__init__(ap_args)	

		self["input_sheet"] = ap_args.input_sheet
		self["output_dir"] = ap_args.output_dir
		self["project_prefix"] = ap_args.project_prefix
		self["config"] = ap_args.config
		self["report_only"] = ap_args.report_only
		self["force"] = ap_args.force
		self["enterobase_groups"] = ap_args.enterobase_groups
		self["unlock"] = ap_args.unlock
		self["runmode"] = ap_args.runmode

		if self["runmode"] == "survey":			
			self["no_normalization"] = ap_args.no_normalization
			self["no_packaging"] = ap_args.no_packaging
			self["full_qaa_analysis"] = ap_args.full_qaa_analysis
			self["minimum_survey_assembly_size"] = ap_args.minimum_survey_assembly_size

		if self["runmode"] == "assembly":
			self["assembler"] = ap_args.assembler
			self["contig_minlen"] = ap_args.contig_minlen
			self["is_final_step"] = ap_args.is_final_step
			self["no_packaging"] = ap_args.no_packaging

		if self["runmode"] == "annotation":
			self["ratt_reference"] = ap_args.ratt_reference
			self["no_packaging"] = ap_args.no_packaging
			self["prokka_package_style"] = ap_args.prokka_package_style

		if self["runmode"] == "package":
			self["package_mode"] = ap_args.package_mode
			self["prokka_package_style"] = ap_args.prokka_package_style
			self["make_ratt_data_tarballs"] = ap_args.make_ratt_data_tarballs



"""
no_normalization        False
no_packaging    False
full_qaa_analysis       False
minimum_survey_assembly_size    1000000.0
input_sheet     samplesheet.csv.20
output_dir      Analysis
project_prefix  bgrrl_test
config  Analysis/config/bgrrl_config.yaml
report_only     False
force   False
enterobase_groups       None
unlock  False
runmode survey
"""	
		
