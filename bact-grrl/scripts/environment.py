class ExecutionEnvironment:
	def __init__(self, args=None, now="eipap", job_suffix=None, log_dir="logs"):
		self.use_drmaa = False
		self.use_scheduler = False
		self.partition = ""
		self.max_nodes = 1
		self.max_cores = 4
		self.sub_cmd = ""
		self.res_cmd = ""
		self.hpc_config = ""
		if args:
			scheduler = args.scheduler if args.scheduler and args.scheduler != '' else "SLURM"
			self.partition = args.partition if args.partition else "{cluster.partition}"
			self.hpc_config = args.hpc_config
			self.use_drmaa = not args.no_drmaa
			self.max_nodes = args.max_nodes
			self.max_cores = args.max_cores
			log_prefix = os.path.join(log_dir, now + "_{rule}_%j")
			job_name = "{rule}_" + job_suffix if job_suffix and job_suffix != "" else "{rule}"
			if scheduler.upper() == "LSF":
				self.sub_cmd = "bsub"
				self.res_cmd = " -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads} -q " + self.partition + " -J " + job_name + " -oo " + log_prefix + ".lsf.log"
				self.use_scheduler = True
			elif scheduler.upper() == "PBS":
				self.sub_cmd = "qsub"
				self.res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q " + self.partition + " -N " + job_name + " -o " + log_prefix + ".pbs.stdout -e " + log_prefix + ".pbs.stderr"
				self.use_scheduler = True
			elif scheduler.upper() == "SLURM":
				self.sub_cmd = "sbatch"
				self.res_cmd = 	" -c {threads}" + \
								" -p " + self.partition + \
								" --exclude={cluster.exclude}" + \
								" --mem={cluster.memory}" + \
								" -J " + job_name + \
								" -o " + log_prefix + ".slurm.log" + \
								" --time={cluster.time}"


				self.use_scheduler = True
			elif scheduler == "" or scheduler.upper() == "NONE":
				pass
			else:
				raise ValueError("Unexpected scheduler configuration.  Check settings.")
