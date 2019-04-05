import sys
import os
import argparse
import pathlib
import shutil

CONFIG_FILES = [
	"bgrrl_config.yaml", 
	"hpc_config.json", 
	"multiqc_config.yaml", 
	"enterobase.yaml"
]

def main():

	ap = argparse.ArgumentParser()
	ap.add_argument(		
		"--output-dir", "-o",
		type=str,
		default="Analysis",
		help="Path to output directory. [Analysis]"
	)

	ap.add_argument(		
		"--force-overwrite", "-f", 
		action="store_true",
		help="Force overwriting existing config directory. [False]"
	)
	args = ap.parse_args()

	print("Output directory is {} (exists={})".format(args.output_dir, os.path.exists(args.output_dir)))

	config_dir = os.path.join(args.output_dir, "config")
	if os.path.exists(config_dir):
		print("Configuration directory {} already exists.".format(config_dir))
		if not args.force_overwrite:
			print("Please choose a different output directory or use the -f option to overwrite the configuration.")
			sys.exit(2)
		print("--force-overwrite/-f option found, configuration will be overwritten with fresh copies")			
	else:
		print("Creating configuration directory {} ... ".format(config_dir), end="", flush=True)
		pathlib.Path(config_dir).mkdir(parents=True, exist_ok=True)
		print("done.")

	etc_dir = os.path.join(os.path.dirname(__file__), "..", "etc")
	print("Configuration templates are located in " + etc_dir)
	print("Copying configuration files to {} ... ".format(config_dir), end="", flush=True)

	for f in CONFIG_FILES:
		try:
			shutil.copyfile(os.path.join(etc_dir, f), os.path.join(config_dir, f))
		except:
			raise
		if f == "bgrrl_config.yaml":
			with open(os.path.join(config_dir, f), "at") as fh:
				print("", file=fh)
				print("enterobase_criteria: {}".format(os.path.join(os.path.abspath(config_dir), "enterobase.yaml")), file=fh)

	print("done.")

		

	pass



if __name__ == "__main__":
	main()







































