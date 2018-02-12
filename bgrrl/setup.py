# coding: utf-8

import sys
import os
import shutil

from setuptools import setup, find_packages, Command
from subprocess import check_call

name = "bgrrl"
version = "0.1"
release = version + ".0"


if sys.version_info.major != 3:
	raise EnvironmentError("""bgrrl is a python module that requires python3, and is not compatible with python2.""")


class InstallAll(Command):
	description = "Build bgrrl, documentation and set python paths"
	user_options = [("prefix=", "d", "directory to install the files to")]

	def initialize_options(self):
		self.prefix = None
	def finalize_options(self):
		pass
	def run(self):
		pythonpath = None
		env = os.environ.copy()
		mqc_cmd = [sys.executable, "setup.py", "install"]

		if self.prefix:
			pythonpath = self.prefix + "/lib/python" + str(sys.version_info.major) + "." + str(sys.version_info.minor) + "/site-packages"
			print()	
			print("Custom install prefix provided: " + self.prefix)
			print("PYTHONPATH=" + pythonpath)
			os.makedirs(pythonpath, exist_ok=True)
			env["PYTHONPATH"] = pythonpath
			mqc_cmd.append("--prefix=" + self.prefix)

		if os.path.exists(".git"):
			print()
			print("Updating git submodules")
			check_call(["git", "submodule", "init"])
			check_call(["git", "submodule", "update"])

			print()
			print("Building qaa")

			check_call(mqc_cmd, cwd=os.path.join("deps", "qaa"), env=env)
		
		print()
		print("Building bgrrl")
		check_call(mqc_cmd, env=env)

		print()

# External python modules that can be gathered from
install_requires = [
    'snakemake',
    'biopython',
    'drmaa'
    ]

setup(
	name=name,
	version=release,
	description="The Bacterial Genome Reconstruction and Recognition Pipeline (BGRR|) for the Earlham Institute",
	long_description='''The Bacterial Genome Reconstruction and Recognition Pipeline (BGRR|) is a bacterial genome assembly and annotation pipeline for Illumina paired end sequencing data.''',
	author="Christian Schudoma",
	author_email="christian.schudoma@earlham.ac.uk",
	license="MIT",
	zip_safe=False,
	keywords="bacterial genomics illumina sequencing assembly annotation",
	packages=find_packages(exclude=["test"]),
	entry_points={"console_scripts": ["bgrrl=bgrrl.__main__:main", 
                                          "qc_eval=bgrrl.bin.qc_eval:main",
                                          "asm_report=bgrrl.bin.asm_report:main"]},
	test_suite = "nose.collector",
	install_requires=install_requires,
	tests_require = [
		'nose',
    ],
	package_data={
		"bgrrl.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[("etc", ["etc/bgrrl_config.yaml", "etc/busco-config.ini", "etc/hpc_config.json", "etc/multiqc_config.yaml"]), ("etc/wrappers", ["etc/wrappers/unicycler_wrapper", "etc/wrappers/asm_wrapper", "etc/wrappers/prokka_wrapper", "etc/wrappers/ratt_wrapper"]), ("etc/util", ["etc/util/busco_init_dir"])],
	cmdclass={
        	'install_all': InstallAll
	}
)

