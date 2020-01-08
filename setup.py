# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

name="bgrrl"
version = "1.0"

if sys.version_info.major != 3:
	raise EnvironmentError("""bgrrl is a python module that requires python3, and is not compatible with python2.""")


setup(
	name=name,
	version=version,
	description=description,
	long_description=long_description,
	url="https://github.com/EI-CoreBioinformatics/bgrrl",
	author="Christian Schudoma",
	author_email="christian.schudoma@earlham.ac.uk",
	license="MIT",
	classifiers=[
		"Development Status :: 4 - Beta",
		"Topic :: Scientific Engineering :: Bio/Informatics",
		"License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        'Programming Language :: Python :: 3.4',
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6"
	],
	zip_safe=False,
	keywords="bacterial genomics illumina sequencing assembly annotation",
	packages=find_packages(exclude=["test"]),
	scripts=[
		path.join("bgrrl/bin/slurm", script) for script in ["bgsurvey_sub", "bgasm_sub", "bgann_sub", "bgpack_sub"]
	] + ["qaa/bin/qaa_sub"],
	install_requires=[
		"snakemake>=4.4.0",
		"drmaa",
		"sphinx"
	],
	entry_points={
		"console_scripts": [
			"bgrrl=bgrrl.__main__:main",
			"bgsurvey=bgrrl.__main__:main",
			"bgasm=bgrrl.__main__:main",
			"bgann=bgrrl.__main__:main",
			"bgpack=bgrrl.__main__:main",
			"bginit=bgrrl.bin.bginit:main",
			"asm_report=bgrrl.bin.asm_report:main",
			"bg_datasum=bgrrl.bin.bg_datasum:main",
			"create_samplesheet=bgrrl.bin.create_samplesheet:main",
			"prokka_report=bgrrl.bin.prokka_report:main",
			"asm_wrapper=bgrrl.bin.wrappers.asm_wrapper:main",
			"prokka_wrapper=bgrrl.bin.wrappers.prokka_wrapper:main",
			"ratt_wrapper=bgrrl.bin.wrappers.ratt_wrapper:main",
			"qc2asm=bgrrl.bin.qc2asm:main",
			"qaa=qaa.__main__:main"
		]
	},
	package_data={
		"bgrrl.zzz": ["*.smk*"],
		"qaa.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[
		("etc", glob.glob("etc/*.*")),
		("etc/busco", glob.glob("etc/busco/*")),
		("etc/container_definitions", glob.glob("etc/container_definitions/*")),
	]
)
