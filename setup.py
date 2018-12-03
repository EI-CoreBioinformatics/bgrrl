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
version = "0.4.4"

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
	scripts=[path.join("bgrrl/bin", script) for script in ["bgqc", "bgasm", "bgann", "bgfin", "bginit"]],
	install_requires=[
		"snakemake>=4.4.0",
		"drmaa"
	],
	entry_points={
		"console_scripts": [
			"bgrrl=bgrrl.__main__:main",
			"qc_eval=bgrrl.bin.qc_eval:main",
			"asm_report=bgrrl.bin.asm_report:main",
			"ann_report=bgrrl.bin.ann_report:main",
			"bg_datasum=bgrrl.bin.bg_datasum:main",
			"asm_stage_report=bgrrl.bin.asm_stage_report:main",
			"ann_cmp=bgrrl.bin.annocmp:main",
			"create_samplesheet=bgrrl.bin.create_samplesheet:main",
			"prokka_report=bgrrl.bin.prokka_report:main",
			"asm_wrapper=bgrrl.bin.wrappers.asm_wrapper:main",
			"prokka_wrapper=bgrrl.bin.wrappers.prokka_wrapper:main"
		]
	},
	package_data={
		"bgrrl.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[
		("bgrrl/etc", glob.glob("bgrrl/etc/*.*")),
		("bgrrl/etc/wrappers", glob.glob("bgrrl/etc/wrappers/*")),
		("bgrrl/etc/util", glob.glob("bgrrl/etc/util/*"))
	]
)

