from setuptools import setup, find_packages

name = "bgrrl"
version = "0.1"
release = version

# External python modules that can be gathered from
install_requires = [
    'snakemake',
    'tabulate',
    'biopython',
    'recordclass',
    'ete3',
    'drmaa'
    ]

import os 
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join(path, filename))
    return paths


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
	# packages=["bgrrl"],
	packages=find_packages(exclude=["test"]),
	entry_points={"console_scripts": ["bgrrl=bgrrl.__main__:main", "qc_eval=bgrrl.bin.qc_eval:main"]},
	# entry_points={"console_scripts": ["eipap = eipap.__main__:main",
	#                                  "eipap_cd=eipap.eicd:main",
	#                                  "centrifuge_mm=eipap.illumina.centrifuge_mm:main",
	#								  "centrifuge_classify=eipap.illumina.centrifuge_classify:main",
	#                                  "bcl2fastq_report=eipap.illumina.demux_stats:main",
	#								  "eitransfer=eipap.illumina.transfer:main",
	#								  "eianalysis=eipap.illumina.analysis:main"]},
	# install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
	# test_suite="nose.collector",

	install_requires=[
		'snakemake',
		'tabulate',
		'biopython',
		'recordclass',
		'ete3',
		'drmaa'
    ],
	tests_require = [
		'nose',
    ],
	package_data={
		"bgrrl.zzz": ["*.smk.py"]
	},
	# package_data={
	#	"eipap.illumina": ["*.smk", "*.css"]
	#},
	include_package_data=True,
	# data_files=[("etc", ["etc/bgrrl_config.yaml", "etc/busco-config.ini", "etc/hpc_config.json", "etc/multiqc_config.yaml"]), ("etc/wrappers", ["etc/wrappers/unicycler_wrapper"]), ("etc/util", ["etc/util/busco_init_dir"]), ("etc/data", package_files("etc/data"))],
	data_files=[("etc", ["etc/bgrrl_config.yaml", "etc/busco-config.ini", "etc/hpc_config.json", "etc/multiqc_config.yaml"]), ("etc/wrappers", ["etc/wrappers/unicycler_wrapper"]), ("etc/util", ["etc/util/busco_init_dir"])],
	#scripts=["bin/eicd", "bin/eipap_sub", "bin/eitransfer_sub", "bin/eianalysis_sub"],
	#cmdclass={
	#	'build': build,
        #'install_all': InstallAll
	#}
)

