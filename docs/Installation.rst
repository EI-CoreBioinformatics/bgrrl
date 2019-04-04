.. bgrrl documentation master file, created by
   sphinx-quickstart on Thu Apr  4 12:14:34 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
=================================

Requirements
------------

bgrr| is a Python3 application, so Python3.5+ is necessary. 

Python dependencies
^^^^^^^^^^^^^^^^^^^

* snakemake >= 4.4.0
* drmaa
* sphinx
* qaa ()

Third party software
^^^^^^^^^^^^^^^^^^^^

* bbmap
* fastqc
* kat
* unicycler/spades
* velvet-optimizer
* prokka
* emboss
* ratt

**For qaa:**

* quast
* blobtools
* qualimap
* busco
* multiqc ()
* picardtools
* samtools
* bwa and/or bowtie2

**We provide Singularity recipes (one for bgrr| and one for qaa) that cover all the requirements.**

**If you plan to rely on your own installations, please make sure you have Python3-compatible versions, e.g. of quast and blobtools.**

Other resources
^^^^^^^^^^^^^^^

* bacteria.odb9 busco database

  http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz

* a comprehensive nucleotide sequence database (including taxonomy ids) for blobtools

  e.g. NCBI nt


