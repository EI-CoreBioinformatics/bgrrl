.. bgrrl documentation master file, created by
   sphinx-quickstart on Thu Apr  4 12:14:34 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



The survey module
======================

.. image:: ../bgrrl_workflow.png
    :align: center

The survey module has two main functions. 1) Preprocessing of the library data and 2) assessment of the assemble-ability of each sample/library.

In detail, the steps performed by the survey module are:


1. Preprocessing with ``bbduk``

   Each paired-end library is preprocessed with bbduk with the following operations:

    * k-mer based adapter trimming against ``bbduk``'s default adapter list

    * gentle quality trimming on either side, keeping only bases with at least phred=3

    * length filtering (100bp)

    * quality filtering (maq=20)

    * tpe/tbo (s. ``bbduk`` documentation)

   These operations can be adjusted via the ``bgrrl_config.yaml``.


2. Normalization with ``bbnorm`` (optional, switch off with ``--no-normalization``)

   Libraries are normalized to the range of 2x-100x.

   The range can be adjusted via the ``bgrrl_config.yaml``.


3. Read quality assessment with ``fastqc``

4. Read sequence feature assessment with ``kat``

5. Survey assembly with ``tadpole``

   Each library is subjected to a simple and quick assembly with ``tadpole``. This determines if the library can be assembled or not.
   Library assembly statistics are then calculated with ``quast`` via ``qaa``.

6. Filtering by assemble-ability

   Libraries that could be assembled with ``tadpole`` and that contain at least 1000 reads and have an assembly size above 1Mbp (or user-specified, s. below) 
   will be automatically passed on to the assembly stage. Libraries that fail these checks will be filtered out. However, the user can manually add them
   to the assembly samplesheet (s. below).

7. Assembly samplesheet generation

   Samples that passed the filtering stage will be written to a new samplesheet (``samplesheet.qc_pass.csv``) in the ``<outdir>/reports/samplesheets`` directory. 
   This samplesheet can then be used as input for the assembly stage.

8. Read packaging

   Preprocessed reads will be automatically packaged into the ``<outdir>/Data_Package`` directory unless suppressed by the ``--no-packaging`` option.



Command line arguments
----------------------

The command ``bgrrl -h`` or ``bgrrl <stage> -h`` (or ``--help`` instead of ``-h``) will display a list of command line options.

**Remember each bgrr| run requires at the very least the following three command line parameters:**

* ``input_sheet``
* ``--config``
* ``--hpc_config``

survey options:
^^^^^^^^^^^^^^^

* ``--no-normalization``

  Disable read normalization. [False]

* ``--no-packaging``

  Disable automatic packaging. [False]

* ``--full-qaa-analysis``

  Perform full qaa-analysis on survey assemblies. [False]

* ``--minimum-survey-assembly-size MINIMUM_SURVEY_ASSEMBLY_SIZE``

  Minimum size (in bp) for tadpole assembly to pass
  survey stage [1Mbp] Setting this option to a smaller size allows plasmid-specific libraries to be
  processed.















