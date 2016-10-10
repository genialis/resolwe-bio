##########
Change Log
##########

All notable changes to this project are documented in this file.


==========
Unreleased
==========

Changed
-------
* Add option to save expression JSON to file before saving it to Storage
* Update ``upload-expression`` process
* No longer treat ``resolwe_bio/tools`` as a Python package
* Move processes' test files to the ``resolwe_bio/tests/files`` directory
  to generalize and simplify handling of tests' files
* Update differential expression (DE) processors
* Update ``generate_diffexpr_cuffdiff`` django-admin command
* Save gene_id source to ``output.source`` for DE, expression and related objects
* Refactor ``upload-diffexp`` processor
* Update sample descriptor schema
* Remove obsolete descriptor schemas
* Update processes trimmomatic, gff_to_gtf and featureCount

Added
-----
* Process definition documentation style and text improvements
* Add tests to ensure generators produce the same results
* Upload Gene sets (``data:geneset``)
* Add ``generate_geneset`` django-admin command
* Add ``generate_diffexpr_deseq`` django-admin command
* Add 'Generate GO gene sets' processor
* Add generic file upload processors
* Add upload processor for common image file types (.jpg/.tiff/.png/.gif)
* Add upload processor for tabular file formats (.tab/.tsv/.csv/.txt/.xls/.xlsx)
* Add featureCounts process
* Add Subread process
* Add process for hierarchical clusteing of samples
* Add gff3 to gtf file converter
* Add microarray data descriptor schema
* Add process for differential expression edgeR

Fixed
-----
* Set ``presample`` to ``False`` for Samples created on Sample endpoint


================
1.2.1 2016-07-27
================

Changed
-------
* Update ``resolwe`` requirement


================
1.2.0 2016-07-27
================

Changed
-------
* Decorate all tests that currently fail on Docker with ``skipDockerFailure``
* Require Resolwe's ``master`` git branch
* Put packaging tests in a separate Tox testing environment
* Rename DB user in test project
* Change PostgreSQL port in test project
* Add ROSE2 results parser
* Compute index for HISAT2 aligner on genome upload
* Updated Cuffquant/Cuffnorm tools
* Change ROSE2 enhancer rank plot labels
* Refactor processor syntax
* Move processes tests into ``processes`` subdirectory
* Split ``sample`` API endpoint to ``sample`` for annotated ``Samples``
  and ``presample`` for unannotated ``Samples``
* Rename test project's data and upload directories to ``.test_data`` and
  ``.test_upload``
* Save fastq files to ``lists:basic:file`` field. Refactor related processors.
* Reference genome-index path when running aligners.
* Add pre-computed genome-index files when uploading reference fasta file.
* Include all necessary files for running the tests in source distribution
* Exclude tests from built/installed version of the package
* Move testing utilities from ``resolwe_bio.tests.processes.utils`` to
  ``resolwe_bio.utils.test``
* Update Cuffdiff processor inputs and results table parsing
* Refactor processes to use the updated ``resolwe.flow.executors.run`` command
* Refactor STAR aligner - export expressions as separate objects

Fixed
-----
* Make Tox configuration more robust to different developer environments
* Set ``required: false`` in processor input/output fields where necessary
* Add ``Sample``'s ``Data objects`` to ``Collection`` when ``Sample`` is added
* Fixed/renamed Cufflinks processor field names

Added
-----
* ``skipDockerFailure`` test decorator
* Expand documentation on running tests
* Use Travis CI to run the tests
* Add ``Sample`` model and corresponding viewset and filter
* Add docker-compose command for PostgreSQL
* API endpoint for adding ``Samples`` to ``Collections``
* HISAT2 aligner
* Use Git Large File Storage (LFS) for large test files
* Test for ``generate_samples`` django-admin command
* django-admin command: ``generate_diffexpr``


================
1.1.0 2016-04-18
================

Changed
-------
* Remove obsolete utilities superseded by resolwe-runtime-utils
* Require Resolwe 1.1.0

Fixed
-----
* Update sample descriptor schema
* Include all source files and supplementary package data in sdist

Added
-----
* ``flow_collection: sample`` to processes
* MACS14 processor
* Initial Tox configuration for running the tests
* Tox tests for ensuring high-quality Python packaging
* ROSE2 processor
* django-admin command: ``generate_samples``


================
1.0.0 2016-03-31
================

Changed
-------
* Renamed assertFileExist to assertFileExists
* Restructured processes folder hierarchy
* Removed re-require and hard-coded tools' paths

Fixed
-----
* Different line endings are correctly handled when opening gzipped files
* Fail gracefully if the field does not exist in assertFileExists
* Enabled processor tests (GO, Expression, Variant Calling)
* Enabled processor test (Upload reads with old Illumina QC encoding)
* Made Resolwe Bioinformatics work with Resolwe and Docker

Added
-----
* Import expressions from tranSMART
* Limma differential expression (tranSMART)
* VC filtering tool (Chemical mutagenesis)
* Additional analysis options to Abyss assembler
* API endpoint for Sample
* Initial documentation
