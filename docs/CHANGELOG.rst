##########
Change Log
##########

All notable changes to this project are documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


==========
Unreleased
==========

Added
-----
* Support bioinformatics process test case based on Resolwe's
  ``TransactionProcessTestCase``
* Custom version of Resolwe's ``with_resolwe_host`` test decorator which
  skips the decorated tests on non-Linux systems
* Add optimal leaf ordering and simulated annealing to gene and sample
  hierarchical clustering
* Add ``resolwebio/chipseq`` docker image
* Add Odocoileus virginianus texanus (deer) organism to sample
  descriptor
* Add test for import-sra process
* Add Resolwe API to test url configuration
* Add TransactionBioProcessTestCase to support writing tests for
  processes which use resdk
* Add RNA-seq DSS test
* Add custom Cutadapt process

Changed
-------
* Improve variant table name in amplicon report
* Prepend ``api/`` to all URL patterns in the Django test project
* Set hisat2 process resources (memory) to 16GB
* Filter LoFreq output VCF files to remove overlapping indels
* Set core resource requirement in Hisat2 process to 1
* Add `Non-canonical splice sites penalty`, `Disallow soft clipping`
  and `Report alignments tailored specifically for Cufflinks`
  parameters to hisat2 process
* Remove `threads` input from cuffquant and rna-seq workfows
* Set core resource requirement in Cuffquant process to 1
* Use ``resolwebio/chipseq`` docker image in ChIP-Seq processes

Fixed
-----
* Correctly handle paired-end parameters in ``featureCount``
* Fix NaN in explained variance in PCA. When PC1 alone explained more
  than 99% of variance, explained variance for PC2 was not returned
* Fix input sanitization error in ``dss-rna-seq``
* Fix gene source check in hierarchical clustering and PCA
* Enable network access for all import processes
* Fix multiple adapter sequnce (list) input and minimal lenght trimming
  in `cutadapt` process and change the name of the process from
  'cutadatp' to 'cutadapt'
* Fix import-sra process to work with input sanitization fix
* Fix RNA-seq DSS adapters bug
* Allow whitespaces in Cutadapt process input reads names
* Fix sample hierarchical clustering output for a single sample case


================
1.4.1 2017-07-20
================

Changed
-------
* Optionally report all amplicons in Amplicon table

Fixed
-----
* Remove remaining references to calling ``pip`` with
  ``--process-dependency-links`` argument


================
1.4.0 2017-07-04
================

Added
-----
* Amplicon workflow
* Amplicon descriptor schemas
* Amplicon report generator
* Add Rattus norvegicus organism choice to sample schema
* Transforming form Phred 64 to Phred 33 when uploading fastq reads
* Add primertrim process
* RNA-Seq experiment descriptor schema
* iCount sample and reads descriptor schemas
* iCount demultiplexing and sample annotation
* ICount QC
* Add MM8, RN4 and RN6 options to rose2 process
* Add RN4 and RN6 options to bamplot process
* Archive-samples process
* Add bamliquidator
* CheMut workflow
* Dicty primary analysis descriptor schema
* IGV session to Archive-samples process
* Use Resolwe's field projection mixins for knowledge base endpoints
* ``amplicon-table`` process
* Add C. griseus organism choice to Sample descriptor schema
* Add S. tuberosum organism choice to Sample descriptor schema
* Add log2 to gene and sample hierarchical clustering
* Add new inputs to import SRA, add read type selection process
* Set memory resource requirement in jbrowse annotation gff3 and gtf
  processes to 16GB
* Set memory resource requirement in star alignment and index processes
  to 32GB
* Add C. elegans organism choice to Sample descriptor schema
* Add D. melanogaster organism choice to Sample descriptor schema
* Set core resource requirement in Bowtie process to 1
* Set memory resource requirement in amplicon BWA trim process to 32GB
* Add new master file choices to amplicon panel descriptor schema
* Add S. tuberosum organism choice to RNA-seq workflow
* Add Cutadapt process
* Add leaf ordering to gene and sample hierarchical clustering

Fixed
-----
* Use new import paths in ``resolwe.flow``
* Upload reads (paired/single) containing whitespace in the file name
* Fix reads filtering processes for cases where input read file names
  contain whitespace
* Add additional filtering option to STAR aligner
* Fix bbduk-star-htseq_count workflow
* Fix cuffnorm process: Use sample names as labels (boxplot, tables),
  remove group labels input, auto assign group labels, add outputs for
  Rscript output files which were only available compressed
* Derive output filenames in hisat2 from the first reads filename
* Correctly fetch KB features in ``goea.py``
* Append JBrowse tracks to sample
* Replace the BAM MD tag in `align-bwa-trim` process to correct for an
  issue with the primerclip tool
* Fix typo in trimmomatic and bbduk processes
* Use re-import in `etc` and `hmmer_database` processes

Changed
-------
* Support Resolwe test framework
* Run tests in parallel with Tox
* Use Resolwe's new ``FLOW_DOCKER_COMMAND`` setting in test project
* Always run Tox's ``docs``, ``linters`` and ``packaging`` environments
  with Python 3
* Add ``extra`` Tox testing environment with a check that there are no
  large test files in ``resolwe_bio/tests/files``
* Replace Travis CI with Genialis' Jenkins for running the tests
* Store compressed and uncompressed .fasta files in
  ``data:genome:fasta`` objects
* Change sample_geo descriptor schema to have strain option available
  for all organisms
* More readable rna-seq-quantseq schema, field stranded
* Remove obsolete Gene Info processes
* Change log2(fc) default from 2 to 1 in diffexp descriptor schema
* Change Efective genome size values to actual values in macs14 process
* Change variable names in bowtie processes
* Remove iClip processes, tools, files and tests


================
1.3.0 2017-01-28
================

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
* Add stitch parameter to rose2 processor
* Add filtering to DESeq2
* Set Docker Compose's project name to ``resolwebio`` to avoid name clashes
* GO enrichment analysis: map features using gene Knowledge base
* Add option to upload .gff v2 files with upload-gtf processor
* Replace Haystack with Resolwe Elastic Search API
* Require Resolwe 1.4.1+
* Update processes to be compatible with Resolwe 1.4.0

Added
-----
* Process definition documentation style and text improvements
* Add ``resolwe_bio.kb`` app, Resolwe Bioinformatics Knowledge Base
* Add tests to ensure generators produce the same results
* Upload Gene sets (``data:geneset``)
* Add ``generate_geneset`` django-admin command
* Add ``generate_diffexpr_deseq`` django-admin command
* Add 'Generate GO gene sets' processor
* Add generic file upload processors
* Add upload processor for common image file types (.jpg/.tiff/.png/.gif)
* Add upload processor for tabular file formats (.tab/.tsv/.csv/.txt/.xls/.xlsx)
* Add Trimmomatic process
* Add featureCounts process
* Add Subread process
* Add process for hierarchical clusteing of samples
* Add gff3 to gtf file converter
* Add microarray data descriptor schema
* Add process for differential expression edgeR
* ``BioCollectionFilter`` and ``BidDataFilter`` to support filtering collections
  and data by samples on API
* Added processes for automatically downloading single and paired end SRA files
  from NCBI and converting them to FASTQ
* Added process for automatically downloading SRA files from NCBI and converting
  them to FASTQ
* Add HEAT-Seq pipeline tools
* Add HEAT-Seq workflow
* Add ``create-geneset``, ``create-geneset-venn``  processors
* Add ``source`` filter to feature search endpoint
* Add bamplot process
* Add gene hiererhical clustering
* Add cuffquant workflow
* Support Django 1.10 and versionfield 0.5.0
* django-admin commands ``insert_features`` and ``insert_mappings`` for
  importing features and mappings to the Knowledge Base
* Add bsmap and mcall to analyse WGBS data
* Vaccinesurvey sample descriptor schema
* Add RNA-Seq single and paired-end workflow

Fixed
-----
* Set ``presample`` to ``False`` for Samples created on Sample endpoint
* Fix FastQC report paths in processors
* Fix ``htseq_count`` and ``featureCounts`` for large files
* Fix ``upload gtf annotation``
* Fix gene_id field type for differential expression storage objects
* Order data objects in ``SampleViewSet``
* Fix sample hiererhical clustering
* Fix name in gff to gtf process
* Fix clustering to read expressed genes as strings
* Fix protocol labels in ``rna-seq-quantseq`` descriptor schema


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
