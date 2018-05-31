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
- Add to ``resolwebio/chipseq`` Docker image:

  - ``Bedops (v2.4.32)``
  - ``Tabix (v1.8)``
  - ``python3-pandas``
  - ``bedGraphToBigWig (kent-v365)``
  - ``bedToBigBed (kent-v365)``
- Add to ``resolwebio/rnaseq:3.1.0`` Docker image:

  - ``genometools (1.5.9)``
  - ``igvtools (v2.3.98)``
  - ``jbrowse (v1.12.0)``
- Add to ``resolwebio/rnaseq:3.2.0`` Docker image:

  - ``Bowtie (v1.2.2)``
  - ``Bowtie2 (v2.3.4.1)``
  - ``BWA (0.7.17-r1188)``
  - ``TopHat (v2.1.1)``
  - ``Picard Tools (v2.18.5)``
  - ``bedGraphToBigWig (kent-v365)``
- Support filtering by type on feature API endpoint

Changed
-------
- **BACKWARD INCOMPATIBLE:** Drop support for Python 3.4 and 3.5
- **BACKWARD INCOMPATIBLE:** Require Resolwe 10.x
- Update ``wigToBigWig`` to kent-v365 version  in ``resolwebio/chipseq``
  Docker image
- Change paths in HTML amplicon report template in ``resolwebio/dnaseq`` Docker
  image

Fixed
-----
- Make management commands use a private random generator instance
- Fix output ``covplot_html`` of ``coveragebed`` process


==================
9.0.0 - 2018-05-15
==================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Simplify the ``amplicon-report`` process inputs
  by using Latex report template from the ``resolwebio/latex`` Docker image assets
- **BACKWARD INCOMPATIBLE:** Simplify the ``coveragebed`` process inputs
  by using Bokeh assets from the ``resolwebio/dnaseq`` Docker image
- **BACKWARD INCOMPATIBLE:** Require Resolwe 9.x
- Update ``wigToBigWig`` tool in ``resolwebio/chipseq`` Docker image
- Use ``resolwebio/rnaseq:3.1.0`` Docker image in the following
  processes:

  - ``cufflinks``
  - ``cuffnorm``
  - ``cuffquant``
- Remove ``differentialexpression-limma`` process
- Use ``resolwebio/rnaseq:3.1.0`` docker image and expand error
  messages in:

  - ``cuffdiff``
  - ``differentialexpression-deseq2``
  - ``differentialexpression-edger``
- Update ``workflow-bbduk-star-htseq``
- Update ``quantseq`` descriptor schema
- Assert species and build in ``htseq-count-normalized`` process
- Set amplicon report template in ``resolwebio/latex`` Docker image to
  landscape mode

Added
-----
- Support Python 3.6
- Add ``template_amplicon_report.tex`` to ``resolwebio/latex`` Docker image
  assets
- Add SnpEff tool and bokeh assets to ``resolwebio/dnaseq`` Docker image
- Add automated library strand detection to ``feature_counts`` quantification process
- Add FastQC option ``nogroup`` to ``bbduk-single`` and ``bbduk-paired`` processes
- Add CPM normalization to ``htseq-count-raw`` process
- Add ``workflow-bbduk-star-htseq-paired``
- Add legend to amplicon report template in ``resolwebio/latex`` Docker image

Fixed
-----
- Fix manual installation of packages in Docker images to handle dots and
  spaces in file names correctly
- Fix COSMIC url template in ``amplicon-table`` process
- Fix Create IGV session in Archive samples process
- Fix ``source`` tracking in ``cufflinks`` and ``cuffquant`` processes
- Fix amplicon master file validation script. Check and report error if
  duplicated amplicon names are included. Validation will now pass also
  for primer sequences in lowercase.
- Fix allele frequency (AF) calculation in ``snpeff`` process
- Fix bug in script for calculating FPKM. Because genes of raw counts from
  ``featureCounts`` were not lexicographically sorted, division of normalized counts
  was done with values from other, incorrect, genes. Results from ``featureCounts``,
  but not ``HTSeq-count`` process, were affected.


==================
8.1.0 - 2018-04-13
==================

Changed
-------
- Use the latest versions of the following Python packages in
  ``resolwebio/rnaseq`` docker image: Cutadapt 1.16, Apache Arrow 0.9.0, pysam
  0.14.1, requests 2.18.4, appdirs 1.4.3, wrapt 1.10.11, PyYAML 3.12
- Bump tools version in ``resolwebio/rnaseq`` docker image:

  - Salmon to 0.9.1
  - FastQC to 0.11.7
- Generalize the no-extraction-needed use-case in ``resolwebio/base`` Docker
  image ``download_and_verify`` script

Added
-----
- Add the following Python packages to ``resolwebio/rnaseq`` docker image: six
  1.11.0, chardet 3.0.4, urllib3 1.22, idna 2.6, and certifi 2018.1.18
- Add ``edgeR`` R library to ``resolwebio/rnaseq`` docker image
- Add Bedtools to ``resolwebio/rnaseq`` docker image

Fixed
-----
- Handle filenames with spaces in the following processes:

  - ``alignment-star-index``
  - ``alignment-tophat2``
  - ``cuffmerge``
  - ``index-fasta-nucl``
  - ``upload-fasta-nucl``
- Fix COSMIC url template in (multisample) amplicon reports


==================
8.0.0 - 2018-04-11
==================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Refactor ``trimmomatic-single``,
  ``trimmomatic-paired``, ``bbduk-single``, and ``bbduk-paired`` processes
- **BACKWARD INCOMPATIBLE:** Merge ``align-bwa-trim`` and ``align-bwa-trim2``
  process functionality. Retain only the refactored process under slug
  ``align-bwa-trim``
- **BACKWARD INCOMPATIBLE:** In processes handling VCF files, the output
  VCF files are stored in bgzip-compressed form. Tabix index is not referenced
  to an original VCF file anymore, but stored in a separate ``tbi`` output
  field
- **BACKWARD INCOMPATIBLE:** Remove an obsolete ``workflow-accel-2`` workflow
- **BACKWARD INCOMPATIBLE:** Use Elasticsearch version 5.x
- **BACKWARD INCOMPATIBLE:** Parallelize execution of the following processes:

  - ``alignment-bowtie2``
  - ``alignment-bwa-mem``
  - ``alignment-hisat2``
  - ``alignment-star``
  - ``alignment-tophat2``
  - ``cuffdiff``
  - ``cufflinks``
  - ``cuffquant``
- Require Resolwe 8.x
- Bump STAR aligner version in ``resolwebio/rnaseq`` docker image to 2.5.4b
- Bump Primerclip version in ``resolwebio/dnaseq`` docker image
- Use ``resolwebio/dnaseq`` Docker image in ``picard-pcrmetrics`` process
- Run ``vc-realign-recalibrate`` process using multiple cpu cores to optimize
  the processing time
- Use ``resolwebio/rnaseq`` Docker image in ``alignment-star`` process

Added
-----
- Add CNVKit, LoFreq and GATK to ``resolwebio/dnaseq`` docker image
- Add BaseSpace files download tool
- Add process to import a file from BaseSpace
- Add process to convert files to single-end reads
- Add process to convert files to paired-end reads
- Add ``vc-gatk4-hc`` process which implements GATK4 HaplotypeCaller variant
  calling tool
- Add ``workflow-accel-gatk4`` pipeline that uses GATK4 HaplotypeCaller as an
  alternative to GATK3 used in ``workflow-accel`` pipeline
- Add ``amplicon-master-file`` descriptor schema
- Add ``workflow-bbduk-star-featurecounts`` pipeline
- Add ``rna-seq-bbduk-star-featurecounts`` RNA-seq descriptor schema

Fixed
-----
- Fix iterative trimming in ``bowtie`` and ``bowtie2`` processes
- Fix ``archive-samples`` to use sample names for headers when merging
  expressions
- Improve ``goea.py`` tool to handle duplicated mapping results
- Handle filenames with spaces in the following processes:

  - ``alignment-hisat2``
  - ``alignment-bowtie``
  - ``prepare-geo-chipseq``
  - ``prepare-geo-rnaseq``
  - ``cufflinks``
  - ``cuffquant``


==================
7.0.1 - 2018-03-27
==================

Fixed
-----
* Use name-ordered BAM file for counting reads in ``HTSeq-count`` process by
  default to avoid buffer overflow with large BAM files


==================
7.0.0 - 2018-03-13
==================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Remove Ubuntu 17.04 base Docker image since it has
  has reached its end of life and change all images to use the new ubuntu 17.10
  base image
- **BACKWARD INCOMPATIBLE:** Require ``species`` and ``build`` inputs in the
  following processes:

  - ``upload-genome``
  - ``upload-gtf``
  - ``upload-gff3``
  - ``upload-bam``
  - ``upload-bam-indexed``
- **BACKWARD INCOMPATIBLE:** Track ``species`` and ``build`` information in the
  following processes:

  - ``cuffmerge``
  - alignment processes
  - variant calling processes
  - JBrowse processes
- **BACKWARD INCOMPATIBLE:** Track ``species``, ``build`` and ``feature_type``
  in the following processes:

  - ``upload-expression-star``
  - quantification processes
  - differential expression processes
- **BACKWARD INCOMPATIBLE:** Track ``species`` in gene set (Venn) and
  ``goenrichment`` processes
- **BACKWARD INCOMPATIBLE:** Rename ``genes_source`` input to ``source`` in
  hierarchical clustering and PCA processes
- **BACKWARD INCOMPATIBLE:** Remove the following obsolete processes:

  - Dictyostelium-specific ncRNA quantification
  - ``go-geneset``
  - bayseq differential expression
  - ``cuffmerge-gtf-to-gff3``
  - ``transdecoder``
  - ``web-gtf-dictybase``
  - ``upload-rmsk``
  - ``snpdat``
- **BACKWARD INCOMPATIBLE:** Unify output fields of processes of type
  ``data:annotation``
- **BACKWARD INCOMPATIBLE:** Rename the organism field names to species in
  ``rna-seq`` and ``cutadapt-star-htseq`` descriptor schemas
- **BACKWARD INCOMPATIBLE:** Rename the ``genome_and_annotation`` field name
  to ``species`` in ``bcm-*`` descriptor schemas and use the full species name
  for the ``species`` field values
- **BACKWARD INCOMPATIBLE:** Refactor ``featureCounts`` process
- **BACKWARD INCOMPATIBLE:** Change ``import-sra`` process to work with
  ``resolwebio/utils`` Docker image and refactor its inputs
- Require Resolwe 7.x
- Add environment export for Jenkins so that the manager will use a
  globally-unique channel name
- Set ``scheduling_class`` of gene and sample hierarchical clustering processes
  to ``interactive``
- Change base Docker images of ``resolwebio/rnaseq`` and ``resolwebio/dnaseq``
  to ``resolwebio/base:ubuntu-18.04``
- Use the latest versions of the following Python packages in
  ``resolwebio/rnaseq`` Docker image: Cutadapt 1.15, Apache Arrow 0.8.0,
  pysam 0.13, and xopen 0.3.2
- Use the latest versions of the following Python packages in
  ``resolwebio/dnaseq`` Docker image: Bokeh 0.12.13, pandas 0.22.0,
  Matplotlib 2.1.2, six 1.11.0, PyYAML 3.12, Jinja2 2.10, NumPy 1.14.0,
  Tornado 4.5.3, and pytz 2017.3
- Use the latest version of ``wigToBigWig`` tool in ``resolwebio/chipseq``
  Docker image
- Use ``resolwebio/rnaseq:3.0.0`` Docker image in ``goenrichment``,
  ``upload-gaf`` and ``upload-obo`` processes
- Use ``resolwebio/dnaseq:3.0.0`` Docker image in ``filtering_chemut`` process
- Change ``cuffnorm`` process type to ``data:cuffnorm``
- Set type of ``coverage-garvan`` process to ``data:exomecoverage``
- Remove ``gsize`` input from ``macs14`` process and automate genome size
  selection
- Adjust ``bam-split`` process so it can be included in workflows
- Make ID attribute labels in ``featureCounts`` more informative
- Change 'source' to 'gene ID database' in labes and descriptions
- Change ``archive-samples`` process to create different IGV session files for
  ``build`` and ``species``
- Expose advanced parameters in Chemical Mutagenesis workflow
- Clarify some descriptions in the ``filtering_chemut`` process and ``chemut``
  workflow
- Change expected genome build formatting for hybrid genomes in ``bam-split``
  process
- Set the ``cooksCutoff`` parameter to ``FALSE`` in ``deseq.R`` tool
- Rename 'Expressions (BCM)' to 'Dicty expressions'

Added
-----
- Mechanism to override the manager's control channel prefix from the
  environment
- Add Ubuntu 17.10 and Ubuntu 18.04 base Docker images
- Add ``resolwebio/utils`` Docker image
- Add ``BBMap``, ``Trimmomatic``, ``Subread``, ``Salmon``, and
  ``dexseq_prepare_annotation2`` tools and ``DEXSeq`` and ``loadSubread`` R
  libraries to ``resolwebio/rnaseq`` Docker image
- Add abstract processes that ensure that all processes that inherit from them
  have the input and output fields that are defined in them:

  - ``abstract-alignment``
  - ``abstract-annotation``
  - ``abstract-expression``
  - ``abstract-differentialexpression``
  - ``abstract-bed``
- Add miRNA workflow
- Add ``prepare-geo-chipseq`` and ``prepare-geo-rnaseq`` processes that produce
  a tarball with necessary data and folder structure for GEO upload
- Add ``library-strandedness`` process which uses the ``Salmon`` tool built-in
  functionality to detect the library strandedness information
- Add ``species`` and ``genome build`` output fields to ``macs14`` process
- Expose additional parameters in ``alignment-star``, ``cutadapt-single`` and
  ``cutadapt-paired`` processes
- Add ``merge expressions`` to ``archive-samples`` process
- Add description of batch mode to Expression aggregator process
- Add error and warning messages to the ``cuffnorm`` process
- Add optional ``species`` input to hierarchical clustering and PCA processes
- Add Rattus norvegicus species choice to the ``rna-seq`` descriptor schema
  to allow running RNA-seq workflow for this species from the Recipes

Fixed
-----
- Fix custom argument passing script for ``Trimmomatic`` in
  ``resolwebio/rnaseq`` Docker image
- Fix installation errors for ``dexseq-prepare-annotation2`` in
  ``resolwebio/rnaseq`` Docker image
- Fix ``consensus_subreads`` input option in Subread process
- Limit variant-calling process in the chemical mutagenesis workflow and the
  Picard tools run inside to 16 GB of memory to prevent them from crashing
  because they try to use too much memory
- The chemical mutagenesis workflow was erroneously categorized as
  ``data:workflow:rnaseq:cuffquant`` type. This is switched to
  ``data:workflow:chemut`` type.
- Fix handling of NA values in Differential expression results table. NA values
  were incorrectly replaced with value 0 instead of 1
- Fix ``cuffnorm`` process to work with samples containing dashes in
  their name and dispense prefixing sample names starting with numbers
  with 'X' in the ``cuffnorm`` normalization outputs
- Fix ``cuffnorm`` process' outputs to correctly track species and
  build information
- Fix typos and sync parameter description common to ``featureCounts``
  and ``miRNA`` workflow


==================
6.2.2 - 2018-02-21
==================

Fixed
-----
- Fix ``cuffnorm`` process to correctly use sample names as labels in output
  files and expand ``cuffnorm`` tests


==================
6.2.1 - 2018-01-28
==================

Changed
-------
- Update description text of ``cutadapt-star-htseq`` descriptor schema to
  better describe the difference between gene/transcript-type analyses
- Speed-up management command for inserting mappings


==================
6.2.0 - 2018-01-17
==================

Added
-----
- Add R, tabix, and CheMut R library to ``resolwebio/dnaseq`` Docker image
- Add ``SRA Toolkit`` to ``resolwebio/rnaseq`` Docker image

Changed
-------
- Require Resolwe 6.x
- Extend pathway map with species and source field
- Move template and logo for multi-sample report into ``resolwebio/latex``
  Docker image
- Refactor ``amplicon-report`` process to contain all relevant inputs for
  ``amplicon-archive-multi-report``
- Refactor ``amplicon-archive-multi-report``
- Use ``resolwebio/dnaseq:1.2.0`` Docker image in ``filtering_chemut`` process

Fixed
-----
- Enable DEBUG setting in tests using Django's ``LiveServerTestCase``
- Wait for ElasticSeach to index the data in ``KBBioProcessTestCase``
- Remove unused parameters in TopHat (2.0.13) process and Chip-seq workflow


==================
6.1.0 - 2017-12-12
==================

Added
-----
- Add ``amplicon-archive-multi-report`` process
- Add ``upload-metabolic-pathway`` process
- Add memory-optimized primerclip as a separate ``align-bwa-trim2`` process
- Add ``workflow-accel-2`` workflow

Changed
-------
- Improve ``PCA`` process performance
- Use ``resolwebio/chipseq:1.1.0`` Docker image in ``macs14`` process
- Change formatting of ``EFF[*].AA`` column in ``snpeff`` process
- Save unmapped reads in ``aligment-hisat2`` process
- Turn off test profiling

Fixed
-----
- Fix pre-sorting in ``upload-master-file`` process
- Revert ``align-bwa-trim`` process to use non-memory-optimized primerclip
- Fix file processing in ``cutadapt-custom-paired`` process


==================
6.0.0 - 2017-11-28
==================

Added
-----
- Add AF filter to amplicon report
- Add number of samples to the output of expression aggregator
- Add ``ChIP-Rx``, ``ChIPmentation`` and ``eClIP`` experiment types to
  ``reads`` descriptor schema
- Add ``pandas`` Python package to ``resolwebio/latex`` Docker image
- Add primerclip, samtools, picard-tools and bwa to ``resolwebio/dnaseq``
  Docker image
- Add ``cufflinks``, ``RNASeqT`` R library, ``pyarrow`` and ``sklearn`` Python
  packages to ``resolwebio/rnaseq`` Docker image
- Add ``wigToBigWig`` tool to ``resolwebio/chipseq`` Docker image

Changed
-------
- **BACKWARD INCOMPATIBLE:** Drop Python 2 support, require Python 3.4 or 3.5
- **BACKWARD INCOMPATIBLE:** Make species part of the feature primary key
- **BACKWARD INCOMPATIBLE:** Substitute Python 2 with Python 3 in
  ``resolwebio/rnaseq`` Docker image. The processes to be updated to this
  version of the Docker image should also have their Python scripts updated to
  Python 3.
- Require Resolwe 5.x
- Set maximum RAM requirement in ``bbduk`` process
- Move *Assay type* input parameter in RNA-Seq descriptor schema from advanced
  options to regular options
- Use ``resolwebio/rnaseq`` Docker image in Cutadapt processes
- Use additional adapter trimming option in ``cutadapt-custom-single/paired``
  processes
- Show antibody information in ``reads`` descriptor for ``ChIP-Seq``,
  ``ChIPmentation``, ``ChIP-Rx``,  ``eClIP``, ``MNase-Seq``, ``MeDIP-Seq``,
  ``RIP-Seq`` and ``ChIA-PET`` experiment types
- Use ``resolwebio/dnaseq`` Docker image in ``align-bwa-trim`` process
- Refactor ``resolwebio/chipseq`` Docker image
- Use Resolwe's Test Runner for running tests and add ability to only run a
  partial test suite based on what proceses have Changed
- Configure Jenkins to only run a partial test suite when testing a pull
  request
- Make tests use the live Resolwe API host instead of external server

Fixed
-----
- Fix merging multiple expressions in DESeq process
- Fix ``resolwebio/rnaseq`` Docker image's README
- Handle multiple ALT values in amplicon report
- Fix BAM file input in ``rsem`` process


==================
5.0.1 - 2017-11-14
==================

Fixed
-----
- Update Features and Mappings ElasticSearch indices building to be compatible
  with Resolwe 4.0


==================
5.0.0 - 2017-10-25
==================

Added
-----
- Add automatic headers extractor to ``bam-split`` process
- Add HTML amplicon plot in ``coveragebed`` process
- Add raw RSEM tool output to `rsem` process output
- Add support for transcript-level differential expression
  in ``deseq2`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Bump Django requirement to version 1.11.x
- **BACKWARD INCOMPATIBLE:** Make ``BioProcessTestCase`` non-transactional
- Require Resolwe 4.x
- Add the advanced options checkbox to the ``rna-seq`` descriptor schema
- Remove static amplicon plot from ``coveragebed`` and ``amplicon-report``
  processes
- Update Dockerfile for ``resolwebio/latex`` with newer syntax and add some
  additional Python packages


==================
4.2.0 - 2017-10-05
==================

Added
-----
- Add ``resolwebio/base`` Docker image based on Ubuntu 17.04
- Add ``resolwebio/dnaseq`` Docker image
- Add ``DESeq2`` tool to ``resolwebio/rnaseq`` docker image
- Add input filename regex validator for ``upload-master-file`` process

Changed
-------
- Remove obsolete mongokey escape functionality
- Report novel splice-site junctions in HISAT2
- Use the latest stable versions of the following bioinformatics
  tools in ``resolwebio/rnaseq`` docker image: Cutadapt 1.14,
  FastQC 0.11.5, HTSeq 0.9.1, and SAMtools 1.5


==================
4.1.0 - 2017-09-22
==================

Added
-----
- Add Mus musculus to all BCM workflows' schemas
- Add ``bam-split`` process with supporting processes
  ``upload-bam-primary``, ``upload-bam-secondary`` and
  ``upload-header-sam``

Changed
-------
- Enable Chemut workflow and process tests

Fixed
-----
- Fix chemut ``intervals`` input option


==================
4.0.0 - 2017-09-14
==================

Added
-----
- New base and legacy Docker images for processes, which support non-root
  execution as implemented by Resolwe

Changed
-------
- **BACKWARD INCOMPATIBLE:** Modify all processes to explicitly use the new Docker images
- **BACKWARD INCOMPATIBLE:** Remove ``clustering-hierarchical-genes-etc`` process
- Require Resolwe 3.x


================
3.2.0 2017-09-13
================

Added
-----
- Add ``index-fasta-nucl`` and ``rsem`` process
- Add custom Cutadapt - STAR - RSEM workflow


================
3.1.0 2017-09-13
================

Added
-----
- Add statistics of logarithmized expressions to ``expression-aggregator``
- Add input field description to ``cutadapt-star-htseq`` descriptor schema
- Add ``HISAT2`` and ``RSEM`` tool to ``resolwebio/rnaseq`` docker image

Changed
-------
- Remove ``eXpress`` tool from ``resolwebio/rnaseq`` docker image
- Use system packages of RNA-seq tools in ``resolwebio/rnaseq`` docker image
- Set ``hisat2`` process' memory resource requirement to 32GB
- Use ``resolwebio/rnaseq`` docker image in ``hisat2`` process


================
3.0.0 2017-09-07
================

Added
-----
- Add custom Cutadapt - STAR - HT-seq workflow
- Add expression aggregator process
- Add ``resolwebio/rnaseq`` docker image
- Add ``resolwebio/latex`` docker image
- Add access to sample field of data objects in processes via ``sample`` filter

Changed
-------
- **BACKWARD INCOMPATIBLE:** Remove ``threads`` input in STAR aligner process
  and replace it with the ``cores`` resources requirement
- **BACKWARD INCOMPATIBLE:** Allow upload of custom amplicon master files (make
  changes to ``amplicon-panel`` descriptor schema, ``upload-master-file`` and
  ``amplicon-report`` processes and ``workflow-accel`` workflow)
- **BACKWARD INCOMPATIBLE:** Remove ``threads`` input in ``cuffnorm`` process
  and replace it with the ``cores`` resources requirement
- Add sample descriptor to ``prepare_expression`` test function
- Prettify amplicon report

Fixed
-----
- Fix ``upload-expression-star`` process to work with arbitrary file names
- Fix STAR aligner to work with arbitrary file names
- Fix ``cuffnorm`` group analysis to work correctly
- Do not crop Amplicon report title as this may result in malformed LaTeX
  command
- Escape LaTeX's special characters in ``make_report.py`` tool
- Fix validation error in ``Test sleep progress`` process


================
2.0.0 2017-08-25
================

Added
-----
- Support bioinformatics process test case based on Resolwe's
  ``TransactionProcessTestCase``
- Custom version of Resolwe's ``with_resolwe_host`` test decorator which skips
  the decorated tests on non-Linux systems
- Add optimal leaf ordering and simulated annealing to gene and sample
  hierarchical clustering
- Add ``resolwebio/chipseq`` docker image and use it in ChIP-Seq processes
- Add Odocoileus virginianus texanus (deer) organism to sample descriptor
- Add test for ``import-sra`` process
- Add RNA-seq DSS test
- Add Cutadapt and custom Cutadapt processes

Changed
-------
- Require Resolwe 2.0.x
- Update processes to support new input sanitization introduced in Resolwe
  2.0.0
- Improve variant table name in amplicon report
- Prepend ``api/`` to all URL patterns in the Django test project
- Set ``hisat2`` process' memory resource requirement to 16GB and cores
  resource requirement to 1
- Filter LoFreq output VCF files to remove overlapping indels
- Add `Non-canonical splice sites penalty`, `Disallow soft clipping` and
  `Report alignments tailored specifically for Cufflinks` parameters to
  ``hisat2`` process
- Remove ``threads`` input from ``cuffquant`` and ``rna-seq`` workfows
- Set core resource requirement in ``cuffquant`` process to 1

Fixed
-----
- Correctly handle paired-end parameters in ``featureCount``
- Fix ``NaN`` in explained variance in PCA. When PC1 alone explained more than
  99% of variance, explained variance for PC2 was not returned
- Fix input sanitization error in ``dss-rna-seq`` process
- Fix gene source check in hierarchical clustering and PCA
- Enable network access for all import processes
- Fix RNA-seq DSS adapters bug
- Fix sample hierarchical clustering output for a single sample case


================
1.4.1 2017-07-20
================

Changed
-------
- Optionally report all amplicons in Amplicon table

Fixed
-----
- Remove remaining references to calling ``pip`` with
  ``--process-dependency-links`` argument


================
1.4.0 2017-07-04
================

Added
-----
- Amplicon workflow
- Amplicon descriptor schemas
- Amplicon report generator
- Add Rattus norvegicus organism choice to sample schema
- Transforming form Phred 64 to Phred 33 when uploading fastq reads
- Add primertrim process
- RNA-Seq experiment descriptor schema
- iCount sample and reads descriptor schemas
- iCount demultiplexing and sample annotation
- ICount QC
- Add MM8, RN4 and RN6 options to rose2 process
- Add RN4 and RN6 options to bamplot process
- Archive-samples process
- Add bamliquidator
- CheMut workflow
- Dicty primary analysis descriptor schema
- IGV session to Archive-samples process
- Use Resolwe's field projection mixins for knowledge base endpoints
- ``amplicon-table`` process
- Add C. griseus organism choice to Sample descriptor schema
- Add S. tuberosum organism choice to Sample descriptor schema
- Add log2 to gene and sample hierarchical clustering
- Add new inputs to import SRA, add read type selection process
- Set memory resource requirement in jbrowse annotation gff3 and gtf
  processes to 16GB
- Set memory resource requirement in star alignment and index processes
  to 32GB
- Add C. elegans organism choice to Sample descriptor schema
- Add D. melanogaster organism choice to Sample descriptor schema
- Set core resource requirement in Bowtie process to 1
- Set memory resource requirement in amplicon BWA trim process to 32GB
- Add new master file choices to amplicon panel descriptor schema
- Add S. tuberosum organism choice to RNA-seq workflow
- Add Cutadapt process
- Add leaf ordering to gene and sample hierarchical clustering

Fixed
-----
- Use new import paths in ``resolwe.flow``
- Upload reads (paired/single) containing whitespace in the file name
- Fix reads filtering processes for cases where input read file names
  contain whitespace
- Add additional filtering option to STAR aligner
- Fix bbduk-star-htseq_count workflow
- Fix cuffnorm process: Use sample names as labels (boxplot, tables),
  remove group labels input, auto assign group labels, add outputs for
  Rscript output files which were only available compressed
- Derive output filenames in hisat2 from the first reads filename
- Correctly fetch KB features in ``goea.py``
- Append JBrowse tracks to sample
- Replace the BAM MD tag in `align-bwa-trim` process to correct for an
  issue with the primerclip tool
- Fix typo in trimmomatic and bbduk processes
- Use re-import in `etc` and `hmmer_database` processes

Changed
-------
- Support Resolwe test framework
- Run tests in parallel with Tox
- Use Resolwe's new ``FLOW_DOCKER_COMMAND`` setting in test project
- Always run Tox's ``docs``, ``linters`` and ``packaging`` environments
  with Python 3
- Add ``extra`` Tox testing environment with a check that there are no
  large test files in ``resolwe_bio/tests/files``
- Replace Travis CI with Genialis' Jenkins for running the tests
- Store compressed and uncompressed .fasta files in
  ``data:genome:fasta`` objects
- Change sample_geo descriptor schema to have strain option available
  for all organisms
- More readable rna-seq-quantseq schema, field stranded
- Remove obsolete Gene Info processes
- Change log2(fc) default from 2 to 1 in diffexp descriptor schema
- Change Efective genome size values to actual values in macs14 process
- Change variable names in bowtie processes
- Remove iClip processes, tools, files and tests


================
1.3.0 2017-01-28
================

Changed
-------
- Add option to save expression JSON to file before saving it to Storage
- Update ``upload-expression`` process
- No longer treat ``resolwe_bio/tools`` as a Python package
- Move processes' test files to the ``resolwe_bio/tests/files`` directory
  to generalize and simplify handling of tests' files
- Update differential expression (DE) processors
- Update ``generate_diffexpr_cuffdiff`` django-admin command
- Save gene_id source to ``output.source`` for DE, expression and related objects
- Refactor ``upload-diffexp`` processor
- Update sample descriptor schema
- Remove obsolete descriptor schemas
- Add stitch parameter to rose2 processor
- Add filtering to DESeq2
- Set Docker Compose's project name to ``resolwebio`` to avoid name clashes
- GO enrichment analysis: map features using gene Knowledge base
- Add option to upload .gff v2 files with upload-gtf processor
- Replace Haystack with Resolwe Elastic Search API
- Require Resolwe 1.4.1+
- Update processes to be compatible with Resolwe 1.4.0

Added
-----
- Process definition documentation style and text improvements
- Add ``resolwe_bio.kb`` app, Resolwe Bioinformatics Knowledge Base
- Add tests to ensure generators produce the same results
- Upload Gene sets (``data:geneset``)
- Add ``generate_geneset`` django-admin command
- Add ``generate_diffexpr_deseq`` django-admin command
- Add 'Generate GO gene sets' processor
- Add generic file upload processors
- Add upload processor for common image file types (.jpg/.tiff/.png/.gif)
- Add upload processor for tabular file formats (.tab/.tsv/.csv/.txt/.xls/.xlsx)
- Add Trimmomatic process
- Add featureCounts process
- Add Subread process
- Add process for hierarchical clusteing of samples
- Add gff3 to gtf file converter
- Add microarray data descriptor schema
- Add process for differential expression edgeR
- ``BioCollectionFilter`` and ``BidDataFilter`` to support filtering collections
  and data by samples on API
- Added processes for automatically downloading single and paired end SRA files
  from NCBI and converting them to FASTQ
- Added process for automatically downloading SRA files from NCBI and converting
  them to FASTQ
- Add HEAT-Seq pipeline tools
- Add HEAT-Seq workflow
- Add ``create-geneset``, ``create-geneset-venn``  processors
- Add ``source`` filter to feature search endpoint
- Add bamplot process
- Add gene hiererhical clustering
- Add cuffquant workflow
- Support Django 1.10 and versionfield 0.5.0
- django-admin commands ``insert_features`` and ``insert_mappings`` for
  importing features and mappings to the Knowledge Base
- Add bsmap and mcall to analyse WGBS data
- Vaccinesurvey sample descriptor schema
- Add RNA-Seq single and paired-end workflow

Fixed
-----
- Set ``presample`` to ``False`` for Samples created on Sample endpoint
- Fix FastQC report paths in processors
- Fix ``htseq_count`` and ``featureCounts`` for large files
- Fix ``upload gtf annotation``
- Fix gene_id field type for differential expression storage objects
- Order data objects in ``SampleViewSet``
- Fix sample hiererhical clustering
- Fix name in gff to gtf process
- Fix clustering to read expressed genes as strings
- Fix protocol labels in ``rna-seq-quantseq`` descriptor schema


================
1.2.1 2016-07-27
================

Changed
-------
- Update ``resolwe`` requirement


================
1.2.0 2016-07-27
================

Changed
-------
- Decorate all tests that currently fail on Docker with ``skipDockerFailure``
- Require Resolwe's ``master`` git branch
- Put packaging tests in a separate Tox testing environment
- Rename DB user in test project
- Change PostgreSQL port in test project
- Add ROSE2 results parser
- Compute index for HISAT2 aligner on genome upload
- Updated Cuffquant/Cuffnorm tools
- Change ROSE2 enhancer rank plot labels
- Refactor processor syntax
- Move processes tests into ``processes`` subdirectory
- Split ``sample`` API endpoint to ``sample`` for annotated ``Samples``
  and ``presample`` for unannotated ``Samples``
- Rename test project's data and upload directories to ``.test_data`` and
  ``.test_upload``
- Save fastq files to ``lists:basic:file`` field. Refactor related processors.
- Reference genome-index path when running aligners.
- Add pre-computed genome-index files when uploading reference fasta file.
- Include all necessary files for running the tests in source distribution
- Exclude tests from built/installed version of the package
- Move testing utilities from ``resolwe_bio.tests.processes.utils`` to
  ``resolwe_bio.utils.test``
- Update Cuffdiff processor inputs and results table parsing
- Refactor processes to use the updated ``resolwe.flow.executors.run`` command
- Refactor STAR aligner - export expressions as separate objects

Fixed
-----
- Make Tox configuration more robust to different developer environments
- Set ``required: false`` in processor input/output fields where necessary
- Add ``Sample``'s ``Data objects`` to ``Collection`` when ``Sample`` is added
- Fixed/renamed Cufflinks processor field names

Added
-----
- ``skipDockerFailure`` test decorator
- Expand documentation on running tests
- Use Travis CI to run the tests
- Add ``Sample`` model and corresponding viewset and filter
- Add docker-compose command for PostgreSQL
- API endpoint for adding ``Samples`` to ``Collections``
- HISAT2 aligner
- Use Git Large File Storage (LFS) for large test files
- Test for ``generate_samples`` django-admin command
- django-admin command: ``generate_diffexpr``


================
1.1.0 2016-04-18
================

Changed
-------
- Remove obsolete utilities superseded by resolwe-runtime-utils
- Require Resolwe 1.1.0

Fixed
-----
- Update sample descriptor schema
- Include all source files and supplementary package data in sdist

Added
-----
- ``flow_collection: sample`` to processes
- MACS14 processor
- Initial Tox configuration for running the tests
- Tox tests for ensuring high-quality Python packaging
- ROSE2 processor
- django-admin command: ``generate_samples``


================
1.0.0 2016-03-31
================

Changed
-------
- Renamed assertFileExist to assertFileExists
- Restructured processes folder hierarchy
- Removed re-require and hard-coded tools' paths

Fixed
-----
- Different line endings are correctly handled when opening gzipped files
- Fail gracefully if the field does not exist in assertFileExists
- Enabled processor tests (GO, Expression, Variant Calling)
- Enabled processor test (Upload reads with old Illumina QC encoding)
- Made Resolwe Bioinformatics work with Resolwe and Docker

Added
-----
- Import expressions from tranSMART
- Limma differential expression (tranSMART)
- VC filtering tool (Chemical mutagenesis)
- Additional analysis options to Abyss assembler
- API endpoint for Sample
- Initial documentation
