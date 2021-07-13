##########
Change Log
##########

All notable changes to this project are documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


===================
38.2.0 - 2021-07-13
===================

Added
-----
- Add more information about output to the ``methylation-array-sesame``
  pipeline documentation
- Support filtering by ``subject_information.sample_label``,
  ``subject_information.subject_id``, ``subject_information.batch``,
  ``subject_information.group``, ``disease_information.disease_type``,
  ``disease_information.disease_status``,
  ``immuno_oncology_treatment_type.io_drug``,
  ``immuno_oncology_treatment_type.io_treatment``,
  ``response_and_survival_analysis.confirmed_bor``,
  ``response_and_survival_analysis.pfs_event``, ``general.description``,
  ``general.biosample_source``, and ``general.biosample_treatment``
  fields in sample descriptor on API

Changed
-------
- Improve automatic sample naming in the ``geo-import`` process

Fixed
-----
- Fix stalled sam-to-bam conversion in ``wgs-preprocess`` process
- Return column betas to ``methylation-array-sesame`` pipeline output


===================
38.1.1 - 2021-06-14
===================

Changed
-------
- Remove mapping of probe_ids to ENSEMBL ids and add extra variables in
  ``methylation-array-sesame`` process


===================
38.1.0 - 2021-06-14
===================

Added
-----
- Add ``wgs-preprocess`` process
- Add ``gatk-haplotypecaller-gvcf`` process
- Add ``workflow-wgs-gvcf`` process
- Add ``gatk-genotype-gvcfs`` process
- Add ``gatk-vqsr`` process
- Add ``bamtofastq-paired`` process
- Add ``methylation_array`` docker image
- Add ``methylation-array-sesame`` process
- Add support for Python 3.9
- Support downloading knowledge base features and mappings from S3 bucket
- Cap process memory consumption at 10GB

Changed
-------
- Bump GATK to version 4.2.0.0 in ``resolwebio/dnaseq:6.0.0`` Docker
  image
- Update ``workflow-mirna``
- Add new parameters -maximumlength/-M and -no-indels in processes
  ``cutadapt-single`` and ``cutadatp-paired``
- Add new ``id_attribute`` to ``feature_counts`` process

Fixed
-----
- Remove some duplicated code in ``test_probe_mapping``
- Rename FastQC output bundle in Trimmomatic processes so that the
  reports are correctly sorted/included in MultiQC reports
- Fix method signature for KB feature/mapping filtering


===================
38.0.0 - 2021-05-17
===================

Added
-----
- Add bioservices python package to the ``resolwebio/common:2.8.0``
  Docker image
- Add ``upload-idat`` process
- Add ``upload-microarray-expression`` and
  ``mapped-microarray-expression`` processes
- Add ``map-microarray-probes`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Support microarray expressions upload in
  ``geo-import`` process
- Trigger an error for microarray data in differential expression
  processes ``differentialexpression-edger`` and
  ``differentialexpression-deseq2``


===================
37.0.0 - 2021-04-19
===================

Added
-----
- Add GEOparse to the ``resolwebio/common:2.7.0`` Docker image
- Add fastq file validation in ``import-sra-single`` and
  ``import-sra-paired`` processes
- Add ``geo-import`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 28.x
- Use ``resolwebio/base:ubuntu-20.04`` Docker image for building
  ``resolwebio/sra-tools`` Docker image. Include ``dnaio`` Python
  library in ``resolwebio/sra-tools``.

Fixed
-----
- Fix handling of non-sample data inputs in ``multiqc`` process


===================
36.1.0 - 2021-03-15
===================

Added
-----
- Fail if wrong filtering arguments are used in KB Feature / Mapping
  search endpoints

Changed
-------
- Use Amazon ECR when building ``resolwebio/base`` Docker images
- Use pinned version of the ``resolwebio/base`` Docker image for
  building ``resolwebio/common`` Docker image. Update versions of
  bioinformatic tools installed in the ``resolwebio/common`` image.
- Use only tagged versions of ``resolwebio/base`` Docker images in
  processes
- Save gene-level estimated counts to the ``rc`` output field in the
  ``salmon-quant`` process

Fixed
-----
- Fix file import in processes ``upload-multiplexed-single`` and
  ``upload-multiplexed-paired``
- Fix ``import-sra-single`` and ``import-sra-paired`` to correctly
  determine Illumina 1.5 and 1.3 quality encoding


===================
36.0.0 - 2021-02-22
===================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 27.x
- Move docker images from Docker Hub to Amazon ECR


===================
35.0.0 - 2021-01-20
===================

Added
-----
- Add OncXerna specific clinical descriptor schema ``oncxerna_clinical``

Changed
-------
- **BACKWARD INCOMPATIBLE:** Support new protocol in Resolwe 26.x


===================
34.3.0 - 2020-12-14
===================

Added
-------
- Add initial general clinical descriptor schema ``general_clinical``
- Add ``id`` field to ``Feature`` and ``Mapping`` serializers
- Add ``resolwebio/base:ubuntu-20.04`` Docker image

Changed
-------
- Update the url for the Orange table example template in
  ``upload-orange-metadata``


===================
34.2.1 - 2020-11-17
===================

Fixed
-------
- Fix ``macs2-callpeak`` process version


===================
34.2.0 - 2020-11-13
===================

Added
-------
- Add ``upload-proteomics-sample`` and ``upload-proteomics-sample-set``
  processes for uploading custom tables holding proteomics data

Fixed
-------
- Changed ``scale-bigwig`` output file field label to ``bigwig file``
- Bump memory requirements in processes ``import-sra``,
  ``import-sra-single`` and ``import-sra-paired`` to 8GB


===================
34.1.0 - 2020-10-20
===================

Added
-------
- Add peakcalling to removed duplicates step in species' line of the
  ``workflow-cutnrun`` workflow

Fixed
-------
- Add BigWig timeout and bin size parameters to ``markduplicates``,
  ``alignmentsieve`` and ``workflow-cutnrun``. Add bin size parameter
  to ``alignment-bowtie2``.


===================
34.0.0 - 2020-10-19
===================

Added
-------
- Added parameters ``--normalizeUsing`` and ``--smoothLength`` to
  script ``bamtobigwig.sh`` to be used in ``bamCoverage`` program
- Added parameters ``--no-unal`` and ``--no-overlap`` to process
  ``alignment-bowtie``
- Add ``alignmentsieve`` process
- Add Trim Galore tool to ``resolwebio/rnaseq:4.12.0``
- Add ``trimgalore-paired`` process
- Add ``bedtools-bamtobed`` and ``scale-bigwig`` processes
- Added BigWig timeout input parameter to ``alignment-bowtie2`` process
- Add workflow ``workflow-cutnrun``
- Add ``clustering-hierarchical-etc`` process
- Add ``find-similar`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 25.x
- **BACKWARD INCOMPATIBLE:** Rewrite ``differentialexpression-deseq2``
  to Python
- Add format parameter to ``macs2-callpeak``
- Rewrite ``differentialexpression-edger`` to Python
- Rewrite ``cuffdiff`` to Python
- Alignment processes ``alignment-bowtie``, ``alignment-bowtie2``,
  ``alignment-star``, ``alignment-bwa-mem``, ``alignment-bwa-sw``,
  ``alignment-bwa-aln``, ``alignment-hisat2`` and ``walt`` now issue a
  warning instead of an error when sample and genome species mismatch
- Support automated upload of gene sets in proceses ``cuffdiff``,
  ``differentialexpression-deseq2`` and ``differentialexpression-edger``
- Support the analysis of S. cerevisiea samples in ``macs2-callpeak``
  process


===================
33.0.0 - 2020-09-14
===================

Added
-------
- Add ``resolwebio/sra-tools`` Docker image
- Add ``resolwebio/orange`` Docker image
- Add ``upload-orange-metadata`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 24.x
- **BACKWARD INCOMPATIBLE:** Include feature full names in full-text
  search
- Support automatic species annotation in alignment processes:
  ``alignment-bowtie``, ``alignment-bowtie2``, ``alignment-bwa-mem``,
  ``alignment-bwa-sw``, ``alignment-bwa-aln``, ``alignment-hisat2``,
  ``alignment-star``, ``walt``
- Pin ``XML`` R package to ensure compatibility with R 3.6.3 in
  ``resolwebio/chipseq:4.1.3`` Docker image
- Use ``resolwebio/sra-tools:1.0.0`` Docker image in processes
  ``import-sra``, ``import-sra-single`` and ``import-sra-paired``
- Optionally use sra-tools ``prefetch`` command when downloading and
  converting SRA files to FASTQ format

Fixed
-----
- Bump Docker image version in ``chipqc`` process to fix enrichment
  heatmap plot


===================
32.0.0 - 2020-08-17
===================

Added
-------
- Prepare ``resolwebio/rnaseq:4.11.0`` Docker image:
  Add rnanorm (1.3.0) RNA-seq normalization package. Use
  ``resolwebio/common:1.6.0`` Docker image as a base image. Pin ``XML``
  R package to fix the image build issues. Install BBMap package from
  Google Drive.

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 23.x.
- **BACKWARD INCOMPATIBLE:** Use rnanorm Python package for TPM/CPM
  normalization of RNA-seq data in featureCounts and HTSeq-count tools
- Support Nanostring sample reports in MultiQC
- Support Nanostring analysis results in
  ``differentialexpression-deseq2`` process

Fixed
-----
- Order results on autocomplete API endpoint in knowledge-base by
  relevance
- Support filtering by type on knowledge base Feature API
- Attach ``rose2`` Data object to the input sample


===================
31.0.0 - 2020-07-10
===================

Added
-------
- Add Sample QC information fields to the ``sample`` descriptor schema

Changed
-------
- **BACKWARD INCOMPATIBLE:** Disable editing capabilities of Knowledge
  Base API endpoints
- Bump Samtools to version 1.10 in ``resolwebio/common:1.6.0`` Docker
  image
- Migrate search for Knowledge Base enpoints from Elasticsearch to
  PostgreSQL
- Use ``resolwebio/common:1.6.0`` for the ``resolwebio/wgbs:1.3.0``
  Docker image
- Support samtools markdup report in ``walt`` process when removing
  duplicates
- Support samtools markdup report from ``walt`` in MultiQC
- Support samtools markdup report in ``workflow-wgbs-single`` and in
  ``workflow-wgbs-paired`` workflows
- Bump memory requirements to 32GB in processes: ``feature_counts``,
  ``coveragebed``, ``library-strandedness``, ``qorts-qc``,
  ``salmon-quant`` and ``vc-realign-recalibrate``
- Rename ``workflow-slamdunk-paired`` process

Fixed
-------
- Fix read length estimation in ``chipqc``


===================
30.0.0 - 2020-06-15
===================

Added
-----
- Add ``workflow-subsample-bwa-aln-single`` and
  ``workflow-subsample-bwa-aln-paired`` workflows

Changed
-------
- **BACKWARD INCOMPATIBLE:** Use Salmon 1.2.1 in ``salmon-quant`` and
  ``salmon-index`` processes
- Salmon quant 1.2.1 is not backwards compatible with indices generated
  with Salmon index prior to version 1.0.0, thus Salmon tool is updated
  to version 1.2.1 in processes that utilize Salmon to detect library
  strandedness type.
- Expose additional limit options in ``alignment-star`` process
- Bump SRA toolkit to 2.10.0 in ``resolwebio/common:1.5.0`` Docker image
- Use SRA tookit 2.10.0 in ``import-sra``, ``import-sra-single`` and
  ``import-sra-paired`` processes
- Format floats to 2 decimal places in custom ChIP-seq pre/post-peak
  MultiQC reports


===================
29.0.0 - 2020-05-18
===================

Added
-----
- Add filtered BAM output to ``macs2-callpeak`` process
- Add an option to use filtered BAM files from ``macs2-callpeak`` to
  ``rose2``, ``workflow-macs-rose``, and ``macs2-rose2-batch``
- Add ChIPQC to the ``resolwebio/chipseq:4.1.0`` Docker image
- Add ``chipqc`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 22.x
- **BACKWARD INCOMPATIBLE:** Remove processes ``alignment-subread`` and
  ``subread-index``
- **BACKWARD INCOMPATIBLE:** Remove process ``upload-genome``. Refactor
  processes and workflows that required ``data:genome:fasta`` type of
  object on the input to work with ``data:seq:nucleotide`` or dedicated
  aligner index files instead.
- Change ``macs2-batch`` and ``macs2-rose2-batch`` to use tagAlign
  files by default
- Bump Salmon to version 1.2.1 in ``resolwebio/rnaseq:4.10.0`` Docker
  image. Fix build issues affecting ``jpeg`` and ``png`` R packages.
- Support ``chipqc`` process outputs in MultiQC
- Support ``chipqc`` in ``workflow-macs-rose``, ``workflow-macs2``,
  ``macs2-batch`` and ``macs2-rose2-batch`` processes
- Bump memory requirements for process ``upload-fasta-nucl`` to 8 GB

Fixed
-------
- Fix Data name in ``bowtie-index``, ``bowtie2-index``, ``bwa-index``,
  ``hisat2-index`` and ``walt-index``
- Fix filtering of empty VCF files in ``lofreq`` process


===================
28.0.0 - 2020-04-10
===================

Added
-----
- Add ``workflow-wgs-paired`` workflow
- Add processes: ``bowtie-index``, ``bowtie2-index``, ``bwa-index``,
  ``hisat2-index``, ``subread-index`` and ``walt-index``.
- Add ``Dictyostelium purpureum`` species choice to ``sample``
  descriptor schema

Changed
-------
- **BACKWARD INCOMPATIBLE:** Refactor ``upload-fasta-nucl`` process:
  ``species`` and ``build`` input information on FASTA file upload are
  now mandatory, while ``source`` input has been removed.
- **BACKWARD INCOMPATIBLE:** Change the ``alignment-star-index`` process
  type to ``data:index:star``. The process now accepts only
  ``upload-fasta-nucl`` objects on input.
- Add trimming with Trimmomatic in ``workflow-wgbs-single`` and
  ``workflow-wgbs-paired`` workflows
- Make intervals an optional input in ``bqsr`` process
- Make intervals an optional input in ``vc-gatk4-hc`` process
- Bump memory requirements in ``walt`` process to 32 GB

Fixed
-------
- Fix data type of adapters input field in ``alignment-summary`` process
- Fix handling of multiple adapters in ``alignment-summary`` process


===================
27.0.0 - 2020-03-13
===================

Added
-----
- Add ``merge-fastq-single`` and ``merge-fastq-paired`` processes that
  merge multiple ``data:reads:fastq`` data objects into a single
  ``data:reads:fastq`` data object (and consequently a single sample)
- Add ``bs-conversion-rate`` process
- Add support for Python 3.8

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 21.x
- **BACKWARD INCOMPATIBLE:** Split ``workflow-wgbs`` into
  ``workflow-wgbs-single`` and ``workflow-wgbs-paired`` workflows
- Extend the ``workflow-wgbs-single`` and ``workflow-wgbs-paired`` with
  the ``markduplicates``, ``insert-size`` and ``bs-conversion-rate``
  QC processes
- Support detection and separation of control spike-in-derived reads
  from endogenous sequencing reads in ``walt`` process
- Replace duplicate-remover in ``walt`` to unify both (.mr and .bam)
  output alignment files
- Support ``markduplicates`` and ``bs-conversion-rate`` process outputs
  in ``multiqc`` reports
- Enable multiple SRR numbers as inputs in processes ``import-sra``,
  ``import-sra-single``, and ``import-sra-paired``
- Bump memory requirements in ``rrbs-metrics`` process
- Improve process test input data for the ``alignment-star`` process
- Bump Bedtools to v2.29.2 in ``resolwebio/common:1.3.2`` Docker image

Fixed
-----
- Fix Jbrowse track creation in ``upload-genome`` process. When
  gzip input was used in ``prepare-refseqs.pl``, not all sequence chunks
  were created for some inputs.
- Fix ``macs2-callpeak`` process to work with paired-end reads when
  not using tagAlign files
- Fix ``bed_file_corrections_genome_browsers.py`` script to handle cases
  where the input file is empty


===================
26.0.0 - 2020-02-14
===================

Added
-----
- Add ``alignment-summary`` process
- Add ``insert-size`` process
- Add ``wgs-metrics`` process
- Add ``rrbs-metrics`` process
- Add ``workflow-macs2`` workflow

Changed
-------
- **BACKWARD INCOMPATIBLE:** Use featureCounts instead of Stringtie in
  the  ``workflow-corall-single`` and ``workflow-corall-paired``
  workflows
- **BACKWARD INCOMPATIBLE:** Remove ``stringtie`` and
  ``upload-metabolic-pathway`` processes
- **BACKWARD INCOMPATIBLE:** Refactor ``walt`` process to support
  Picard quality metrics and update ``methcounts`` process and to match
  the new outputs
- **BACKWARD INCOMPATIBLE:** Support MultiQC report in ``wgbs`` workflow
- Remove Stringtie tool from ``resolwebio/rnaseq`` Docker image
- Remove ``resolwe/base:ubuntu-14.04`` and ``resolwe/base:ubuntu-17.10``
  Docker images
- Use pigz for output file compression in ``bbduk-single`` and
  ``bbduk-paired`` processes
- Use ``resolwebio/rnaseq:4.9.0`` Docker image in processes
  ``bbduk-single``, ``bbduk-paired``, ``trimmomatic-single``,
  ``trimmomatic-paired``, ``alignment-bowtie``, ``alignment-bowtie2``,
  ``alignment-hisat2``, ``alignment-subread``, ``cuffmerge``, ``pca``,
  ``cuffdiff``, ``differentialexpression-edger``, ``cufflinks``,
  ``cuffnorm``, ``cuffquant``, ``expression-aggregator``,
  ``htseq-count``, ``htseq-count-raw``, ``index-fasta-nucl``, ``rsem``,
  ``upload-bam``, ``upload-bam-indexed``, ``upload-bam-secondary``,
  ``upload-expression``, ``upload-expression-cuffnorm``,
  ``upload-expression-star``, ``upload-genome``,
  ``upload-gaf``, ``upload-obo``, ``upload-fasta-nucl``,
  ``regtools-junctions-annotate``, ``cutadapt-custom-single``,
  ``cutadapt-custom-paired``, ``bam-split``, ``gff-to-gtf``,
  ``spikein-qc``, ``differentialexpression-shrna``, ``feature_counts``,
  ``salmon-index``, ``salmon-quant``, ``library-strandedness``,
  ``qorts-qc``, ``alignment-star``, ``alignment-star-index``,
  ``cutadapt-3prime-single``, ``cutadapt-single``, ``cutadapt-paired``,
  ``differentialexpression-deseq2``, ``cutadapt-corall-single``,
  ``cutadapt-corall-paired``, ``umi-tools-dedup`` and ``shrna-quant``.
- Use ``resolwebio/common:1.3.1`` Docker image in processes
  ``amplicon-table``, ``mergeexpressions``, ``upload-bedpe``,
  ``upload-bam-scseq-indexed``, ``upload-diffexp``, ``upload-etc``,
  ``upload-sc-10x``, ``upload-multiplexed-single``,
  ``upload-multiplexed-paired``, ``archive-samples``,
  ``samtools-idxstats``, ``seqtk-sample-single``,
  ``seqtk-sample-paired``, ``basespace-file-import``,
  ``clustering-hierarchical-samples``,
  ``clustering-hierarchical-genes``, ``import-sra``,
  ``import-sra-single``, ``import-sra-paired``.
- Compute TPM values and map gene_ids to gene symbols in
  ``alleyoop-collapse`` process output
- Rewrite ``multiqc`` process to Python
- Save ``lib_format_counts.json`` in a separate output field in the
  ``salmon-quant`` process
- Use ``resolwebio/common:1.3.1`` as a base Docker image for the
  ``resolwebio/wgbs:1.2.0`` Docker image
- Support MultiQC reports in ChIP-seq workflows

Fixed
-----
- Fix Mapping search for ``source_id`` / ``target_id``
- Fix handling of input file names in processes: ``cellranger-count``,
  ``cutadapt-3prime-single``, ``cutadapt-corall-single``,
  ``cutadapt-corall-paired``, ``salmon-quant``, ``umi-tools-dedup``,
  ``upload-sc-10x`` and ``upload-bam-scseq-indexed``
- Fix handling of chimeric alignments in ``alignment-star``


===================
25.1.0 - 2020-01-14
===================

Added
-----

Changed
-------
- Extend the MultiQC report so that the Sample summary table is created
  for the compatible Data objects
- Bump CPU and memory requirements for the ``alignment-bowtie2`` process
- Move upload test files of differential expression to its own folder

Fixed
-----
- Fix typo in ``scheduling_class`` variable in several Python processes
- Handle cases of improper tags passed to ``read_group`` argument of
  the ``bqsr`` process
- When processing differential expression files, a validation is
  performed for numeric columns


===================
25.0.0 - 2019-12-17
===================

Added
-----
- Add ``alleyoop-rates`` process
- Add ``alleyoop-utr-rates`` process
- Add ``alleyoop-summary`` process
- Add ``alleyoop-snpeval`` process
- Add ``alleyoop-collapse`` process
- Add ``slam-count`` process
- Add ``workflow-slamdunk-paired`` workflow

Changed
-------
- **BACKWARD INCOMPATIBLE:** Refactor ``slamdunk-all-paired`` process
  to support genome browser visualization and add additional output
  fields
- Append sample and genome reference information to the summary output
  file in the ``filtering-chemut`` process
- Bigwig output field in ``bamclipper``, ``bqsr`` and ``markduplicates``
  processes is no longer required
- Support Slamdunk/Alleyoop processes in MultiQC
- Enable sorting of files in ``alignment-star`` process using Samtools
- Support merging of multi-lane sequencing data into a single (pair) of
  FASTQ files in the ``upload-fastq-single``, ``upload-fastq-paired``,
  ``files-to-fastq-single`` and ``files-to-fastq-paired`` processes


===================
24.0.0 - 2019-11-15
===================

Added
-----
- Add ``resolwebio/slamdunk`` Docker image
- Add Tabix (1.7-2) to ``resolwebio/bamliquidator:1.2.0`` Docker image
- Add ``seqtk-rev-complement-single`` and
  ``seqtk-rev-complement-paired`` process
- Add ``slamdunk-all-paired`` process

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 20.x
- Make BaseSpace file download more robust
- Bump ``rose2`` to 1.1.0, ``bamliquidator`` to 1.3.8, and use
  ``resolwebio/base:ubuntu-18.04`` Docker image as a base image in
  ``resolwebio/bamliquidator:1.1.0`` Docker image
- Use ``resolwebio/bamliquidator:1.2.0`` in ``rose2`` process
- Bump CPU, memory and Docker image (``resolwebio/rnaseq:4.9.0``)
  requirements in ``alignment-bwa-mem``, ``alignment-bwa-sw`` and
  ``alignment-bwa-aln`` processes
- Use multi-threading option in Samtools commands in
  ``alignment-bwa-mem``, ``alignment-bwa-sw`` and ``alignment-bwa-aln``
  processes


===================
23.1.1 - 2019-10-11
===================

Changed
-------
- Renamed ``workflow-trim-align-quant`` workflow to make the name more
  informative


===================
23.1.0 - 2019-09-30
===================

Added
-----
- Add ``Macaca mulatta`` species choice to the ``sample`` descriptor
  schema
- Add ``workflow-cutadapt-star-fc-quant-wo-depletion-single`` process

Changed
-------
- Test files improved for ``workflow-wes``, ``bamclipper``,
  ``markduplicates`` and ``bqsr``
- Fix typo in ``differentialexpression-shrna`` process docstring

Fixed
-----
- Fix transcript-to-gene_id mapping for Salmon expressions in
  ``differentialexpression-deseq2`` process. Transcript versions are
  now ignored when matching IDs using the transcript-to-gene_id mapping
  table.
- Fix ``workflow-cutadapt-star-fc-quant-single`` process description


===================
23.0.0 - 2019-09-17
===================

Changed
-------
- Update order of QC reports in MultiQC configuration file. The updated
  configuration file is part of the ``resolwebio/common:1.3.1``
  Docker image.
- Bump Jbrowse to version 1.16.6 in ``resolwebio/rnaseq:4.9.0`` Docker
  image
- Use JBrowse ``generate-names.pl`` script to index GTF/GFF3 features
  upon annotation file upload
- Support Salmon reports in MultiQC and expose ``dirs_depth`` parameter
- Expose transcript-level expression file in the ``salmon-quant``
  process

Added
-----
- Add ``workflow-bbduk-salmon-qc-single`` and
  ``workflow-bbduk-salmon-qc-paired`` workflows

Fixed
-----
- Give process ``upload-bedpe`` access to network


===================
22.0.0 - 2019-08-20
===================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 19.x
- **BACKWARD INCOMPATIBLE:** Unify ``cutadapt-single`` and
  ``cutadapt-paired`` process inputs and refactor to use Cutadapt v2.4
- Expose BetaPrior parameter in ``differentialexpression-deseq2``
  process
- Install R from CRAN-maintained repositories in Docker images build
  from the ``resolwebio/base:ubuntu-18.04`` base image
- Prepare ``resolwebio/common:1.3.0`` Docker image:

  - Install R v3.6.1
  - Bump Resdk to v10.1.0
  - Install gawk package
  - Fix Docker image build issues
- Use ``resolwebio/common:1.3.0`` as a base image for
  ``resolwebio/rnaseq:4.8.0``
- Update StringTie to v2.0.0 in ``resolwebio/rnaseq:4.8.0``
- Support StringTie analysis results in DESeq2 tool

Added
-----
- Add ``cutadapt-3prime-single`` process
- Add ``workflow-cutadapt-star-fc-quant-single`` process
- Add argument ``skip`` to ``bamclipper`` which enables skipping of
  the said process
- Add ``cutadapt-corall-single`` and ``cutadapt-corall-paired``
  processes for pre-processing of reads obtained using Corall Total
  RNA-seq library prep kit
- Add ``umi-tools-dedup`` process
- Add ``stringtie`` process
- Add ``workflow-corall-single`` and ``workflow-corall-paired``
  workflows optimized for Corall Total RNA-seq library prep kit data

Fixed
-----
- Fix warning message in hierarchical clustering of genes. Incorrect
  gene names were reported in the warning message about removed
  genes. Computation of hierarchical clustering was correct.


===================
21.0.1 - 2019-07-26
===================

Changed
-------
- Bump Cutadapt to v2.4 and use ``resolwebio/common:1.2.0`` as a base
  image in ``resolwebio/rnaseq:4.6.0``

Added
-----
- Add pigz package to ``resolwebio/common:1.2.0`` Docker image
- Add StringTie and UMI-tools to ``resolwebio/rnaseq:4.7.0`` Docker
  image

Fixed
-----
- Fix ``spikeins-qc`` process to correctly handle the case where all
  expressions are without spikeins
- Fix an error in ``macs2-callpeak`` process that prevented correct
  reporting of build/species mismatch between inputs
- Support UCSC annotations in ``feature_counts`` process by assigning
  empty string gene_ids to the "unknown" gene


===================
21.0.0 - 2019-07-16
===================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 18.x
- Bump the number of allocated CPU cores to 20 in ``alignment-bwa-mem``
  process
- Bump memory requirements in ``seqtk-sample-single`` and
  ``seqtk-sample-paired`` processes
- Bump Salmon to v0.14.0 in ``resolwebio/rnaseq:4.5.0`` Docker image
- Expose additional inputs in ``salmon-index`` process
- Use ``resolwebio/rnaseq:4.5.0`` Docker image in processes that call
  Salmon tool (``library-strandedness``, ``feature_counts`` and
  ``qorts-qc``)
- Implement dropdown menu for ``upload-bedpe`` process
- Add validation stringency parameter to ``bqsr`` process and propagate
  it to the ``workflow-wes`` as well
- Add LENIENT value to validation stringency parameter of the
  ``markduplicates`` process
- Improve performance of RPKUM normalization in ``featureCounts`` process

Added
-----
- Add ``salmon-quant`` process

Fixed
-----
- Fix genome upload process to correctly handle filenames with dots
- Fix merging of expressions in ``archive-samples`` process. Previously
  some genes were missing in the merged expression files. The genes that
  were present had expression values correctly assigned. The process was
  optimized for performance and now supports parallelization.


=================
20.0.0 2019-06-19
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 17.x
- **BACKWARD INCOMPATIBLE:** Use Elasticsearch version 6.x
- **BACKWARD INCOMPATIBLE:** Bump Django requirement to version 2.2
- **BACKWARD INCOMPATIBLE:** Remove obsolete RNA-seq workflows
  ``workflow-bbduk-star-featurecounts-single``,
  ``workflow-bbduk-star-featurecounts-paired``,
  ``workflow-cutadapt-star-featurecounts-single`` and
  ``workflow-cutadapt-star-featurecounts-paired``
- **BACKWARD INCOMPATIBLE:** Remove obsolete descriptor schemas:
  ``rna-seq-bbduk-star-featurecounts``, ``quantseq``,
  ``rna-seq-cutadapt-star-featurecounts`` and
  ``kapa-rna-seq-bbduk-star-featurecounts``
- **BACKWARD INCOMPATIBLE:** In ``upload-fasta-nucl`` process, store
  compressed and uncompressed FASTA files in ``fastagz`` and ``fasta``
  ouput fields, respectively
- Allow setting the Java memory usage flags for the QoRTs tool in
  ``resolwebio/common:1.1.3`` Docker image
- Use ``resolwebio/common:1.1.3`` Docker image as a base image for
  ``resolwebio/rnaseq:4.4.2``
- Bump GATK4 version to 4.1.2.0 in ``resolwebio/dnaseq:4.2.0``
- Use MultiQC configuration file and prepend directory name to sample
  names by default in ``multiqc`` process
- Bump ``resolwebio/common`` to 1.1.3 in ``resolwebio/dnaseq:4.2.0``
- Process ``vc-gatk4-hc`` now also accepts BED files through parameter
  ``intervals_bed``

Added
-----
- Support Python 3.7
- Add Tabix (1.7-2) to ``resolwebio/wgbs`` docker image
- Add JBrowse index output to ``hmr`` process
- Add ``bamclipper`` tool and ``parallel`` package to ``resolwebio/dnaseq:4.2.0`` image
- Support ``hg19_mm10`` hybrid genome in ``bam-split`` process
- Support mappability-based normalization (RPKUM) in featureCounts
- Add BEDPE upload process
- Add ``bamclipper`` process
- Add ``markduplicates`` process
- Add ``bqsr`` (BaseQualityScoreRecalibrator) process
- Add whole exome sequencing (WES) pipeline

Fixed
-----
- Fix building problems of ``resolwebio/dnaseq`` docker
- Fix handling of no-adapters input in workflows
  ``workflow-bbduk-star-featurecounts-qc-single`` and
  ``workflow-bbduk-star-featurecounts-qc-paired``


=================
19.0.1 2019-05-13
=================

Fixed
-----
- Use ``resolwebio/rnaseq:4.4.2`` Docker image that enforces the memory limit
  and bump memory requirements for ``qorts-qc`` process
- Bump memory requirements for ``multiqc`` process


=================
19.0.0 2019-05-07
=================

Changed
-------
- Use Genialis fork of MultiQC 1.8.0b in ``resolwebio/common:1.1.2``
- Support Samtools idxstats and QoRTs QC reports in ``multiqc`` process
- Support ``samtools-idxstats`` QC step in workflows:

  - ``workflow-bbduk-star-featurecounts-qc-single``
  - ``workflow-bbduk-star-featurecounts-qc-paired``
  - ``workflow-bbduk-star-fc-quant-single``
  - ``workflow-bbduk-star-fc-quant-paired``
- Simplify ``cellranger-count`` outputs folder structure
- Bump STAR aligner to version 2.7.0f in ``resolwebio/rnaseq:4.4.1``
  Docker image
- Use ``resolwebio/rnaseq:4.4.1`` in ``alignment-star`` and
  ``alignment-star-index`` processes
- Save filtered count-matrix output file produced by DESeq2 differential
  expression process

Added
-----
- Add ``samtools-idxstats`` process
- Improve ``cellranger-count`` and ``cellranger-mkref`` logging
- Add FastQC report to ``upload-sc-10x`` process

Fixed
-----
- Fix ``archive-samples`` to work with ``data:chipseq:callpeak:macs2``
  data objects when downloading only peaks without QC reports
- Fix parsing gene set files with empty lines to avoid saving gene sets
  with empty string elements


=================
18.0.0 2019-04-16
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 16.x
- **BACKWARD INCOMPATIBLE:** Rename and improve descriptions of
  processes specific to CATS RNA-seq kits. Remove related
  ``cutadapt-star-htseq`` descriptor schema.
- **BACKWARD INCOMPATIBLE:** Remove ``workflow-accel-gatk4`` pipeline.
  Remove ``amplicon-panel``, ``amplicon-panel-advanced`` and
  ``amplicon-master-file`` descriptor schemas.
- **BACKWARD INCOMPATIBLE:** Remove obsolete processes and descriptor
  schemas: ``rna-seq-quantseq``, ``bcm-workflow-rnaseq``,
  ``bcm-workflow-chipseq``, ``bcm-workflow-wgbs``, ``dicty-align-reads``,
  ``dicty-etc``, ``affy`` and ``workflow-chip-seq``
- Expose additional parameters of ``bowtie2`` process
- Support strandedness auto detection in ``qorts-qc`` process

Added
-----
- Add shRNAde (v1.0) R package to the ``resolwebio/rnaseq:4.4.0`` Docker image
- Add ``resolwebio/scseq`` Docker image
- Add shRNA differential expression process. This is a two-step process which
  trims, aligns and quantifies short hairpin RNA species. These are then used
  in a differential expression.
- Add ``sc-seq`` processes:

  - ``cellranger-mkref``
  - ``cellranger-count``
  - ``upload-sc-10x``
  - ``upload-bam-scseq-indexed``

Fixed
-----
- Bump memory requirements in ``seqtk-sample-single`` and
  ``seqtk-sample-paired`` processes
- Fix ``cellranger-count`` html report
- Mark spliced-alignments with XS flags in ``workflow-rnaseq-cuffquant``
- Fix whitespace handling in ``cuffnorm`` process


=================
17.0.0 2019-03-19
=================

Added
-----
- Add ``qorts-qc`` (Quality of RNA-seq Tool-Set QC) process
- Add ``workflow-bbduk-star-fc-quant-single`` and
  ``workflow-bbduk-star-fc-quant-paired`` processes
- Add independent gene filtering and gene filtering based on Cook's distance
  in ``DESeq2`` differential expression process

Changed
-------
- **BACKWARD INCOMPATIBLE**: Move gene filtering by expression count
  input to ``filter.min_count_sum`` in ``DESeq2`` differential expression
  process
- **BACKWARD INCOMPATIBLE:** Require Resolwe 15.x
- Update ``resolwebio/common:1.1.0`` Docker image:

  - add QoRTs (1.3.0) package
  - bump MultiQC to 1.7.0
  - bump Subread package to 1.6.3
- Expose ``maxns`` input parameter in ``bbduk-single`` and
  ``bbduk-paired`` processes. Make this parameter available in workflows
  ``workflow-bbduk-star-featurecounts-qc-single``,
  ``workflow-bbduk-star-featurecounts-qc-paired``,
  ``workflow-bbduk-star-featurecounts-single`` and
  ``workflow-bbduk-star-featurecounts-paired``.
- Save CPM-normalized expressions in ``feature_counts`` process. Control
  the default expression normalization type (``exp_type``) using the
  ``normalization_type`` input.
- Bump MultiQC to version 1.7.0 in ``multiqc`` process
- Use ``resolwebio/rnaseq:4.3.0`` with Subread/featureCounts version
  1.6.3 in ``feature_counts`` process


=================
16.3.0 2019-02-19
=================

Changed
-------
- Bump STAR aligner version to 2.7.0c in ``resolwebio/rnaseq:4.2.2``
- Processes ``alignment-star`` and ``alignment-star-index`` now use Docker
  image ``resolwebio/rnaseq:4.2.2`` which contains STAR version ``2.7.0c``
- Persistence of ``basespace-file-import`` process changed from ``RAW`` to
  ``TEMP``

Added
-----
- Make ``prepare-geo-chipseq`` work with both
  ``data:chipseq:callpeak:macs2`` and
  ``data:chipseq:callpeak:macs14`` as inputs

Fixed
-----
- Report correct total mapped reads and mapped reads percentage in
  prepeak QC report for ``data:alignment:bam:bowtie2`` inputs in
  ``macs2-callpeak`` process


=================
16.2.0 2019-01-28
=================

Changed
-------
- Enable multithreading mode in ``alignment-bwa-aln`` and
  ``alignment-bwa-sw``
- Lineary lower the timeout for BigWig calculation when running on
  multiple cores

Fixed
-----
- Remove ``pip`` ``--process-dependency-links`` argument in testenv
  settings
- Fix walt getting killed when ``sort`` runs out of memory. The ``sort``
  command buffer size was limited to the process memory limit.


=================
16.1.0 2019-01-17
=================

Changed
-------

Added
-----
- Add the ``FASTQ`` file validator script to the ``upload-fastq-single``,
  ``upload-fastq-paired``, ``files-to-fastq-single`` and
  ``files-to-fastq-paired`` processes
- Add ``spikein-qc`` process
- Add to ``resolwebio/rnaseq:4.1.0`` Docker image:

  - ``dnaio`` Python library
- Add to ``resolwebio/rnaseq:4.2.0`` Docker image:

  - ERCC table
  - common Genialis fonts and css file
  - spike-in QC report template
- Set ``MPLBACKEND`` environment variable to ``Agg`` in
  ``resolwebio/common:1.0.1`` Docker image

Fixed
-----
- Fix the format of the output ``FASTQ`` file in the ``demultiplex.py``
  script
- Fix NSC and RSC QC metric calculation for ATAC-seq and paired-end
  ChIP-seq samples in ``macs2-callpeak`` and ``qc-prepeak`` processes


=================
16.0.0 2018-12-19
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 14.x
- **BACKWARD INCOMPATIBLE:** Remove obsolete processes ``findsimilar``
- **BACKWARD INCOMPATIBLE:** Include ENCODE-proposed QC analysis metrics
  methodology in the ``macs2-callpeak`` process. Simplified MACS2
  analysis inputs now allow the use of sample relations
  (treatment/background) concept to trigger multiple MACS2 jobs
  automatically using the ``macs2-batch`` or ``macs2-rose2-batch``
  processes.
- **BACKWARD INCOMPATIBLE:** Update ``workflow-atac-seq`` inputs to
  match the updated ``macs2-callpeak`` process
- Use ``resolwebio/rnaseq:4.0.0`` Docker image in
  ``alignment-star-index``, ``bbduk-single``, ``bbduk-paired``,
  ``cuffdiff``, ``cufflinks``, ``cuffmerge``, ``cuffnorm``,
  ``cuffquant``, ``cutadapt-custom-single``, ``cutadapt-custom-paired``,
  ``cutadapt-single``, ``cutadapt-paired``,
  ``differentialexpression-deseq2``, ``differentialexpression-edger``,
  ``expression-aggregator``, ``feature_counts``, ``goenrichment``,
  ``htseq-count``, ``htseq-count-raw``, ``index-fasta-nucl``,
  ``library-strandedness``, ``pca``, ``regtools-junctions-annotate``,
  ``rsem``, ``salmon-index``, ``trimmomatic-single``,
  ``trimmomatic-paired``, ``upload-expression``,
  ``upload-expression-cuffnorm``, ``upload-expression-star``,
  ``upload-fasta-nucl``, ``upload-fastq-single``,
  ``upload-fastq-paired``, ``files-to-fastq-single``,
  ``files-to-fastq-paired``, ``upload-gaf``, ``upload-genome``,
  ``upload-gff3``, ``upload-gtf`` and ``upload-obo``
- Order statistical groups in expression aggregator output by sample
  descriptor field value
- Use ``resolwebio/biox:1.0.0`` Docker image in ``etc-bcm``,
  ``expression-dicty`` and ``mappability-bcm`` processes
- Use ``resolwebio/common:1.0.0`` Docker image in ``amplicon-table``,
  ``mergeexpressions``, ``upload-diffexp``, ``upload-etc``,
  ``upload-multiplexed-single`` and ``upload-multiplexed-paired``
  processes
- Use ``resolwebio/base:ubuntu-18.04`` Docker image in
  ``create-geneset``, ``create-geneset-venn``,  ``mergeetc``,
  ``prepare-geo-chipseq``, ``prepare-geo-rnaseq``, ``upload-cxb``,
  ``upload-geneset``, ``upload-header-sam``, ``upload-mappability``,
  ``upload-snpeff`` and ``upload-picard-pcrmetrics`` processes
- Update GATK4 to version 4.0.11.0 in ``resolwebio/dnaseq:4.1.0`` Docker
  image. Install and use JDK v8 by default to ensure compatibility with
  GATK4 package.
- Use ``resolwebio/dnaseq:4.1.0`` Docker image in ``align-bwa-trim``,
  ``coveragebed``, ``filtering-chemut``, ``lofreq``,
  ``picard-pcrmetrics``, ``upload-master-file``, ``upload-variants-vcf``
  and ``vc-gatk4-hc`` processes
- Expose reads quality filtering (q) parameter, reorganize inputs and
  rename the stats output file in ``alignment-bwa-aln`` process
- Use ``resolwebio/chipseq:4.0.0`` Docker image in ``chipseq-genescore``,
  ``chipseq-peakscore``, ``macs14``, ``upload-bed`` and ``qc-prepeak``
  processes
- Use ``resolwebio/bamliquidator:1.0.0`` Docker image in
  ``bamliquidator`` and ``bamplot`` processes

Added
-----
- Add biosample source field to ``sample`` descriptor schema
- Add ``background_pairs`` Jinja expressions filter that accepts a list of
  data objects and orders them in a list of pairs (case, background) based on
  the background relation between corresponding samples
- Add ``chipseq-bwa`` descriptor schema. This schema specifies the
  default inputs for BWA ALN aligner process as defined in ENCODE
  ChIP-Seq experiments.
- Add support for MACS2 result files to MultiQC process
- Add ``macs2-batch``, ``macs2-rose2-batch`` and ``workflow-macs-rose``
  processes
- Add feature symbols to expressions in ``archive-samples`` process

Fixed
-----
- Make ChIP-seq fields in ``sample`` descriptor schema visible when
  ChIPmentation assay type is selected
- Fix handling of whitespace in input BAM file name in script
  ``detect_strandedness.sh``
- Set available memory for STAR aligner to 36GB. Limit the available
  memory for STAR aligner ``--limitBAMsortRAM`` parameter to 90% of the
  Docker requirements setting
- Set ``bbduk-single`` and ``bbduk-paired`` memory requirements to 8GB
- Fix wrong file path in ``archive-samples`` process


=================
15.0.0 2018-11-20
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Remove obsolete processes: ``bsmap``,
  ``mcall``, ``coverage-garvan``, ``igv``, ``jbrowse-bed``,
  ``jbrowse-gff3``, ``jbrowse-gtf``, ``jbrowse-bam-coverage``,
  ``jbrowse-bam-coverage-normalized``, ``jbrowse-refseq``,
  ``fastq-mcf-single``, ``fastq-mcf-paired``, ``hsqutils-trim``,
  ``prinseq-lite-single``, ``prinseq-lite-paired``,
  ``sortmerna-single``, ``sortmerna-paired``, ``bam-coverage``,
  ``hsqutils-dedup``, ``vc-samtools``, ``workflow-heat-seq`` and
  ``alignment-tophat2``
- **BACKWARD INCOMPATIBLE:** Remove ``jbrowse-bam-coverage`` process
  step from the ``workflow-accel`` workflow. The bigwig coverage track
  is computed in ``align-bwa-trim`` process instead.
- **BACKWARD INCOMPATIBLE:** Remove ``resolwebio/utils`` Docker image.
  This image is replaced by the ``resolwebio/common`` image.
- **BACKWARD INCOMPATIBLE:** Use ``resolwebio/common`` Docker image
  as a base image for the ``resolwebio/biox``, ``resolwebio/chipseq``,
  ``resolwebio/dnaseq`` and ``resolwebio/rnaseq`` images
- **BACKWARD INCOMPATIBLE:** Remove ``resolwebio/legacy`` Docker image.
- Use sample name as the name of the data object in:

  - ``alignment-bwa-aln``
  - ``alignment-bowtie2``
  - ``qc-prepeak``
  - ``macs2-callpeak``
- Attach ``macs2-callpeak``, ``macs14`` and ``rose2`` process data to
  the case/treatment sample
- Use ``resolwebio/dnaseq:4.0.0`` docker image in ``align-bwa-trim``
  process
- Use ``resolwebio/rnaseq:4.0.0`` docker image in aligners:
  ``alignment-bowtie``, ``alignment-bowtie2``, ``alignment-bwa-mem``,
  ``alignment-bwa-sw``, ``alignment-bwa-aln``, ``alignment-hisat2``,
  ``alignment-star`` and ``alignment-subread``.
- Set memory limits in ``upload-genome``, ``trimmomatic-single`` and
  ``trimmomatic-paired`` processes
- Improve error messages in differential expression process ``DESeq2``

Added
-----
- Add ``makedb (WALT 1.01)`` - callable as ``makedb-walt``, tool to
  create genome index for WALT aligner, to ``resolwebio/rnaseq`` docker
  image
- Add ``resolwebio/wgbs`` docker image including the following tools:

  - ``MethPipe (3.4.3)``
  - ``WALT (1.01)``
  - ``wigToBigWig (kent-v365)``
- Add ``resolwebio/common`` Docker image. This image includes common
  bioinformatics utilities and can serve as a base image for other,
  specialized ``resolwebio`` Docker images: ``resolwebio/biox``,
  ``resolwebio/chipseq``, ``resolwebio/dnaseq``
  and ``resolwebio/rnaseq``.
- Add ``shift`` (user-defined cross-correlation peak strandshift) input
  to ``qc-prepeak`` process
- Add ATAC-seq workflow
- Compute index for ``WALT`` aligner on genome upload and support
  uploading the index together with the genome
- Add ``Whole genome bisulfite sequencing`` workflow and related WGBS
  processes:

  - ``WALT``
  - ``methcounts``
  - ``HMR``
- Add bedClip to `resolwebio/chipseq:3.1.0` docker image
- Add ``resolwebio/biox`` Docker image. This image is based on the
  ``resolwebio/common`` image and includes Biox Python library for
  Dictyostelium RNA-Seq analysis support.
- Add ``resolwebio/snpeff`` Docker image. The image includes
  SnpEff (4.3K) tool.
- Add spike-in names, rRNA and globin RNA cromosome names in
  ``resolwebio/common`` image
- Add UCSC bedGraphtoBigWig tool for calculating BigWig in
  ``bamtobigwig.sh`` script. In ``align-bwa-trim`` processor set this
  option (that BigWig is calculated by UCSC tool instead of deepTools),
  because it is much faster for amplicon files. In other processors update
  the input parameters for ``bamtobigwig.sh``: ``alignment-bowtie``,
  ``alignment-bowtie2``, ``alignment-bwa-mem``, ``alignment-bwa-sw``,
  ``alignment-bwa-aln``, ``alignment-hisat2``, ``alignment-star``
  ``alignment-subread``, ``upload-bam``, ``upload-bam-indexed`` and
  ``upload-bam-secondary``.
- In ``bamtobigwig.sh`` don't create BigWig when bam file was aligned on
  globin RNA or rRNA (this are QC steps and BigWig is not needed)

Fixed
-----
- **BACKWARD INCOMPATIBLE:** Use user-specificed distance metric in
  hierarchical clustering
- Handle integer expression values in hierarchical clustering
- Fix Amplicon table gene hyperlinks for cases where multiple genes
  are associated with detected variant
- Handle empty gene name in expression files in PCA
- Fix PBC QC reporting  in ``qc-prepeak`` process for a case where
  there are no duplicates in the input bam
- Fix ``macs2-callpeak`` process so that user defined fragment lenth
  has priority over the ``qc-prepeak`` estimated fragment length when
  shifting reads for post-peakcall QC
- Fix ``macs2-callpeak`` to prevent the extension of intervals beyond
  chromosome boundaries in MACS2 bedgraph outputs
- Fix warning message in hierarchical clustering of genes to display gene
  names


=================
14.0.2 2018-10-23
=================

Fixed
-----
- Fix ``htseq-count-raw`` process to correctly map features with associated
  feature symbols.


=================
14.0.1 2018-10-23
=================

Fixed
-----
- Handle missing gene expression in hierarchical clustering of genes. If one or
  more genes requested in gene filter are missing in selected expression files
  a warning is issued and hierarchical clustering of genes is computed with the
  rest of the genes instead of failing.
- Fix PCA computation for single sample case


=================
14.0.0 2018-10-09
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 13.x
- **BACKWARD INCOMPATIBLE:** Remove ``gsize`` input from
  ``macs2-callpeak`` process and automate genome size selection
- **BACKWARD INCOMPATIBLE:** Set a new default ``sample`` and ``reads``
  descriptor schema. Change slug from ``sample2`` to ``sample``, modify group
  names, add ``cell_type`` field to the new ``sample`` descriptor schema, and
  remove the original ``sample``, ``sample-detailed``, and ``reads-detailed``
  descriptor schemas.
- **BACKWARD INCOMPATIBLE:** Unify types of ``macs14`` and
  ``macs2-callpeak`` processes and make ``rose2`` work with both
- **BACKWARD INCOMPATIBLE:** Remove ``replicates`` input in ``cuffnorm``
  process. Use sample relation information instead.
- Use ``resolwebio/chipseq:3.0.0`` docker image in the following processes:

  - ``macs14``
  - ``macs2-callpeak``
  - ``rose2``
- Downgrade primerclip to old version (v171018) in ``resolwebio/dnaseq:3.3.0``
  docker image and move it to google drive.
- Move ``bam-split`` process to ``resolwebio/rnaseq:3.7.1`` docker image
- Count unique and multimmaping reads in ``regtools-junctions-annotate``
  process

Added
-----
- Add ``qc-prepeak`` process that reports ENCODE3 accepted ChIP-seq and
  ATAC-seq QC metrics
- Add QC report to ``macs2-callpeak`` process
- Add combining ChIP-seq QC reports in ``archive-samples`` process
- Add detection of globin-derived reads as an additional QC step in the
  ``workflow-bbduk-star-featurecounts-qc-single`` and
  ``workflow-bbduk-star-featurecounts-qc-paired`` processes.
- Add mappings from ENSEMBL or NCBI to UCSC chromosome names and deepTools
  (v3.1.0) to ``resolwebio/dnaseq:3.3.0`` docker image
- Add BigWig output field to following processors:

  - ``align-bwa-trim``
  - ``upload-bam``
  - ``upload-bam-indexed``
  - ``upload-bam-secondary``
- Add ``replicate_groups`` Jinja expressions filter that accepts a list of
  data objects and returns a list of labels determining replicate groups.
- Add 'Novel splice junctions in BED format' output to
  ``regtools-junctions-annotate`` process, so that user can visualize only
  novel splice juntions in genome browsers.

Fixed
-----
- Fix handling of numerical feature_ids (NCBI source) in
  ``create_expression_set.py`` script
- Make ``chipseq-peakscore`` work with gzipped narrowPeak input from
  ``macs2-callpeak``
- Use uncompressed FASTQ files as input to STAR aligner to prevent
  issues on (network) filesystems without FIFO support


=================
13.0.0 2018-09-18
=================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 12.x
- **BACKWARD INCOMPATIBLE:** Remove obsolete processes: ``assembler-abyss``,
  ``cutadapt-amplicon``, ``feature_location``, ``microarray-affy-qc``,
  ``reads-merge``, ``reference_compatibility``, ``transmart-expressions``,
  ``upload-hmmer-db``, ``upload-mappability-bigwig``,
  ``upload-microarray-affy``.
- **BACKWARD INCOMPATIBLE:** Remove obsolete descriptor schema: ``transmart``.
- **BACKWARD INCOMPATIBLE:** Remove tools which are not used by any process:
  ``clustering_leaf_ordering.py``, ``go_genesets.py``, ``VCF_ad_extract.py``,
  ``volcanoplot.py``, ``xgff.py``, ``xgtf2gff.py``.
- **BACKWARD INCOMPATIBLE:** Management command for inserting features and
  mappings requires PostgreSQL version 9.5 or newer
- Update the meta data like name, description, category, etc. of most of the
  processes
- Speed-up management command for inserting mappings
- Change location of cufflinks to Google Drive for resolwebio/rnaseq Docker
  build
- Calculate alignment statistics for the uploaded alignment (.bam) file in the
  ``upload-bam``, ``upload-bam-indexed`` and ``upload-bam-secondary`` processes.
- Annotation (GTF/GFF3) file input is now optional for the creation of the
  STAR genome index files. Annotation file can be used at the alignment stage
  to supplement the genome indices with the set of known features.
- Trigger process warning instead of process error in the case when
  ``bamtobigwig.sh`` scripts detects an empty .bam file.
- Set the default reads length filtering parameter to 30 bp in the
  ``rna-seq-bbduk-star-featurecounts`` and ``kapa-rna-seq-bbduk-star-featurecounts``
  experiment descriptor schema. Expand the kit selection choice options in the
  latter descriptor schema.

Added
-----
- Add ``MultiQC (1.6.0)`` and ``Seqtk (1.2-r94)`` to the
  ``resolwebio/utils:1.5.0`` Docker image
- Add ``sample2`` descriptor schema which is the successor of the original
  ``sample`` and ``reads`` descriptor schemas
- Add bedToBigBed and Tabix to resolwebio/rnaseq:3.7.0 docker image
- Add ``HS Panel`` choice option to the ``amplicon-master-file`` descriptor
  schema
- Add MultiQC process
- Add process for the Seqtk tool ``sample`` sub-command. This process allows
  sub-sampling of ``.fastq`` files using either a fixed number of reads or the
  ratio of the input file.
- Add MultiQC analysis step to the ``workflow-bbduk-star-featurecounts-single``
  and ``workflow-bbduk-star-featurecounts-single`` processes.
- Add ``workflow-bbduk-star-featurecounts-qc-single`` and
  ``workflow-bbduk-star-featurecounts-qc-paired`` processes which support
  MultiQC analysis, input reads down-sampling (using Seqtk) and rRNA
  sequence detection using STAR aligner.
- Add to ``resolwebio/chipseq`` Docker image:

  - ``bedtools (2.25.0-1)``
  - ``gawk (1:4.1.3+dfsg-0.1)``
  - ``picard-tools (1.113-2)``
  - ``run_spp.R (1.2) (as spp)``
  - ``SPP (1.14)``
- Add ``regtools-junctions-annotate`` process that annotates novel splice
  junctions.
- Add ``background`` relation type to fixtures

Fixed
-----
- Track ``source`` information in the ``upload-fasta-nucl`` process.
- When STAR aligner produces an empty alignment file, re-sort the alignment
  file to allow successful indexing of the output ``.bam`` file.
- Create a symbolic link to the alignment file in the ``feature_counts`` process,
  so that relative path is used in the quantification results. This prevent the
  FeatureCounts output to be listed as a separate sample in the MultiQC reports.
- Fix handling of expression objects in ``archive-samples`` process


===================
12.0.0 - 2018-08-13
===================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Require Resolwe 11.x
- **BACKWARD INCOMPATIBLE:** Use read count instead of sampling rate
  in strandedness detection
- **BACKWARD INCOMPATIBLE:** Remove ``genome`` input from ``rose2``
  process and automate its selection
- **BACKWARD INCOMPATIBLE:** Refactor ``cutadapt-paired`` process
- **BACKWARD INCOMPATIBLE:** Improve leaf ordering performance in gene and
  sample hierarchical clustering. We now use exact leaf ordering which has
  been recently added to ``scipy`` instead of an approximate in-house
  solution based on nearest neighbor algorithm. Add informative warning
  and error messages to simplify troubleshooting with degenerate datasets.
- Remove ``igvtools`` from ``resolwebio/utils`` Docker image
- Improve helper text and labels in processes used for sequencing data upload
- Allow using custom adapter sequences in the
  ``workflow-bbduk-star-featurecounts-single`` and
  ``workflow-bbduk-star-featurecounts-paired`` processes
- Change chromosome names from ENSEMBL / NCBI to UCSC (example: "1" to
  "chr1") in BigWig files. The purpose of this is to enable viewing BigWig
  files in UCSC genome browsers for files aligned with ENSEBML or NCBI genome.
  This change is done by adding script bigwig_chroms_to_ucsc.py to
  bamtobigwig.sh script.
- Reduce RAM requirement in SRA import processes

Added
-----
- Add two-pass mode to ``alignment-star`` process
- Add ``regtools (0.5.0)`` to ``resolwebio/rnaseq`` Docker image
- Add KAPA experiment descriptor schema
- Add ``resdk`` Python 3 package to ``resolwebio/utils`` Docker image
- Add to ``cutadapt-single`` process an option to discard reads having more
  'N' bases than specified.
- Add workflows for single-end ``workflow-cutadapt-star-featurecounts-single``
  and paired-end reads ``workflow-cutadapt-star-featurecounts-paired``.
  Both workflows consist of preprocessing with Cutadapt, alignment
  with STAR two pass mode and quantification with featureCounts.
- Add descriptor schema ``rna-seq-cutadapt-star-featurecounts``

Fixed
-----
- **BACKWARD INCOMPATIBLE:** Fix the ``stitch`` parameter handling in
  ``rose2``
- fix ``upload-gtf`` to create JBrowse track only if GTF file is ok
- Pin ``sra-toolkit`` version to 2.9.0 in ``resolwebio/utils`` Docker image.
- Fix and improve ``rose2`` error messages
- Fail gracefully if bam file is empty when producing bigwig files
- Fail gracefully if there are no matches when mapping chromosome names


===================
11.0.0 - 2018-07-17
===================

Changed
-------
- **BACKWARD INCOMPATIBLE:** Remove management command module
- **BACKWARD INCOMPATIBLE:** Remove filtering of genes with low expression
  in PCA analysis
- **BACKWARD INCOMPATIBLE:** Remove obsolete RNA-seq DSS process
- Expand error messages in ``rose2`` process
- Check for errors during download of FASTQ files and use
  ``resolwebio/utils:1.3.0`` Docker image in import SRA process
- Increase Feature's full name's max length to 350 to support a long full
  name of "Complement C3 Complement C3 beta chain C3-beta-c Complement C3
  alpha chain C3a anaphylatoxin Acylation stimulating protein Complement C3b
  alpha' chain Complement C3c alpha' chain fragment 1 Complement C3dg
  fragment Complement C3g fragment Complement C3d fragment Complement C3f
  fragment Complement C3c alpha' chain fragment 2" in Ensembl

Added
-----
- Add `exp_set` and `exp_set_json` output fields to expression processes:

  - ``feature_counts``
  - ``htseq-count``
  - ``htseq-count-raw``
  - ``rsem``
  - ``upload-expression``
  - ``upload-expression-cuffnorm``
  - ``upload-expression-star``
- Add 'Masking BED file' input to ``rose2`` process which allows
  masking reagions from the analysis
- Add ``filtering.outFilterMismatchNoverReadLmax`` input to
  ``alignment-star`` process
- Add mappings from ENSEMBL or NCBI to UCSC chromosome names to
  ``resolwebio/rnaseq:3.5.0`` docker image

Fixed
-----
- Fix peaks BigBed output in ``macs14`` process
- Remove duplicated forward of ``alignIntronMax`` input field in
  BBDuk - STAR - featureCounts workflow
- Make ``cuffnorm`` process attach correct expression data objects to
  samples
- Fix ``upload-gtf`` in a way that GTF can be shown in JBrowse. Because
  JBrowse works only with GFF files, input GTF is converted to GFF from
  which JBrowse track is created.


===================
10.0.1 - 2018-07-06
===================

Fixed
-----
- Fix ``bamtobigwig.sh`` to timeout the ``bamCoverage`` calculation after
  defined time


===================
10.0.0 - 2018-06-19
===================

Added
-----
- Add to ``resolwebio/chipseq`` Docker image:

  - ``Bedops (v2.4.32)``
  - ``Tabix (v1.8)``
  - ``python3-pandas``
  - ``bedGraphToBigWig (kent-v365)``
  - ``bedToBigBed (kent-v365)``
- Add to ``resolwebio/rnaseq:3.2.0`` Docker image:

  - ``genometools (1.5.9)``
  - ``igvtools (v2.3.98)``
  - ``jbrowse (v1.12.0)``
  - ``Bowtie (v1.2.2)``
  - ``Bowtie2 (v2.3.4.1)``
  - ``BWA (0.7.17-r1188)``
  - ``TopHat (v2.1.1)``
  - ``Picard Tools (v2.18.5)``
  - ``bedGraphToBigWig (kent-v365)``
- Add Debian package ``file`` to ``resolwebio/rnaseq:3.3.0`` Docker image
- Support filtering by type on feature API endpoint
- Add BigWig output field to following processes:

  - ``alignment-bowtie``
  - ``alignment-bowtie2``
  - ``alignment-tophat2``
  - ``alignment-bwa-mem``
  - ``alignment-bwa-sw``
  - ``alignment-bwa-aln``
  - ``alignment-hisat2``
  - ``alignment-star``
- Add Jbrowse track output field to ``upload-genome`` processor.
- Use ``reslowebio/rnaseq`` Docker image and add Jbrowse track and IGV
  sorting and indexing to following processes:

  - ``upload-gff3``
  - ``upload-gtf``
  - ``gff-to-gtf``
- Add Tabix index for Jbrowse to ``upload-bed`` processor and use
  ``reslowebio/rnaseq`` Docker image
- Add BigWig, BigBed and JBrowse track outputs to ``macs14`` process
- Add Species and Build outputs to ``rose2`` process
- Add Species, Build, BigWig, BigBed and JBrowse track outputs to ``macs2``
  process
- Add ``scipy`` (v1.1.0) Python 3 package to ``resolwebio/utils`` Docker image

Changed
-------
- **BACKWARD INCOMPATIBLE:** Drop support for Python 3.4 and 3.5
- **BACKWARD INCOMPATIBLE:** Require Resolwe 10.x
- **BACKWARD INCOMPATIBLE:** Upgrade to Django Channels 2
- **BACKWARD INCOMPATIBLE:** Count fragments (or templates) instead of reads
  by default in ``featureCounts`` process and
  ``BBDuk - STAR - featureCounts`` pipeline. The change applies only to
  paired-end data.
- **BACKWARD INCOMPATIBLE:** Use ``resolwebio/rnaseq:3.2.0`` Docker image
  in the following processes that output reads:

  - ``upload-fastq-single``
  - ``upload-fastq-paired``
  - ``files-to-fastq-single``
  - ``files-to-fastq-paired``
  - ``reads-merge``
  - ``bbduk-single``
  - ``bbduk-paired``
  - ``cutadapt-single``
  - ``cutadapt-paired``
  - ``cutadapt-custom-single``
  - ``cutadapt-custom-paired``
  - ``trimmomatic-single``
  - ``trimmomatic-paired``.

  This change unifies the version of ``FastQC`` tool (0.11.7) used for
  quality control of reads in the aforementioned processes. The new Docker
  image comes with an updated version of Cutadapt (1.16) which affects
  the following processes:

  - ``cutadapt-single``
  - ``cutadapt-paired``
  - ``cutadapt-custom-single``
  - ``cutadapt-custom-paired``.

  The new Docker image includes also an updated version of Trimmomatic (0.36)
  used in the following processes:

  - ``upload-fastq-single``
  - ``upload-fastq-paired``
  - ``files-to-fastq-single``
  - ``files-to-fastq-paired``
  - ``trimmomatic-single``
  - ``trimmomatic-paired``.
- **BACKWARD INCOMPATIBLE:** Change Docker image in ``alignment-subread``
  from ``resolwebio/legacy:1.0.0`` with Subread (v1.5.1) to
  ``resolwebio/rnaseq:3.2.0`` with Subread (v1.6.0). ``--multiMapping`` option
  was added instead of ``--unique_reads``. By default aligner report uniquely
  mapped reads only.
- Update ``wigToBigWig`` to kent-v365 version  in ``resolwebio/chipseq``
  Docker image
- Change paths in HTML amplicon report template in ``resolwebio/dnaseq``
  Docker image
- Move assay type input in BBDuk - STAR - featureCounts pipeline descriptor
  schema to advanced options
- Use ``resolwebio/rnaseq:3.2.0`` Docker image with updated versions of tools
  instead of ``resolwebio/legacy:1.0.0`` Docker image in following processes:

  - ``alignment-bowtie`` with Bowtie (v1.2.2) instead of Bowtie (v1.1.2)
  - ``alignment-bowtie2`` with Bowtie2 (v2.3.4.1) instead of Bowtie2 (v2.2.6)
  - ``alignment-tophat2`` with TopHat (v2.1.1) instead of TopHat (v2.1.0)
  - ``alignment-bwa-mem``, ``alignment-bwa-sw` and ``alignment-bwa-aln``
    with BWA (v0.7.17-r1188) instead of BWA (v0.7.12-r1039)
  - ``alignment-hisat2`` with HISAT2 (v2.1.0) instead of HISAT2 (v2.0.3-beta)
  - ``upload-genome``
- Use ``resolwebio/base:ubuntu-18.04`` Docker image as a base image in
  ``resolwebio/utils`` Docker image
- Update Python 3 packages in ``resolwebio/utils`` Docker image:

  - ``numpy`` (v1.14.4)
  - ``pandas`` (v0.23.0)
- Replace ``bedgraphtobigwig`` with ``deepTools`` in ``resolwebio/rnaseq``
  Docker image, due to faster performance
- Use ``resolwebio/rnaseq:3.3.0`` Docker image in ``alignment-star-index``
  with STAR (v2.5.4b)

Fixed
-----
- Make management commands use a private random generator instance
- Fix output ``covplot_html`` of ``coveragebed`` process
- Fix process ``archive-samples`` and ``amplicon-archive-multi-report`` to
  correctly handle nested file paths
- Change ``rose2`` and ``chipseq-peakscore`` to work with ``.bed`` or
  ``.bed.gz`` input files
- Fix the ``expression-aggregator`` process so that it tracks the
  ``species`` of the input expression data
- Fix ``bamtobigwig.sh`` to use ``deepTools`` instead of ``bedtools`` with
  ``bedgraphToBigWig`` due to better time performance


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
