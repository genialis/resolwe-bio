"""RNA-seq Variant Calling Workflow."""

from resolwe.process import (
    BooleanField,
    Data,
    DataField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    StringField,
)
from resolwe.process.models import Process as BioProcess


class WorkflowRnaseqVariantCalling(Process):
    """Identify variants in RNA-seq data.

    This pipeline follows GATK best practices recommendantions for
    variant calling with RNA-seq data.

    The pipeline steps include read alignment (STAR), data cleanup
    (MarkDuplicates), splitting reads that contain Ns in their cigar string
    (SplitNCigarReads), base quality recalibration (BaseRecalibrator,
    ApplyBQSR), variant calling (HaplotypeCaller), variant filtering
    (VariantFiltration) and variant annotation (SnpEff). The last step of the
    pipeline is process Mutations table which prepares variants for ReSDK
    VariantTables.

    There is also possibility to run the pipeline directly from BAM file.
    In this case, it is recommended that you use two-pass mode in STAR
    alignment as well as turn the option '--outSAMunmapped Within' on.
    """

    slug = "workflow-rnaseq-variantcalling"
    name = "RNA-seq Variant Calling Workflow"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1",
            },
        },
    }
    data_name = "RNA-seq Variants ({{ reads|name|default('?') if reads else bam|name|default('?') }})"
    version = "2.6.0"
    process_type = "data:workflow:rnaseq:variants"
    category = "Pipeline"
    entity = {
        "type": "sample",
    }

    class Input:
        """Input fields."""

        bam = DataField(
            data_type="alignment:bam:star",
            label="Input BAM file",
            description="Input BAM file that was computed with STAR aligner. "
            "It is highly recommended that two-pass mode was used for the "
            "alignment as well as '--outSAMunmapped Within' option if you want "
            "to use BAM file as an input.",
            disabled="reads",
            required=False,
        )
        reads = DataField(
            data_type="reads:fastq",
            label="Input sample (FASTQ)",
            description="Input data in FASTQ format.",
            disabled="bam",
            required=False,
        )
        preprocessing = BooleanField(
            label="Perform reads processing with BBDuk",
            default=True,
            description="If your reads have not been processed, set this to True.",
            disabled="bam",
        )
        ref_seq = DataField(
            data_type="seq:nucleotide", label="Reference FASTA sequence"
        )
        genome = DataField(
            data_type="index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
            disabled="bam",
            required=False,
        )
        dbsnp = DataField(
            data_type="variants:vcf",
            label="dbSNP file",
            description="File with known variants.",
        )
        indels = ListField(
            DataField(data_type="variants:vcf"),
            label="Known INDEL sites",
            required=False,
        )
        intervals = DataField(
            data_type="bed",
            label="Intervals (from BED file)",
            description="Use this option to perform the analysis over only part "
            "of the genome.",
            required=False,
        )
        clinvar = DataField(
            data_type="variants:vcf",
            label="ClinVar VCF file",
            description="[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a "
            "freely available, public archive of human genetic variants and "
            "interpretations of their significance to disease.",
            required=False,
        )
        geneset = DataField(
            data_type="geneset",
            label="Gene set",
            description="Select a gene set with genes you are interested in. "
            "Only variants of genes in the selected gene set will be in the "
            "output.",
            required=False,
            disabled="mutations",
        )
        mutations = ListField(
            StringField(),
            label="Gene and its mutations",
            description="Insert the gene you are interested in, together "
            "with mutations. First enter the name of the gene and then "
            "the mutations. Seperate gene from mutations with ':' and mutations "
            "with ','. Example of an input: 'KRAS: Gly12, Gly61'. Press enter "
            "after each input (gene + mutations). NOTE: Field only accepts "
            "three character amino acid symbols. If you use this option, "
            "the selected geneset will not be used for Mutations table process.",
            required=False,
            disabled="geneset",
        )

        class Bbduk:
            """Preprocessing options (BBDuk)."""

            adapters = ListField(
                DataField(data_type="seq:nucleotide"),
                label="Adapters",
                required=False,
                description="Provide a list of sequencing adapters files (.fasta) "
                "to be removed by BBDuk.",
            )
            custom_adapter_sequences = ListField(
                StringField(),
                label="Custom adapter sequences",
                required=False,
                default=[],
                description="Custom adapter sequences can be specified by inputting them "
                "one by one and pressing Enter after each sequence.",
            )
            kmer_length = IntegerField(
                label="K-mer length [k=]",
                default=23,
                description="Kmer length used for finding contaminants. "
                "Contaminants shorter than kmer length will not be found. "
                "Kmer length must be at least 1.",
            )
            min_k = IntegerField(
                label="Minimum k-mer length at right end of reads used for trimming [mink=]",
                default=11,
                disabled="bbduk.adapters.length === 0 && bbduk.custom_adapter_sequences.length === 0",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for k-mers [hammingdistance=]",
                default=1,
                description="Hamming distance i.e. the number of mismatches allowed in the kmer.",
            )
            maxns = IntegerField(
                label="Max Ns after trimming [maxns=]",
                default=-1,
                description="If non-negative, reads with more Ns than this (after trimming) will be discarded.",
            )
            trim_quality = IntegerField(
                label="Average quality below which to trim region [trimq=]",
                default=28,
                description="Phred algorithm is used, which is more accurate than naive trimming.",
            )
            min_length = IntegerField(
                label="Minimum read length [minlength=]",
                default=30,
                description="Reads shorter than minimum read length after trimming are discarded.",
            )
            quality_encoding_offset = StringField(
                label="Quality encoding offset",
                choices=[
                    ("33", "Sanger / Illumina 1.8+"),
                    ("64", "Illumina up to 1.3+, 1.5+"),
                    ("auto", "Auto"),
                ],
                default="auto",
                description="Quality encoding offset for input FASTQ files.",
            )
            ignore_bad_quality = BooleanField(
                label="Ignore bad quality",
                default=False,
                description="Don't crash if quality values appear to be incorrect.",
            )

        class Alignment:
            """Alignment options (STAR)."""

            two_pass_mode = BooleanField(
                label="Use two pass mode [--twopassMode]",
                default=True,
                description="Use two-pass maping instead of first-pass only. In two-pass mode we "
                "first perform first-pass mapping, extract junctions, insert them into genome "
                "index, and re-map all reads in the second mapping pass.",
            )
            out_unmapped = BooleanField(
                label="Output unmapped reads (SAM) [--outSAMunmapped Within]",
                default=True,
                description="Output of unmapped reads in the SAM format.",
            )
            align_end_alignment = StringField(
                label="Read ends alignment [--alignEndsType]",
                choices=[
                    ("Local", "Local"),
                    ("EndToEnd", "EndToEnd"),
                    ("Extend5pOfRead1", "Extend5pOfRead1"),
                    ("Extend5pOfReads12", "Extend5pOfReads12"),
                ],
                description="Type of read ends alignment (default: Local). Local: standard local "
                "alignment with soft-clipping allowed. EndToEnd: force end-to-end read alignment, "
                "do not soft-clip. Extend5pOfRead1: fully extend only the 5p of the read1, all "
                "other ends: local alignment. Extend5pOfReads12: fully extend only the 5' of the "
                "both read1 and read2, all other ends use local alignment.",
                default="Local",
            )

        class BAMProcessing:
            """BAM file processing options."""

            read_group = StringField(
                label="Replace read groups in BAM",
                description="Replace read groups in a BAM file. This argument enables the "
                "user to replace all read groups in the INPUT file with a single new read "
                "group and assign all reads to this read group in the OUTPUT BAM file. "
                "Addition or replacement is performed using Picard's AddOrReplaceReadGroups "
                "tool. Input should take the form of -name=value delimited by a "
                '";", e.g. "-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1". '
                "See tool's documentation for more information on tag names. Note that "
                "PL, LB, PU and SM are required fields. See caveats of rewriting read groups "
                "in the documentation.",
                default="-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1",
            )

        class HaplotypeCaller:
            """GATK HaplotypeCaller options."""

            stand_call_conf = IntegerField(
                label="Min call confidence threshold",
                default=20,
                description="The minimum phred-scaled confidence threshold at which "
                "variants should be called.",
            )
            soft_clipped = BooleanField(
                label="Do not analyze soft clipped bases in the reads",
                default=True,
                description="Suitable option for RNA-seq variant calling.",
            )
            interval_padding = IntegerField(
                label="Interval padding",
                default=100,
                description="Amount of padding (in bp) to add to each interval "
                "you are including. The recommended value is 100. Set to 0 if you "
                "want to turn it off.",
                hidden="!intervals",
            )

        class VariantFiltration:
            """GATK VariantFiltration options."""

            filter_expressions = ListField(
                StringField(),
                label="Expressions used with INFO fields to filter",
                description="VariantFiltration accepts any number of JEXL expressions "
                "(so you can have two named filters by using --filter-name One "
                "--filter-expression 'X < 1' --filter-name Two --filter-expression 'X > 2'). "
                "It is preferable to use multiple expressions, each specifying an individual "
                "filter criteria, to a single compound expression that specifies multiple "
                "filter criteria. Input expressions one by one and press ENTER after each "
                "expression. Examples of filter expression: 'FS > 30', 'DP > 10'.",
                default=["FS > 30.0", "QD < 2.0"],
            )
            filter_name = ListField(
                StringField(),
                label="Names to use for the list of filters",
                description="This name is put in the FILTER field for variants that get "
                "filtered. Note that there must be a 1-to-1 mapping between filter expressions "
                "and filter names. Input expressions one by one and press ENTER after each name. "
                "Warning: filter names should be in the same order as filter expressions. "
                "Example: you specified filter expressions 'FS > 30' and 'DP > 10', now "
                "specify filter names 'FS' and 'DP'.",
                default=["FS", "QD"],
            )
            genotype_filter_expressions = ListField(
                StringField(),
                label="Expressions used with FORMAT field to filter",
                description="Similar to the INFO field based expressions, but used on the FORMAT "
                "(genotype) fields instead. VariantFiltration will add the sample-level FT tag to "
                "the FORMAT field of filtered samples (this does not affect the record's FILTER tag). "
                "One can filter normally based on most fields (e.g. 'GQ < 5.0'), but the GT "
                "(genotype) field is an exception. We have put in convenience methods so that "
                "one can now filter out hets ('isHet == 1'), refs ('isHomRef == 1'), or homs "
                "('isHomVar == 1'). Also available are expressions isCalled, isNoCall, isMixed, "
                "and isAvailable, in accordance with the methods of the Genotype object. "
                "To filter by alternative allele depth, use the expression: 'AD.1 < 5'. This "
                "filter expression will filter all the samples in the multi-sample VCF file.",
                default=["AD.1 < 5.0"],
            )
            genotype_filter_name = ListField(
                StringField(),
                label="Names to use for the list of genotype filters",
                description="Similar to the INFO field based expressions, but used on the FORMAT "
                "(genotype) fields instead. Warning: filter names should be in the same order as "
                "filter expressions.",
                default=["AD"],
            )
            mask = DataField(
                data_type="variants:vcf",
                label="Input mask",
                description="Any variant which overlaps entries from the provided "
                "mask file will be filtered.",
                required=False,
            )
            mask_name = StringField(
                label="The text to put in the FILTER field if a 'mask' is provided",
                description="When using the mask file, the mask name will be annotated in "
                "the variant record.",
                required=False,
                disabled="!variant_filtration.mask",
            )

        class SnpEff:
            """SnpEff options."""

            filtering_options = StringField(
                label="SnpEff filtering expressions",
                description="Filter annotated VCF file using arbitraty expressions."
                "Examples of filtering expressions: '(ANN[*].GENE = 'PSD3')' "
                "or '( REF = 'A' )' or "
                "'(countHom() > 3) | (( exists INDEL ) & (QUAL >= 20)) | (QUAL >= 30 )'."
                "For more information checkout the official documentation of [SnpSift]"
                "(https://pcingola.github.io/SnpEff/ss_filter/)",
                required=False,
            )

        class MutationsTable:
            """Mutations table options."""

            vcf_fields = ListField(
                StringField(),
                label="Select VCF fields",
                description="The name of a standard VCF field or an "
                "INFO field to include in the output table. "
                "The field can be any standard VCF column (e.g. CHROM, ID, QUAL) "
                "or any annotation name in the INFO field (e.g. AC, AF). "
                "Required fields are CHROM, POS, ID, REF and ANN. If your variants "
                "file was annotated with clinvar information then fields CLNDN, "
                "CLNSIG and CLNSIGCONF might be of your interest.",
                default=[
                    "CHROM",
                    "POS",
                    "ID",
                    "QUAL",
                    "REF",
                    "ALT",
                    "FILTER",
                    "ANN",
                    "CLNDN",
                    "CLNSIG",
                    "QD",
                ],
            )
            ann_fields = ListField(
                StringField(),
                label="ANN fields to use",
                description="Only use specific fields from the SnpEff ANN "
                "field. All available fields: Allele | Annotation | Annotation_Impact "
                "| Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType "
                "| Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length "
                "| AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ."
                "Fields are seperated by '|'. For more information, follow this [link]"
                "(https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files).",
                default=[
                    "Allele",
                    "Annotation",
                    "Annotation_Impact",
                    "Gene_Name",
                    "Feature_ID",
                    "HGVS.p",
                ],
            )
            split_alleles = BooleanField(
                label="Split multi-allelic records into multiple lines",
                description="By default, a variant record with multiple "
                "ALT alleles will be summarized in one line, with per "
                "alt-allele fields (e.g. allele depth) separated by commas."
                "This may cause difficulty when the table is loaded by "
                "an R script, for example. Use this flag to write multi-allelic "
                "records on separate lines of output.",
                default=True,
            )
            show_filtered = BooleanField(
                label="Include filtered records in the output",
                default=True,
                description="Include filtered records in the output of the GATK "
                "VariantsToTable.",
            )
            gf_fields = ListField(
                StringField(),
                label="Include FORMAT/sample-level fields. Note: If you specify DP "
                "from genotype field, it will overwrite the original DP field. "
                "By default fields GT (genotype), GQ (genotype quality) AD (allele depth), "
                "DP (depth at the sample level), FT (sample-level filter) "
                "are included in the analysis.",
                default=["GT", "GQ", "AD", "DP", "FT"],
            )

        class Advanced:
            """Advanced options."""

            multiqc = BooleanField(
                label="Trigger MultiQC",
                default=False,
                hidden="!bam",
                description="If the input for the pipeline is BAM file "
                "that has been computed by the RNA-seq gene expression pipeline, than "
                "MultiQC object already exists for this sample, so there is no need "
                "for an additional MultiQC process. If the input for this pipeline is "
                "FASTQ, than MultiQC cannot be disabled.",
            )
            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of "
                "the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        bbduk = GroupField(Bbduk, label="Preprocessing with BBDuk", hidden="bam")
        alignment = GroupField(Alignment, label="Alignment with STAR", hidden="bam")
        bam_processing = GroupField(BAMProcessing, label="Processing of BAM file")
        haplotype_caller = GroupField(
            HaplotypeCaller, label="Options for HaplotypeCaller"
        )
        variant_filtration = GroupField(
            VariantFiltration, label="Options for GATK VariantFiltration"
        )
        snpeff = GroupField(SnpEff, label="SnpEff options")
        mutations_table = GroupField(
            MutationsTable, label="Options for Mutations table"
        )
        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields."""

        # Workflow does not have outputs.

    def run(self, inputs, outputs):
        """Run the workflow."""

        if not inputs.bam and not inputs.reads:
            self.error("Please select input BAM or FASTQ file.")

        if inputs.reads and not inputs.genome:
            self.error("Please select STAR index for the alignment.")

        if inputs.genome:
            if (
                inputs.ref_seq.output.species != inputs.genome.output.species
                or inputs.dbsnp.output.species != inputs.genome.output.species
            ):
                self.error(
                    "All input files must be from the same species. "
                    f"FASTA reference sequence is from {inputs.ref_seq.output.species}, "
                    f"STAR index is from {inputs.genome.output.species} and "
                    f"DBSNP file is from {inputs.dbsnp.output.species}"
                )

            if (
                inputs.ref_seq.output.build != inputs.genome.output.build
                or inputs.dbsnp.output.build != inputs.genome.output.build
            ):
                self.error(
                    "All input files must have the same build. "
                    f"FASTA reference sequence is based on {inputs.ref_seq.output.build}, "
                    f"STAR index is based on {inputs.genome.output.build} and "
                    f"DBSNP file is based on {inputs.dbsnp.output.build}. "
                )

            if inputs.intervals:
                if inputs.intervals.output.species != inputs.genome.output.species:
                    self.error(
                        "All input files must be from the same species. "
                        f"Intervals BED file is from {inputs.intervals.output.species}, "
                        f"while STAR index is from {inputs.genome.output.species}."
                    )
                if inputs.intervals.output.build != inputs.genome.output.build:
                    self.error(
                        "All input files must have the same build. "
                        f"Intervals BED file is based on {inputs.intervals.output.build}, "
                        f"while STAR index is based on {inputs.genome.output.build}."
                    )

            if inputs.indels:
                for indel in inputs.indels:
                    if indel.output.species != inputs.genome.output.species:
                        self.error(
                            "All input files must be from the same species. "
                            f"File with INDELs is from {indel.output.species}, "
                            f"while STAR index is from {inputs.genome.output.species}."
                        )
                    if indel.output.build != inputs.genome.output.build:
                        self.error(
                            "All input file must have the same build. "
                            f"File with INDELs is based on {indel.output.build}, "
                            f"while STAR index is based on {inputs.genome.output.build}."
                        )

        if not inputs.mutations and not inputs.geneset:
            self.error(
                "Mutations or geneset were not specified. You must either enter desired "
                "mutations or select your geneset of interest."
            )

        if not all(
            field in inputs.mutations_table.vcf_fields
            for field in ["CHROM", "POS", "ID", "REF", "ANN"]
        ):
            self.error(
                "Input VCF fields do not contain all required values. "
                "Required fields are CHROM, POS, ID, REF and ANN."
            )

        if inputs.bam:
            if inputs.bam.input.two_pass_mapping.two_pass_mode != True:
                self.warning(
                    "Two-pass mode was not used in alignment with STAR. It is "
                    "highly recommended that you use two-pass mode for RNA-seq "
                    "variant calling."
                )
            if inputs.bam.input.output_options.out_unmapped != True:
                self.warning(
                    "It is recommended that you use parameter '--outSAMunmapped Within' "
                    "in STAR alignment."
                )

        if inputs.reads:
            name = inputs.reads.name
        elif inputs.bam:
            name = inputs.bam.name

        if inputs.reads and inputs.preprocessing:
            # BBDuk options used are the same as the default options
            # in workflow-bbduk-star-featurecounts-qc
            input_preprocessing = {
                "reads": inputs.reads,
                "min_length": inputs.bbduk.min_length,
                "reference": {
                    "sequences": inputs.bbduk.adapters or [],
                    "literal_sequences": inputs.bbduk.custom_adapter_sequences,
                },
                "processing": {
                    "kmer_length": inputs.bbduk.kmer_length,
                    "hamming_distance": inputs.bbduk.hamming_distance,
                },
                "operations": {
                    "quality_trim": "r",
                    "trim_quality": inputs.bbduk.trim_quality,
                    "quality_encoding_offset": inputs.bbduk.quality_encoding_offset,
                    "ignore_bad_quality": inputs.bbduk.ignore_bad_quality,
                    "maxns": inputs.bbduk.maxns,
                },
            }

            if inputs.bbduk.adapters or inputs.bbduk.custom_adapter_sequences:
                input_preprocessing["operations"]["k_trim"] = "r"
            else:
                input_preprocessing["operations"]["k_trim"] = "f"

            if inputs.bbduk.adapters or inputs.bbduk.custom_adapter_sequences:
                input_preprocessing["operations"]["min_k"] = inputs.bbduk.min_k
            else:
                input_preprocessing["operations"]["min_k"] = -1

            if inputs.reads.type.startswith("data:reads:fastq:single:"):
                slug_bbduk = "bbduk-single"
            elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
                input_preprocessing["operations"]["trim_pairs_evenly"] = True
                input_preprocessing["operations"]["trim_by_overlap"] = True
                slug_bbduk = "bbduk-paired"
            else:
                self.error("Wrong reads input type.")

            preprocessing = Data.create(
                process=BioProcess.get_latest(slug=slug_bbduk),
                input=input_preprocessing,
                name=f"Trimmed ({inputs.reads.name})",
            )

        if inputs.reads:
            input_star = {
                "reads": preprocessing if inputs.preprocessing else inputs.reads,
                "genome": inputs.genome,
                "two_pass_mapping": {
                    "two_pass_mode": inputs.alignment.two_pass_mode,
                },
                "output_options": {
                    "out_unmapped": inputs.alignment.out_unmapped,
                },
                "alignment": {
                    "align_end_alignment": inputs.alignment.align_end_alignment
                },
            }

            alignment = Data.create(
                process=BioProcess.get_latest(slug="alignment-star"),
                input=input_star,
                name=f"Aligned ({name})",
            )

        elif inputs.bam:
            alignment = inputs.bam

        input_preprocess = {
            "bam": alignment,
            "ref_seq": inputs.ref_seq,
            "known_sites": [inputs.dbsnp],
            "read_group": inputs.bam_processing.read_group,
        }
        if inputs.indels:
            for indel in inputs.indels:
                input_preprocess["known_sites"].append(indel)

        preprocess = Data.create(
            process=BioProcess.get_latest(slug="rnaseq-vc-preprocess"),
            input=input_preprocess,
            name=f"Preprocessed ({name})",
        )

        input_hc = {
            "alignment": preprocess,
            "genome": inputs.ref_seq,
            "dbsnp": inputs.dbsnp,
            "stand_call_conf": inputs.haplotype_caller.stand_call_conf,
            "mbq": 10,
            "advanced": {
                "soft_clipped": inputs.haplotype_caller.soft_clipped,
                "java_gc_threads": inputs.advanced.java_gc_threads,
                "max_heap_size": inputs.advanced.max_heap_size,
            },
        }

        if inputs.intervals:
            input_hc["intervals_bed"] = inputs.intervals
            input_hc["advanced"][
                "interval_padding"
            ] = inputs.haplotype_caller.interval_padding

        hc = Data.create(
            process=BioProcess.get_latest(slug="vc-gatk4-hc"),
            input=input_hc,
            name=f"Genotypes ({name})",
        )

        input_filtration = {
            "vcf": hc,
            "ref_seq": inputs.ref_seq,
            "filter_expressions": inputs.variant_filtration.filter_expressions,
            "filter_name": inputs.variant_filtration.filter_name,
            "genotype_filter_expressions": inputs.variant_filtration.genotype_filter_expressions,
            "genotype_filter_name": inputs.variant_filtration.genotype_filter_name,
            "advanced": {
                "cluster": 3,
                "window": 35,
                "java_gc_threads": inputs.advanced.java_gc_threads,
                "max_heap_size": inputs.advanced.max_heap_size,
            },
        }
        if inputs.variant_filtration.mask:
            if not inputs.variant_filtration.mask_name:
                self.error(
                    "If you specify a mask file, please specify 'mask name' - the text to "
                    "put in the FILTER field"
                )

            input_filtration["mask"] = inputs.variant_filtration.mask
            input_filtration["mask_name"] = inputs.variant_filtration.mask_name

        filtration = Data.create(
            process=BioProcess.get_latest(slug="gatk-variant-filtration-single"),
            input=input_filtration,
            name=f"Filtered genotypes ({name})",
        )

        snpeff_inputs = {
            "variants": filtration,
        }

        if "GRCh37" in inputs.ref_seq.output.build:
            snpeff_inputs["database"] = "GRCh37.75"
        elif "GRCh38" in inputs.ref_seq.output.build:
            snpeff_inputs["database"] = "GRCh38.109"

        if inputs.snpeff.filtering_options:
            snpeff_inputs["filtering_options"] = inputs.snpeff.filtering_options

        if inputs.clinvar:
            snpeff_inputs["dbsnp"] = inputs.clinvar

        snpeff = Data.create(
            process=BioProcess.get_latest(slug="snpeff-single"),
            input=snpeff_inputs,
            name=f"Annotated genotypes ({name})",
        )

        mutations_inputs = {
            "variants": snpeff,
            "vcf_fields": inputs.mutations_table.vcf_fields,
            "ann_fields": inputs.mutations_table.ann_fields,
            "bam": preprocess,
            "advanced": {
                "split_alleles": inputs.mutations_table.split_alleles,
                "show_filtered": inputs.mutations_table.show_filtered,
                "gf_fields": inputs.mutations_table.gf_fields,
            },
        }

        if inputs.mutations:
            mutations_inputs["mutations"] = inputs.mutations
        elif inputs.geneset:
            mutations_inputs["geneset"] = inputs.geneset

        Data.create(
            process=BioProcess.get_latest(slug="mutations-table"),
            input=mutations_inputs,
            name=f"Mutations table ({name})",
        )

        multiqc_tools = [alignment]

        if inputs.reads:
            multiqc_tools.append(inputs.reads)
            if inputs.preprocessing:
                multiqc_tools.append(preprocessing)

            input_multiqc = {"data": multiqc_tools}

            Data.create(
                process=BioProcess.get_latest(slug="multiqc"),
                input=input_multiqc,
            )

        elif inputs.bam and inputs.advanced.multiqc:
            input_multiqc = {"data": multiqc_tools}

            Data.create(
                process=BioProcess.get_latest(slug="multiqc"),
                input=input_multiqc,
            )
