"""Base quality score recalibration."""

import os
import shutil

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


class BQSR(Process):
    """A two pass process of BaseRecalibrator and ApplyBQSR from GATK.

    See GATK website for more information on BaseRecalibrator.

    It is possible to modify read group using GATK's AddOrReplaceGroups through Replace read groups in BAM
    (``read_group``) input field.
    """

    slug = "bqsr"
    name = "BaseQualityScoreRecalibrator"
    process_type = "data:alignment:bam:bqsr:"
    version = "2.5.1"
    category = "GATK"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
    }
    data_name = "{{ bam|name|default('?') }}"

    class Input:
        """Input fields to perform Base quality score recalibration."""

        bam = DataField("alignment:bam", label="BAM file containing reads")
        reference = DataField("seq:nucleotide", label="Reference genome file")
        known_sites = ListField(
            DataField(
                data_type="variants:vcf",
                description="One or more databases of known polymorphic sites used to exclude regions around known "
                "polymorphisms from analysis.",
            ),
            label="List of known sites of variation",
        )
        intervals = DataField(
            data_type="bed",
            required=False,
            label="One or more genomic intervals over which to operate.",
            description="This field is optional, but it can speed up the process by restricting calculations to "
            "specific genome regions.",
        )
        read_group = StringField(
            label="Replace read groups in BAM",
            description="Replace read groups in a BAM file.This argument enables the user to replace all read groups "
            "in the INPUT file with a single new read group and assign all reads to this read group in "
            "the OUTPUT BAM file. Addition or replacement is performed using Picard's "
            "AddOrReplaceReadGroups tool. Input should take the form of -name=value delimited by a "
            '";", e.g. "-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1". See tool\'s '
            "documentation for more information on tag names. Note that PL, LB, PU and SM are require "
            "fields. See caveats of rewriting read groups in the documentation.",
            default="",
        )
        validation_stringency = StringField(
            label="Validation stringency",
            description="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
            "can improve performance when processing a BAM file in which variable-length data (read, "
            "qualities, tags) do not otherwise need to be decoded. Default is STRICT. This setting is "
            "used in BaseRecalibrator and ApplyBQSR processes.",
            choices=[
                ("STRICT", "STRICT"),
                ("LENIENT", "LENIENT"),
                ("SILENT", "SILENT"),
            ],
            default="STRICT",
        )

        class Advanced:
            """Advanced options."""

            use_original_qualities = BooleanField(
                label="Use the base quality scores from the OQ tag",
                description="This flag tells GATK to use the original base qualities "
                "(that were in the data before BQSR/recalibration) which are stored in the OQ tag, if they are "
                "present, rather than use the post-recalibration quality scores. If no OQ tag is present for a "
                "read, the standard qual score will be used.",
                default=False,
            )
            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to BaseQualityScoreRecalibrator."""

        bam = FileField(label="Base quality score recalibrated BAM file")
        bai = FileField(label="Index of base quality score recalibrated BAM file")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")
        recal_table = FileField(label="Recalibration tabled")

    def run(self, inputs, outputs):
        """Run the analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        # Prepare output file names.
        bam = os.path.basename(inputs.bam.output.bam.path)
        file_name = os.path.splitext(os.path.basename(inputs.bam.output.bam.path))[0]
        bam_rg = f"{file_name}_RG.bam"

        species = inputs.bam.output.species

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        # Parse read_group argument from a string, delimited by a ; and =
        # into a form that will be accepted by AddOrReplaceReadGroups tool.
        # E.g. '-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1' should become
        # ['-LB', 'DAB', '-PL', 'Illumina', '-PU', 'barcode', '-SM', 'sample1']
        # prepended by INPUT and OUTPUT.
        if inputs.read_group:
            arrg = [
                "--java-options",
                f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                "--INPUT",
                f"{inputs.bam.output.bam.path}",
                "--VALIDATION_STRINGENCY",
                f"{inputs.validation_stringency}",
                "--OUTPUT",
                f"{bam_rg}",
                "--TMP_DIR",
                TMPDIR,
            ]

            present_tags = []
            for x in inputs.read_group.split(";"):
                split_tag = x.split("=")
                arrg.extend(split_tag)
                present_tags.append(split_tag[0])

            # Make sure all arguments to read_group are valid.
            all_tags = {
                "-LB",
                "-PL",
                "-PU",
                "-SM",
                "-CN",
                "-DS",
                "-DT",
                "-FO",
                "-ID",
                "-KS",
                "-PG",
                "-PI",
                "-PM",
                "-SO",
            }
            present_tag_set = set(present_tags)
            check_all_tags = present_tag_set.issubset(all_tags)
            if not check_all_tags:
                self.error("One or more read_group argument(s) improperly formatted.")

            # Check that there are no double entries of arguments to read_group.
            if len(present_tag_set) != len(present_tags):
                self.error("You have duplicate tags in read_group argument.")

            # Check that all mandatory arguments to read_group are present.
            mandatory_tags = {"-LB", "-PL", "-PU", "-SM"}
            check_tags = mandatory_tags.issubset(present_tag_set)
            if not check_tags:
                self.error(
                    "Missing mandatory read_group argument(s) (-PL, -LB, -PU and -SM are mandatory)."
                )

            Cmd["gatk"]["AddOrReplaceReadGroups"](arrg)
        else:
            shutil.copy2(inputs.bam.output.bam.path, bam_rg)

        # Make sure the file is indexed.
        Cmd["samtools"]["index"](bam_rg)

        recal_table = f"{file_name}_recalibration.table"
        br_inputs = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--input",
            f"{bam_rg}",
            "--output",
            f"{recal_table}",
            "--reference",
            f"{inputs.reference.output.fasta.path}",
            "--read-validation-stringency",
            f"{inputs.validation_stringency}",
            "--tmp-dir",
            TMPDIR,
        ]
        if inputs.intervals:
            br_inputs.extend(["--intervals", f"{inputs.intervals.output.bed.path}"])

        if inputs.advanced.use_original_qualities:
            br_inputs.append("--use-original-qualities")

        # Add known sites to the input parameters of BaseRecalibrator.
        for site in inputs.known_sites:
            br_inputs.extend(["--known-sites", f"{site.output.vcf.path}"])

        # Prepare bqsr recalibration file.
        Cmd["gatk"]["BaseRecalibrator"](br_inputs)

        self.progress(0.5)

        # Apply base recalibration.
        ab_inputs = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--input",
            f"{bam_rg}",
            "--output",
            f"{bam}",
            "--reference",
            f"{inputs.reference.output.fasta.path}",
            "--bqsr-recal-file",
            f"{recal_table}",
            "--read-validation-stringency",
            f"{inputs.validation_stringency}",
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.advanced.use_original_qualities:
            ab_inputs.append("--use-original-qualities")

        Cmd["gatk"]["ApplyBQSR"](ab_inputs)

        stats = f"{bam}_stats.txt"
        (Cmd["samtools"]["flagstat"][f"{bam}"] > stats)()

        self.progress(0.9)

        outputs.bam = bam
        outputs.bai = file_name + ".bai"
        outputs.stats = stats

        outputs.species = species
        outputs.build = inputs.bam.output.build
        outputs.recal_table = recal_table
