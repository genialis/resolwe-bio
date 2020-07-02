"""QoRTs QC."""
import json

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)

STRAND_CODES = {
    "IU": "non_specific",
    "U": "non_specific",
    "ISF": "forward",
    "OSF": "forward",
    "SF": "forward",
    "ISR": "reverse",
    "OSR": "reverse",
    "SR": "reverse",
}


class QortsQC(Process):
    """QoRTs QC analysis."""

    slug = "qorts-qc"
    name = "QoRTs QC"
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/rnaseq:4.10.0",},},
        "resources": {"cores": 1, "memory": 32768,},
    }
    data_name = "QoRTs QC report ({{alignment|sample_name}})"
    version = "1.4.1"
    process_type = "data:qorts:qc"
    category = "Other"
    entity = {
        "type": "sample",
        "input": "alignment",
    }
    scheduling_class = SchedulingClass.BATCH
    description = "Quality of RNA-seq Tool-Set."

    class Input:
        """Input fields."""

        alignment = DataField("alignment:bam", label="Alignment")
        annotation = DataField("annotation:gtf", label="GTF annotation")

        class Options:
            """Options."""

            stranded = StringField(
                label="Assay type",
                default="non_specific",
                choices=[
                    ("non_specific", "Strand non-specific"),
                    ("forward", "Strand-specific forward"),
                    ("reverse", "Strand-specific reverse"),
                    ("auto", "Detect automatically"),
                ],
            )

            cdna_index = DataField(
                "index:salmon",
                label="cDNA index file",
                required=False,
                hidden="options.stranded != 'auto'",
            )

            n_reads = IntegerField(
                label="Number of reads in subsampled alignment file",
                default=5000000,
                hidden="options.stranded != 'auto'",
            )

            maxPhredScore = IntegerField(label="Max Phred Score", required=False,)

            adjustPhredScore = IntegerField(label="Adjust Phred Score", required=False,)

        options = GroupField(Options, label="Options")

    class Output:
        """Output fields."""

        plot = FileField(label="QC multiplot")
        summary = FileField(label="QC summary")
        qorts_data = FileField(label="QoRTs report data")

    def run(self, inputs, outputs):
        """Run the analysis."""
        lib_strand = ""
        if inputs.options.stranded == "auto":
            detect_strandedness_inputs = [
                inputs.alignment.bam.path,
                inputs.options.n_reads,
                inputs.options.cdna_index.index.path,
                self.requirements.resources.cores,
            ]
            Cmd["detect_strandedness.sh"](detect_strandedness_inputs)
            try:
                lib_strand = STRAND_CODES[
                    json.load(open("results/lib_format_counts.json")).get(
                        "expected_format", ""
                    )
                ]
            except KeyError:
                self.error(
                    "Library strandedness autodetection failed. Use manual selection options instead."
                )

        # Default and required arguments
        args = [
            "--skipFunctions",
            "writeDESeq,writeDEXSeq",
            "--randomSeed",
            42,
            "--generatePlots",
            inputs.alignment.bam.path,
            inputs.annotation.annot.path,
            "qorts_output",
        ]

        optional_args = []

        # Detect if aligned reads in BAM file are of single or paired-end type
        # The samtools view command counts the number of reads with the SAM flag "read paired (0x1)"
        if (
            Cmd["samtools"]("view", "-c", "-f", "1", inputs.alignment.bam.path).strip()
            == "0"
        ):
            optional_args.append("--singleEnded")

        if inputs.options.stranded == "forward" or lib_strand == "forward":
            optional_args.extend(["--stranded", "--stranded_fr_secondstrand"])
        elif inputs.options.stranded == "reverse" or lib_strand == "reverse":
            optional_args.append("--stranded")

        if inputs.options.maxPhredScore or inputs.options.maxPhredScore == 0:
            optional_args.extend(["--maxPhredScore", inputs.options.maxPhredScore])

        if inputs.options.adjustPhredScore or inputs.options.adjustPhredScore == 0:
            optional_args.extend(
                ["--adjustPhredScore", inputs.options.adjustPhredScore]
            )

        memory_limit = "-Xmx{}g".format(self.requirements.resources.memory // 1024)

        # join optional and required arguments
        Cmd["QoRTs"]([memory_limit, "QC"] + optional_args + args)

        # Compress QoRTs output folder
        Cmd["zip"](["-r", "qorts_report.zip", "qorts_output"])

        outputs.plot = "qorts_output/QC.multiPlot.pdf"
        outputs.summary = "qorts_output/QC.summary.txt"
        outputs.qorts_data = "qorts_report.zip"
