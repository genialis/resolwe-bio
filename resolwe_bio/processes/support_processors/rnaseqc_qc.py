"""RNA-SeQC QC."""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Persistence,
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


def detect_strandedness(bam_path, n_reads, cdna_path, resources_cores):
    """Detect strandedness using SALMON tool."""
    detect_strandedness_inputs = [bam_path, n_reads, cdna_path, resources_cores]
    Cmd["detect_strandedness.sh"](detect_strandedness_inputs)
    try:
        lib_strand = STRAND_CODES[
            json.load(open("results/lib_format_counts.json")).get("expected_format", "")
        ]
        return lib_strand

    except KeyError:
        return None


def format_ucsc(annotation_path):
    """Convert UCSC annotation to a format suitable for parsing with collapse_annotation.py script."""

    def add_attributes(x):
        new_attrs = ""
        new_attrs += f'gene_name "{x.gene_id}"; '
        new_attrs += f'gene_id "{x.gene_id}"; '
        new_attrs += f'transcript_name "{x.gene_id}"; '
        new_attrs += f'transcript_id "{x.gene_id}"; '
        new_attrs += f'gene_type "{x.gene_id}"; '
        new_attrs += f'transcript_type "{x.gene_id}"; '

        new_attrs = new_attrs.strip()

        return new_attrs

    df = pd.read_csv(
        annotation_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        low_memory=False,
    )

    df = pd.DataFrame(
        {
            "chromosome": df[0],
            "feature_type": df[2],
            "start": df[3],
            "end": df[4],
            "strand": df[6],
            "gene_id": df[8].str.extract(r'gene_id "(.*?)"', expand=False).values,
        },
    )

    df.replace("", np.nan, inplace=True)
    df = df[df["feature_type"] == "exon"]

    columns = {
        "chromosome": "first",
        "feature_type": "first",
        "start": "min",
        "end": "max",
        "strand": "first",
        "gene_id": "first",
    }

    df["gene_id"] = (
        "chr"
        + df["chromosome"].astype(str)
        + "_"
        + df["strand"].astype(str)
        + "_"
        + df["gene_id"].astype(str)
    )

    chromosome_order = df["chromosome"].unique().tolist()
    hierarchy_order = ["gene", "transcript"]

    genes_df = df.groupby(["gene_id"]).agg(columns).assign(feature_type="gene")
    transcripts_df = (
        df.groupby(["gene_id"]).agg(columns).assign(feature_type="transcript")
    )
    df = pd.concat([genes_df, transcripts_df, df], ignore_index=True)
    df["attributes"] = df.apply(lambda x: add_attributes(x), axis=1)
    df["chromosome_order"] = df["chromosome"].map(lambda x: chromosome_order.index(x))
    df["hierarchy_order"] = df["feature_type"].map(
        lambda x: hierarchy_order.index(x) if x in hierarchy_order else float("inf")
    )

    df = df.sort_values(
        by=["chromosome_order", "gene_id", "hierarchy_order"], ascending=True
    )

    out_df = pd.DataFrame(
        {
            "chromosome": df["chromosome"],
            "annotation_source": "UCSC",
            "feature_type": df["feature_type"],
            "start": df["start"],
            "end": df["end"],
            "score": ".",
            "strand": df["strand"],
            "genomic_phase": ".",
            "attributes": df["attributes"],
        },
    )

    out_filename = "ucsc_formatted_annotation.gtf"

    out_df.to_csv(
        out_filename,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
    )

    return out_filename


class QcRnaseqc(Process):
    """RNA-SeQC QC analysis.

    An efficient new version of RNA-SeQC that computes a comprehensive set of metrics
    for characterizing samples processed by a wide range of protocols.
    It also quantifies gene- and exon-level expression,
    enabling effective quality control of large-scale RNA-seq datasets.

    More information can be found in the
    [GitHub repository](https://github.com/getzlab/rnaseqc)
    and in the [original paper](https://academic.oup.com/bioinformatics/article/37/18/3048/6156810?login=false).
    """

    slug = "rnaseqc-qc"
    name = "RNA-SeQC"
    process_type = "data:rnaseqc:qc"
    version = "2.0.0"
    category = "QC"
    data_name = "{{ alignment|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/qc:1.1.0"}
        },
        "resources": {"cores": 1, "memory": 32768, "storage": 100},
    }

    entity = {
        "type": "sample",
        "input": "alignment",
    }

    class Input:
        """Input fields."""

        alignment = DataField("alignment:bam", label="Input aligned reads (BAM file)")

        annotation = DataField(
            "annotation:gtf",
            label="Annotation file (GTF)",
            description="The input GTF file containing features to check the bam against. "
            "The file should include gene_id in the attributes column for all entries. "
            "During the process the file is formatted so the transcript_id matches the gene_id. "
            "Exons are merged to remove overlaps and exon_id field is then "
            "matched with gene_id including the consecutive exon number.",
        )

        class RnaseqcOptions:
            """RNA-SeQC options."""

            mapping_quality = IntegerField(
                label="Mapping quality [--mapping-quality]",
                default=255,
                description="Set the lower bound on read quality for exon coverage counting. "
                "Reads below this number are excluded from coverage metrics.",
            )

            base_mismatch = IntegerField(
                label="Base mismatch [--base-mismatch]",
                default=6,
                description="Set the maximum number of allowed mismatches "
                "between a read and the reference sequence. "
                "Reads with more than this number of mismatches "
                "are excluded from coverage metrics.",
            )

            offset = IntegerField(
                label="Offset [--offset]",
                default=150,
                description="Set the offset into the gene for the 3' and 5' windows in bias calculation. "
                "A positive value shifts the 3' and 5' windows towards each other, "
                "while a negative value shifts them apart.",
            )

            window_size = IntegerField(
                label="Window size [--window-size]",
                default=100,
                description="Set the offset into the gene for the 3' and 5' windows in bias calculation.",
            )

            gene_length = IntegerField(
                label="Window size [--gene-length]",
                default=600,
                description="Set the minimum size of a gene for bias calculation. "
                "Genes below this size are ignored in the calculation.",
            )

            detection_threshold = IntegerField(
                label="Detection threshold [--detection-threshold]",
                default=5,
                description="Number of counts on a gene to consider the gene 'detected'. "
                "Additionally, genes below this limit are excluded from 3' bias computation.",
            )

            exclude_chimeric = BooleanField(
                label="Exclude chimeric reads [--exclude-chimeric]",
                default=False,
                description="Exclude chimeric reads from the read counts.",
            )

        class StrandDetectionOptions:
            """Strand detection options."""

            stranded = StringField(
                label="Assay type [--stranded]",
                default="non_specific",
                choices=[
                    ("non_specific", "Strand non-specific"),
                    ("reverse", "Strand-specific reverse then forward"),
                    ("forward", "Strand-specific forward then reverse"),
                    ("auto", "Detect automatically"),
                ],
            )

            cdna_index = DataField(
                "index:salmon",
                label="cDNA index file",
                required=False,
                hidden="strand_detection_options.stranded != 'auto'",
            )

            n_reads = IntegerField(
                label="Number of reads in subsampled alignment file. "
                "Subsampled reads will be used in strandedness detection",
                default=5000000,
                hidden="strand_detection_options.stranded != 'auto'",
            )

        rnaseqc_options = GroupField(RnaseqcOptions, label="RNA-SeQC options")
        strand_detection_options = GroupField(
            StrandDetectionOptions, label="Strand detection options"
        )

    class Output:
        """Output fields."""

        metrics = FileField(label="metrics")

    def run(self, inputs, outputs):
        """Run the analysis."""

        bam_filename = Path(inputs.alignment.output.bam.path).name

        args = [
            "--mapping-quality",
            inputs.rnaseqc_options.mapping_quality,
            "--base-mismatch",
            inputs.rnaseqc_options.base_mismatch,
            "--offset",
            inputs.rnaseqc_options.offset,
            "--window-size",
            inputs.rnaseqc_options.window_size,
            "--gene-length",
            inputs.rnaseqc_options.gene_length,
            "--detection-threshold",
            inputs.rnaseqc_options.detection_threshold,
            "--verbose",
            "--verbose",
        ]

        if inputs.rnaseqc_options.exclude_chimeric:
            args.append("--exclude-chimeric")

        # Detect if aligned reads in BAM file are of single or paired-end type
        # The samtools view command counts the number of reads with the SAM flag "read paired (0x1)"
        if (
            Cmd["samtools"](
                "view", "-c", "-f", "1", inputs.alignment.output.bam.path
            ).strip()
            == "0"
        ):
            args.append("--unpaired")  # Required for single-end libraries
        else:
            lib_strand = ""
            if inputs.strand_detection_options.stranded == "auto":
                lib_strand = detect_strandedness(
                    bam_path=inputs.alignment.output.bam.path,
                    n_reads=inputs.strand_detection_options.n_reads,
                    cdna_path=inputs.strand_detection_options.cdna_index.output.index.path,
                    resources_cores=self.requirements.resources.cores,
                )
                if lib_strand is None:
                    self.error(
                        "Library strandedness autodetection failed. Use manual selection options instead."
                    )

            if (
                inputs.strand_detection_options.stranded == "forward"
                or lib_strand == "forward"
            ):
                args.extend(["--stranded", "FR"])
            elif (
                inputs.strand_detection_options.stranded == "reverse"
                or lib_strand == "reverse"
            ):
                args.extend(["--stranded", "RF"])

        self.progress(0.3)

        if inputs.annotation.output.source == "UCSC":
            collapse_args = [
                format_ucsc(annotation_path=inputs.annotation.output.annot.path),
                "collapsed_annotation.gtf",
            ]
        else:
            collapse_args = [
                inputs.annotation.output.annot.path,
                "collapsed_annotation.gtf",
            ]

        if "--stranded" in args:
            collapse_args.append("--collapse_only")

        # Collapsing the annotation with collapse_annotation.py script (included in the qc docker image)
        return_code, _, _ = Cmd["collapse_annotation.py"][collapse_args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("Collapse of GTF file failed.")

        self.progress(0.5)

        args.extend(
            [
                "collapsed_annotation.gtf",
                inputs.alignment.output.bam.path,
                "rnaseqc_output",
            ]
        )

        return_code, _, _ = Cmd["rnaseqc"][args] & TEE(retcode=None)
        if return_code:
            self.error("QC analysis failed.")

        outputs.metrics = f"rnaseqc_output/{bam_filename}.metrics.tsv"
