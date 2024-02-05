"""Normalization of STAR quantification results."""

import gzip
import io
import json
import shutil
from pathlib import Path

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    IntegerField,
    JsonField,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio

STRANDEDNESS_CODES = {
    "IU": 0,
    "U": 0,
    "non_specific": 0,
    "ISF": 1,
    "OSF": 1,
    "SF": 1,
    "forward": 1,
    "ISR": 2,
    "OSR": 2,
    "SR": 2,
    "reverse": 2,
}


def prepare_gene_counts(infile, outfile, summary, strandedness):
    """Extract gene counts from STAR input."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        names=["Geneid", 0, 1, 2],
        index_col="Geneid",
        dtype={"Geneid": str, 0: int, 1: int, 2: int},
    )
    # Raw counts for genes
    gene_rc_df = exp.iloc[4:][[strandedness]]
    gene_rc_df.to_csv(
        outfile,
        index_label="FEATURE_ID",
        header=["EXPRESSION"],
        sep="\t",
    )

    assigned_reads = gene_rc_df.sum()
    assigned_reads = int(assigned_reads.values)
    summary_df = exp.iloc[:4][[strandedness]]
    summary_df.loc["N_assigned"] = assigned_reads

    summary_df.to_csv(
        summary,
        index_label="Status",
        header=["Read count"],
        sep="\t",
    )
    return exp.iloc[4:].index.to_list()


def rename_columns_and_compress(infile, outfile):
    """Rename columns and compress the expression files."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        usecols=["FEATURE_ID", "EXPRESSION"],
        index_col="FEATURE_ID",
        squeeze=True,
        float_precision="round_trip",
    )
    exp.to_csv(
        outfile, index_label="Gene", header=["Expression"], sep="\t", compression="gzip"
    )


def compress_outputs(input_file, output_file):
    """Compress outputs."""
    with open(file=input_file, mode="rb") as f_in:
        with gzip.open(filename=output_file, mode="wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def exp_to_df(exp_file, exp_type):
    """Prepare expression file for gene sets merging."""
    with open(exp_file) as exp:
        df = pd.read_csv(exp, sep="\t", float_precision="round_trip")
        df.rename(
            index=str,
            columns={
                "EXPRESSION": exp_type,
            },
            inplace=True,
        )
        # Cast FEATURE_ID column to string
        df["FEATURE_ID"] = df["FEATURE_ID"].astype("str")

    return df


def prepare_expression_set(rc, tpm, cpm, feature_dict, outfile_name):
    """Prepare expression set output data."""
    rc_exp = exp_to_df(rc, "RAW_COUNT")
    tpm_exp = exp_to_df(tpm, "TPM")
    cpm_exp = exp_to_df(cpm, "CPM")
    rc_exp["GENE_SYMBOL"] = rc_exp["FEATURE_ID"].map(feature_dict)
    input_features = rc_exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(rc_exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    # Merge with normalized expression values
    exp_set = rc_exp.merge(tpm_exp, on="FEATURE_ID")
    exp_set = exp_set.merge(cpm_exp, on="FEATURE_ID")
    # Reorder columns
    columns = ["FEATURE_ID", "GENE_SYMBOL", "RAW_COUNT", "TPM", "CPM"]
    exp_set = exp_set[columns]
    # Replace NaN values with empty string
    exp_set.fillna("", inplace=True)

    # Write to file
    exp_set.to_csv(
        outfile_name + ".txt.gz",
        header=True,
        index=False,
        sep="\t",
        compression="gzip",
    )

    # Write to JSON
    df_dict = exp_set.set_index("FEATURE_ID").to_dict(orient="index")
    with open(outfile_name + ".json", "w") as f:
        json.dump({"genes": df_dict}, f, allow_nan=False)


def expression_to_storage(rc_input, rc_output):
    """Convert expressions file to JSON format."""

    def isfloat(value):
        """Check if value is float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    with io.TextIOWrapper(io.BufferedReader(gzip.open(rc_input))) as f:
        # Split lines by tabs
        # Ignore lines without a number in second column
        # Build a dictionary of gene-expression pairs
        exp = {
            "genes": {
                gene_exp[0]: float(gene_exp[1])
                for gene_exp in (l.split("\t") for l in f)
                if len(gene_exp) == 2 and isfloat(gene_exp[1])
            }
        }

    with open(file=rc_output, mode="wt") as f:
        json.dump(exp, f)

    return rc_output


class NormalizeSTARGeneQuantification(ProcessBio):
    """Normalize STAR quantification results.

    This process is based on the output of STAR aligner 'gene counts'.
    Strandedness is detected with Salmon or it can be specified manually.
    Finally, for normalization of gene counts TPM and CPM are used.
    """

    slug = "star-quantification"
    name = "STAR gene quantification"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.2.0",
            },
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
            "network": True,
        },
    }
    data_name = "{{ aligned_reads|name|default('?') }}"
    version = "1.2.0"
    process_type = "data:expression:star"
    category = "Quantify"
    entity = {
        "type": "sample",
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        aligned_reads = DataField(
            "alignment:bam:star",
            label="Aligned reads",
            description="Make sure aligned object from STAR also "
            "includes gene counts otherwise the process will fail.",
        )

        annotation = DataField(
            "annotation",
            label="Annotation",
            description="GTF and GFF3 annotation formats are supported.",
        )

        normalization_type = StringField(
            label="Normalization type",
            default="TPM",
            choices=[
                ("TPM", "TPM"),
                ("CPM", "CPM"),
            ],
            description="The expression normalization type.",
        )

        assay_type = StringField(
            label="Assay type",
            default="non_specific",
            choices=[
                ("non_specific", "Strand non-specific"),
                ("forward", "Strand-specific forward"),
                ("reverse", "Strand-specific reverse"),
                ("auto", "Detect automatically"),
            ],
            description="Indicate if strand-specific read counting should be performed. "
            "For paired-end reads, strand of the first read is taken as the strand "
            "of the whole fragment. FLAG field is used to tell if a read is "
            "first or second read in a pair. Automated strand detection is enabled "
            "using the [Salmon](https://salmon.readthedocs.io/en/latest/library_type.html) "
            "tool's build-in functionality. To use this option, cDNA (transcriptome) "
            "index file created using the Salmon indexing tool must be provided",
        )

        cdna_index = DataField(
            "index:salmon",
            label="Salmon index file",
            required=False,
            hidden="assay_type != 'auto'",
            description="Transcriptome index file created using the Salmon indexing tool. "
            "cDNA (transcriptome) sequences used for index file creation must be "
            "derived from the same species as the input sequencing reads to "
            "obtain the reliable analysis results.",
        )

        n_reads = IntegerField(
            label="Number of reads in subsampled alignment file",
            default=5000000,
            hidden="assay_type != 'auto'",
            description="Alignment (.bam) file subsample size to detect "
            "strandedness. Increase the number of reads to make automatic "
            "detection more reliable. Decrease the number of reads to "
            "make automatic detection run faster.",
        )

    class Output:
        """Output fields."""

        rc = FileField(label="Read counts")
        tpm = FileField(label="TPM")
        cpm = FileField(label="CPM")
        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        counts_summary = FileField(label="Counts summary")
        strandedness_report = FileField(
            label="Strandedness report file",
            required=False,
        )
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""

        if not inputs.aligned_reads.output.gene_counts:
            self.error(
                "Aligned reads should contain gene count information, but do not."
            )

        if inputs.aligned_reads.output.species != inputs.annotation.output.species:
            self.error(
                f"Species of aligned reads {inputs.aligned_reads.output.species} "
                f"and annotation {inputs.annotation.output.species} do not match. Please provide "
                "aligned reads and annotation with the same species."
            )

        if inputs.aligned_reads.output.build != inputs.annotation.output.build:
            self.error(
                f"Builds of aligned reads {inputs.aligned_reads.output.species} "
                f"and annotation {inputs.annotation.output.species} do not match. Please provide "
                "aligned reads and annotation with the same build."
            )

        if inputs.assay_type == "auto" and not inputs.cdna_index:
            self.error(
                "cDNA sequence index must be provided to automatically detect strandedness."
            )

        if (
            inputs.cdna_index
            and inputs.aligned_reads.output.species != inputs.cdna_index.output.species
        ):
            self.error(
                f"Species of aligned reads {inputs.aligned_reads.output.species} "
                f"and cDNA index {inputs.annotation.output.species} do not match. Please provide "
                "aligned reads and cDNA index with the same species."
            )

        bam_file = Path(inputs.aligned_reads.output.bam.path)

        # Set output file names
        assert bam_file.name.endswith(".bam")
        name = bam_file.name[:-4]

        # check if aligned reads are single or paired-end
        paired_end = True
        if int(Cmd["samtools"]["view"]["-c", "-f", "1", bam_file]().strip()) == 0:
            paired_end = False

        # set strandedness
        if inputs.assay_type == "auto":
            all_reads = int(Cmd["samtools"]["view"]["-c", bam_file]().strip())
            sampling_rate = min(inputs.n_reads / all_reads, 1)
            # subsample the BAM file
            if sampling_rate < 1:
                strand_check_bam = "subsampled_sorted.bam"
                (
                    Cmd["samtools"]["view"][
                        f"-@ {self.requirements.resources.cores}",
                        "-h",
                        f"-s {sampling_rate}",
                        bam_file,
                    ]
                    | Cmd["samtools"]["sort"][
                        f"-@ {self.requirements.resources.cores}", "-n", "-"
                    ]
                    > strand_check_bam
                )()
            else:
                strand_check_bam = "sorted.bam"
                sort_args = [
                    f"-@ {self.requirements.resources.cores}",
                    "-n",
                    "-o",
                    strand_check_bam,
                ]
                return_code, _, _ = Cmd["samtools"]["sort"][sort_args][bam_file] & TEE(
                    retcode=None
                )
                if return_code:
                    self.error("Error while running Samtools sort.")

            # Consider only proper paired-end reads for strandedness detection (-0, -s to /dev/null).
            # Failure to do so will result in improper strandedness detection.
            fastq_args = [f"-@ {self.requirements.resources.cores}", "-N"]
            if paired_end:
                reads_input = ["-1", "mate1.fastq", "-2", "mate2.fastq"]
                fastq_args.extend(["-0", "/dev/null", "-s", "/dev/null"])
            else:
                reads_input = ["-0", "reads.fastq"]

            fastq_args.extend(reads_input)

            return_code, _, _ = Cmd["samtools"]["fastq"][fastq_args][
                strand_check_bam
            ] & TEE(retcode=None)
            if return_code:
                self.error("Samtools fastq command failed.")

            salmon_out_folder = "salmon_output"

            # Run Salmon Quant
            salmon_args = [
                "-i",
                inputs.cdna_index.output.index.path,
                "-l",
                "A",
                reads_input if paired_end else ["-r", "reads.fastq"],
                "-o",
                salmon_out_folder,
                "-p",
                self.requirements.resources.cores,
                "--minAssignedFrags",
                1,
            ]
            return_code, _, _ = Cmd["salmon"]["quant"][salmon_args] & TEE(retcode=None)
            if return_code:
                self.error("Error while running Salmon Quant.")

            # Extract the strandedness code from the JSON report produced by the Salmon tool
            lib_type_report = f"{salmon_out_folder}/lib_format_counts.json"
            outputs.strandedness_report = lib_type_report
            strand_code = json.load(open(lib_type_report)).get("expected_format", "")

            if strand_code:
                try:
                    strandedness = STRANDEDNESS_CODES[strand_code]
                except KeyError:
                    self.error(
                        f"Unsupported strand code detected: {strand_code} "
                        "Please re-run analysis in user-selected strandedness mode "
                        "or try increasing the subsample size."
                    )
            else:
                self.error(
                    "Automated detection of strandedness failed. "
                    "Re-run analysis in user-selected strandedness mode."
                )
        else:
            strandedness = STRANDEDNESS_CODES[inputs.assay_type]

        raw_counts = "rc.txt"
        tpm = "tpm.txt"
        cpm = "cpm.txt"

        # prepare_gene_counts() has a side effect of creating csv files
        # used in 'rnanorm'.
        features = prepare_gene_counts(
            infile=inputs.aligned_reads.output.gene_counts.path,
            outfile=raw_counts,
            summary="summary.txt",
            strandedness=strandedness,
        )

        annotation_file = inputs.annotation.output.annot.path
        if (
            inputs.annotation.output.source == "UCSC"
            and inputs.annotation.type.startswith("data:annotation:gtf")
        ):
            with open(annotation_file, "r") as infile:
                filedata = infile.read()

            # Replace the missing gene_ids
            annot_data = filedata.replace('gene_id "";', 'gene_id "unknown";')

            # Write the output file
            annotation_file = "annotation_modified.gtf"
            with open(annotation_file, "w") as outfile:
                outfile.write(annot_data)

        # Normalize counts
        rnanorm_args = [
            raw_counts,
            "--tpm-output",
            tpm,
            "--cpm-output",
            cpm,
            "--annotation",
            annotation_file,
        ]
        return_code, _, _ = Cmd["rnanorm"][rnanorm_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while normalizing counts using rnanorm.")

        feature_filters = {
            "source": inputs.annotation.output.source,
            "species": inputs.aligned_reads.output.species,
            "feature_id__in": features,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        prepare_expression_set(
            rc=raw_counts,
            tpm=tpm,
            cpm=cpm,
            feature_dict=feature_ids_to_names,
            outfile_name=f"{name}_expressions",
        )

        # rename and compress the expression files
        rename_columns_and_compress(raw_counts, f"{name}_rc.tab.gz")
        rename_columns_and_compress(tpm, f"{name}_tpm.tab.gz")
        rename_columns_and_compress(cpm, f"{name}_cpm.tab.gz")

        exp_output = f"{name}_{inputs.normalization_type.lower()}.tab.gz"

        # Save the abundance estimates to JSON storage
        json_output = "json.txt"
        expression_to_storage(rc_input=exp_output, rc_output=json_output)

        # Save the outputs
        outputs.counts_summary = "summary.txt"
        outputs.rc = f"{name}_rc.tab.gz"
        outputs.tpm = f"{name}_tpm.tab.gz"
        outputs.cpm = f"{name}_cpm.tab.gz"
        outputs.exp = exp_output
        outputs.exp_json = json_output
        outputs.exp_set = f"{name}_expressions.txt.gz"
        outputs.exp_set_json = f"{name}_expressions.json"
        outputs.exp_type = inputs.normalization_type
        outputs.source = inputs.annotation.output.source
        outputs.species = inputs.aligned_reads.output.species
        outputs.build = inputs.aligned_reads.output.build
        outputs.feature_type = "gene"
