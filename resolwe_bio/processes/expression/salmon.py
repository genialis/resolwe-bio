"""Salmon Quant."""

import gzip
import io
import json
import os

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DirField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    JsonField,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


def parse_transcript_exp(infile, outfile):
    """Parse trancript-level expressions from Salmon output."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        usecols=["Name", "TPM"],
        index_col="Name",
        dtype={
            "Name": str,
            "TPM": float,
        },
        squeeze=True,
    )
    return exp.to_csv(
        outfile,
        index_label="Transcript",
        header=["Expression"],
        sep="\t",
        compression="gzip",
    )


def expression_to_storage(exp_input, exp_output):
    """Convert expressions file to JSON format."""

    def isfloat(value):
        """Check if value is float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    with io.TextIOWrapper(io.BufferedReader(gzip.open(exp_input))) as f:
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

    with open(file=exp_output, mode="wt") as f:
        json.dump(exp, f)

    return exp_output


def rename_cols(infile, outfile, abundance_unit):
    """Rename columns in expression file."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        skip_blank_lines=True,
        usecols=["Gene", "Expression"],
        index_col="Gene",
        dtype={
            "Gene": str,
            "Expression": float,
        },
        squeeze=True,
    )
    return exp.to_csv(
        outfile,
        index_label="FEATURE_ID",
        header=[abundance_unit],
        sep="\t",
    )


def prepare_expression_set(infile, abundance_unit, feature_dict, outfile_name):
    """Prepare expression set output data."""
    exp = pd.read_csv(infile, sep="\t", float_precision="round_trip")
    exp["FEATURE_ID"] = exp["FEATURE_ID"].astype("str")
    exp["GENE_SYMBOL"] = exp["FEATURE_ID"].map(feature_dict)
    input_features = exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    columns = ["FEATURE_ID", "GENE_SYMBOL", abundance_unit]
    exp_set = exp[columns]
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


class SalmonQuant(ProcessBio):
    """Perform mapping-based estimation of transcript abundance from RNA-seq reads.

    Final abundance estimates are summarized to the gene-level using
    [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html).
    """

    slug = "salmon-quant"
    name = "Salmon Quant"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0",
            },
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
            "network": True,
        },
    }
    data_name = "{{ reads|name|default('?') }}"
    version = "2.7.1"
    process_type = "data:expression:salmon"
    category = "Quantify"
    entity = {
        "type": "sample",
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        reads = DataField("reads:fastq", label="Input sample(s)")
        salmon_index = DataField("index:salmon", label="Salmon index")
        annotation = DataField("annotation:gtf", label="GTF annotation")

        class Options:
            """Options."""

            stranded = StringField(
                label="Assay type",
                default="A",
                choices=[
                    ("A", "Detect automatically"),
                    ("U", "Strand non-specific (U)"),
                    ("SF", "Strand-specific forward (SF)"),
                    ("SR", "Strand-specific reverse (SR)"),
                    ("IU", "Strand non-specific (paired-end IU)"),
                    ("ISF", "Strand-specific forward (paired-end ISF)"),
                    ("ISR", "Strand-specific reverse (paired-end (ISR)"),
                ],
            )

            seq_bias = BooleanField(
                label="--seqBias",
                default=False,
                description="Perform sequence-specific bias correction.",
            )

            gc_bias = BooleanField(
                label="--gcBias",
                default=False,
                description="[beta for single-end reads] Perform fragment GC bias correction.",
            )

            discard_orphans_quasi = BooleanField(
                label="--discardOrphansQuasi",
                default=False,
                description="Discard orphan mappings in quasi-mapping mode. "
                "If this flag is passed then only paired "
                "mappings will be considered toward "
                "quantification estimates. The default "
                "behavior is to consider orphan mappings "
                "if no valid paired mappings exist.",
            )

            no_length_correction = BooleanField(
                label="--noLengthCorrection",
                default=False,
                description="[Experimental] Entirely disables "
                "length correction when estimating the "
                "abundance of transcripts. The abundance "
                "estimates are reported in CPM (counts per "
                "million) unit. This option can be used "
                "with protocols where one expects that "
                "fragments derive from their underlying "
                "targets without regard to that target's  "
                "length (e.g. QuantSeq).",
            )

            consensus_slack = FloatField(
                label="--consensusSlack",
                required=False,
                description="The amount of slack allowed in the quasi-mapping "
                "consensus mechanism.  Normally, a transcript must "
                "cover all hits to be considered for mapping.  "
                "If this is set to a fraction, X, greater than 0 "
                "(and in [0,1)), then a transcript can fail "
                "to cover up to (100 * X)% of the hits before it "
                "is discounted as a mapping candidate. The default "
                "value of this option is 0.2 in selective alignment mode "
                "and 0 otherwise.",
            )

            min_score_fraction = FloatField(
                label="--minScoreFraction",
                default=0.65,
                description="The fraction of the optimal possible alignment "
                "score that a mapping must achieve in order to be "
                "considered valid - should be in (0,1]",
            )

            incompat_prior = FloatField(
                label="---incompatPrior",
                default=0,
                description="This option sets the prior probability "
                "that an alignment that disagrees with "
                "the specified library type (--libType) "
                "results from the true fragment origin. "
                "Setting this to 0 specifies that "
                "alignments that disagree with the "
                "library type should be impossible, "
                "while setting it to 1 says that "
                "alignments that disagree with the "
                "library type are no less likely than "
                "those that do.",
            )

            range_factorization_bins = IntegerField(
                label="--rangeFactorizationBins",
                default=4,
                description="Factorizes the likelihood used in "
                "quantification by adopting a new notion "
                "of equivalence classes based on the "
                "conditional probabilities with which "
                "fragments are generated from different "
                "transcripts.  This is a more "
                "fine-grained factorization than the "
                "normal rich equivalence classes. The "
                "default value (4) corresponds to the "
                "default used in Zakeri et al. 2017 "
                "and larger values imply a more "
                "fine-grained factorization. If range "
                "factorization is enabled, a common "
                "value to select for this parameter is "
                "4. A value of 0 signifies the use of "
                "basic rich equivalence classes.",
            )

            min_assigned_frag = IntegerField(
                label="--minAssignedFrags",
                default=10,
                description="The minimum number of fragments that "
                "must be assigned to the transcriptome "
                "for quantification to proceed.",
            )
            num_bootstraps = IntegerField(
                label="--numBootstraps",
                description="Salmon has the ability to optionally "
                "compute bootstrapped abundance estimates. This is "
                "done by resampling (with replacement) from the counts "
                "assigned to the fragment equivalence classes, and then "
                "re-running the optimization procedure, either the EM or VBEM, "
                "for each such sample. The values of these different bootstraps "
                "allows us to assess technical variance in the main abundance "
                "estimates we produce. Such estimates can be useful for downstream "
                "(e.g. differential expression) tools that can make use of such "
                "uncertainty estimates. This option takes a positive integer that "
                "dictates the number of bootstrap samples to compute. The more samples "
                "computed, the better the estimates of varaiance, but the more "
                "computation (and time) required.",
                disabled="options.num_gibbs_samples",
                required=False,
            )
            num_gibbs_samples = IntegerField(
                label="--numGibbsSamples",
                description="Just as with the bootstrap procedure above, this option "
                "produces samples that allow us to estimate the variance in abundance "
                "estimates. However, in this case the samples are generated using posterior "
                "Gibbs sampling over the fragment equivalence classes rather than "
                "bootstrapping. We are currently analyzing these different approaches to "
                "assess the potential trade-offs in time / accuracy. The --numBootstraps "
                "and --numGibbsSamples options are mutually exclusive (i.e. in a given run, "
                "you must set at most one of these options to a positive integer.)",
                disabled="options.num_bootstraps",
                required=False,
            )

        options = GroupField(Options, label="Options")

    class Output:
        """Output fields."""

        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        rc = FileField(label="Gene-level estimated counts")
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        variance = FileField(label="Variance of inferential replicates", required=False)
        quant = FileField(label="Salmon quant file")
        transcripts = FileField(label="Transcript-level expressions")
        salmon_output = DirField(label="Salmon output")
        txdb = FileField(label="Transcript to gene mapping")
        strandedness = StringField(label="Strandedness code")
        strandedness_report = FileField(label="Strandedness report file")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if inputs.salmon_index.output.species != inputs.annotation.output.species:
            self.error(
                "Salmon index file species ({}) must match GTF annotation "
                "file species ({})".format(
                    inputs.salmon_index.output.species, inputs.annotation.output.species
                )
            )

        if inputs.salmon_index.output.build != inputs.annotation.output.build:
            self.error(
                "Salmon index file build ({}) must match GTF annotation "
                "file build ({})".format(
                    inputs.salmon_index.output.build, inputs.annotation.output.build
                )
            )

        if inputs.salmon_index.output.source != inputs.annotation.output.source:
            self.error(
                "Salmon index file source ({}) must match GTF annotation "
                "file source ({})".format(
                    inputs.salmon_index.output.source, inputs.annotation.output.source
                )
            )

        if inputs.options.no_length_correction:
            abundance_unit = "CPM"
            output_suffix = "_cpm.txt.gz"
        else:
            abundance_unit = "TPM"
            output_suffix = "_tpm.txt.gz"

        args = [
            "-i",
            inputs.salmon_index.output.index.path,
            "-l",
            inputs.options.stranded,
            "--incompatPrior",
            inputs.options.incompat_prior,
            "--minAssignedFrags",
            inputs.options.min_assigned_frag,
            "--rangeFactorizationBins",
            inputs.options.range_factorization_bins,
            "-p",
            self.requirements.resources.cores,
            "-o",
            "salmon_output",
        ]

        # Prepare .FASTQ file inputs based on the reads input type
        if inputs.reads.type.startswith("data:reads:fastq:single:"):
            args.extend(["-r"] + [lane.path for lane in inputs.reads.output.fastq])
        else:
            args.extend(["-1"] + [lane.path for lane in inputs.reads.output.fastq])
            args.extend(["-2"] + [lane.path for lane in inputs.reads.output.fastq2])

        # Prepare optional inputs
        if inputs.options.seq_bias and not inputs.options.no_length_correction:
            args.append("--seqBias")
        elif inputs.options.seq_bias and inputs.options.no_length_correction:
            self.warning(
                "Since bias correction relies on modifying effective lengths, "
                "you cannot enable bias correction simultaneously with the "
                "--noLengthCorrection option. Skipping --seqBias option."
            )

        if inputs.options.gc_bias and not inputs.options.no_length_correction:
            args.append("--gcBias")
        elif inputs.options.gc_bias and inputs.options.no_length_correction:
            self.warning(
                "Since bias correction relies on modifying effective lengths, "
                "you cannot enable bias correction simultaneously with the "
                "--noLengthCorrection option. Skipping --gcBias option."
            )

        if inputs.options.discard_orphans_quasi:
            args.append("--discardOrphansQuasi")

        if inputs.options.no_length_correction:
            args.append("--noLengthCorrection")

        if inputs.options.min_score_fraction > 0:
            args.extend(["--minScoreFraction", inputs.options.min_score_fraction])

        if inputs.options.consensus_slack is not None:
            args.extend(["--consensusSlack", inputs.options.consensus_slack])

        if inputs.options.num_bootstraps:
            args.extend(["--numBootstraps", inputs.options.num_bootstraps])
        elif inputs.options.num_gibbs_samples:
            args.extend(["--numGibbsSamples", inputs.options.num_gibbs_samples])

        # Run Salmon Quant
        return_code, _, _ = Cmd["salmon"]["quant"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error while running Salmon Quant.")

        # Use tximport to produce gene-level TPM values
        reads_basename = os.path.basename(inputs.reads.output.fastq[0].path)
        assert reads_basename.endswith(".fastq.gz")
        reads_name = reads_basename[:-9]
        annot_basename = os.path.basename(inputs.annotation.output.annot.path)
        assert annot_basename.endswith(".gtf")
        annot_name = annot_basename[:-4]
        tx2gene = "tx2gene_{}.txt".format(annot_name)
        counts = f"{reads_name}_counts.txt"
        counts_gz = counts + ".gz"
        if os.path.exists("salmon_output/quant.sf"):
            tximport_args = [
                "salmon_output/quant.sf",
                inputs.annotation.output.annot.path,
                reads_name,
                counts,
                tx2gene,
            ]
            # Strip feature_id version for non-UCSC annotation source type
            # UCSC annotation type (mm10) contains features with dot in gene names
            if inputs.annotation.output.source != "UCSC":
                tximport_args.append("--ignoreTxVersion")

            if inputs.options.num_bootstraps or inputs.options.num_gibbs_samples:
                tximport_args.append("--variance")
                variance = f"variance_{reads_name}"
                variance_gz = variance + ".txt.gz"

            return_code, _, _ = Cmd["tximport_summarize.R"][tximport_args] & TEE(
                retcode=None
            )
            if return_code:
                self.error("Error while running tximport.")
            # Prepare transcript-level expression file
            transcript_out_file = "{}_transcripts{}".format(reads_name, output_suffix)
            parse_transcript_exp("salmon_output/quant.sf", transcript_out_file)
        else:
            self.error("Salmon Quant results file quant.sf does not exists.")

        # Zip the gene-level abundance estimates
        (Cmd["gzip"]["-c", reads_name] > reads_name + output_suffix)()

        # Zip the gene-level count estimates
        (Cmd["gzip"]["-c", counts] > counts_gz)()

        if inputs.options.num_bootstraps or inputs.options.num_gibbs_samples:
            (Cmd["gzip"]["-c", variance] > variance_gz)()
            outputs.variance = variance_gz

        # Save the abundance estimates to JSON storage
        json_output = "json.txt"
        expression_to_storage(
            exp_input=(reads_name + output_suffix), exp_output=json_output
        )

        # Rename columns of the expression file
        reads_name_renamed = f"{reads_name}_renamed"
        rename_cols(
            infile=reads_name, outfile=reads_name_renamed, abundance_unit=abundance_unit
        )

        # Prepare the expression set outputs
        feature_ids = pd.read_csv(
            reads_name_renamed, sep="\t", index_col="FEATURE_ID"
        ).index.tolist()

        feature_filters = {
            "source": inputs.annotation.output.source,
            "species": inputs.annotation.output.species,
            "feature_id__in": feature_ids,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        prepare_expression_set(
            infile=reads_name_renamed,
            abundance_unit=abundance_unit,
            feature_dict=feature_ids_to_names,
            outfile_name=f"{reads_name}_expressions",
        )

        Cmd["ln"]["-s", "salmon_output/quant.sf", reads_name + ".sf"]()

        lib_type_report = "salmon_output/lib_format_counts.json"
        strandedness = json.load(open(lib_type_report)).get("expected_format", "")

        # Save all the outputs
        outputs.salmon_output = "salmon_output"
        outputs.quant = reads_name + ".sf"
        outputs.transcripts = transcript_out_file
        outputs.txdb = tx2gene
        outputs.rc = counts_gz
        outputs.exp = reads_name + output_suffix
        outputs.exp_json = json_output
        outputs.exp_set = reads_name + "_expressions.txt.gz"
        outputs.exp_set_json = reads_name + "_expressions.json"
        outputs.strandedness = strandedness
        outputs.strandedness_report = lib_type_report
        outputs.exp_type = abundance_unit
        outputs.feature_type = "gene"
        outputs.source = inputs.salmon_index.output.source
        outputs.species = inputs.salmon_index.output.species
        outputs.build = inputs.salmon_index.output.build
