"""Run Cuffdiff 2.2."""
import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    JsonField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class Cuffdiff(Process):
    """Run Cuffdiff 2.2 analysis.

    Cuffdiff finds significant changes in transcript expression, splicing, and
    promoter use.  You can use it to find differentially expressed genes and
    transcripts, as well as genes that are being differentially regulated at
    the transcriptional and post-transcriptional level. See
    [here](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/) and
    [here](https://software.broadinstitute.org/cancer/software/genepattern/modules/docs/Cuffdiff/7)
    for more information.
    """

    slug = "cuffdiff"
    name = "Cuffdiff 2.2"
    process_type = "data:differentialexpression:cuffdiff"
    version = "3.2.0"
    category = "Differential Expression"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/rnaseq:4.9.0"}},
        "resources": {"cores": 10, "memory": 8192},
    }
    data_name = "Cuffdiff results"

    class Input:
        """Input fields to process Cuffdiff."""

        case = ListField(
            DataField("cufflinks:cuffquant"),
            label="Case samples",
        )
        control = ListField(
            DataField("cufflinks:cuffquant"),
            label="Control samples",
        )
        labels = ListField(
            StringField(),
            label="Group labels",
            description="Define labels for each sample group.",
            default=["control", "case"],
        )
        annotation = DataField(
            "annotation",
            label="Annotation (GTF/GFF3)",
            description="A transcript annotation file produced by "
            "cufflinks, cuffcompare, or other tool.",
        )
        genome = DataField(
            "seq:nucleotide",
            label="Run bias detection and correction algorithm",
            required=False,
            description="Provide Cufflinks with a multifasta file "
            "(genome file) via this option to instruct it to run a "
            "bias detection and correction algorithm which can "
            "significantly improve accuracy of transcript abundance "
            "estimates.",
        )
        multi_read_correct = BooleanField(
            label="Do initial estimation procedure to more accurately "
            "weight reads with multiple genome mappings",
            default=False,
        )
        create_sets = BooleanField(
            label="Create gene sets",
            description="After calculating differential gene "
            "expressions create gene sets for up-regulated genes, "
            "down-regulated genes and all genes.",
            default=False,
        )
        gene_logfc = FloatField(
            label="Log2 fold change threshold for gene sets",
            description="Genes above Log2FC are considered as "
            "up-regulated and genes below -Log2FC as down-regulated.",
            default=1.0,
            hidden="!create_sets",
        )
        gene_fdr = FloatField(
            label="FDR threshold for gene sets",
            default=0.05,
            hidden="!create_sets",
        )
        fdr = FloatField(
            label="Allowed FDR",
            description="The allowed false discovery rate. The default is 0.05.",
            default=0.05,
        )
        library_type = StringField(
            label="Library type",
            description="In cases where Cufflinks cannot determine the "
            "platform and protocol used to generate input reads, you "
            "can supply this information manually, which will allow "
            "Cufflinks to infer source strand information with certain "
            "protocols. The available options are listed below. For "
            "paired-end data, we currently only support protocols "
            "where reads point towards each other: fr-unstranded - "
            "Reads from the left-most end of the fragment (in "
            "transcript coordinates) map to the transcript strand and "
            "the right-most end maps to the opposite strand; "
            "fr-firststrand - Same as above except we enforce the rule "
            "that the right-most end of the fragment (in transcript "
            "coordinates) is the first sequenced (or only sequenced "
            "for single-end reads). Equivalently, it is assumed that "
            "only the strand generated during first strand synthesis "
            "is sequenced; fr-secondstrand - Same as above except we "
            "enforce the rule that the left-most end of the fragment "
            "(in transcript coordinates) is the first sequenced (or "
            "only sequenced for single-end reads). Equivalently, it is "
            "assumed that only the strand generated during second "
            "strand synthesis is sequenced.",
            default="fr-unstranded",
            choices=[
                ("fr-unstranded", "fr-unstranded"),
                ("fr-firststrand", "fr-firststrand"),
                ("fr-secondstrand", "fr-secondstrand"),
            ],
        )
        library_normalization = StringField(
            label="Library normalization method",
            description="You can control how library sizes (i.e. "
            "sequencing depths) are normalized in Cufflinks and "
            "Cuffdiff. Cuffdiff has several methods that require "
            "multiple libraries in order to work. Library "
            "normalization methods supported by Cufflinks work on one "
            "library at a time.",
            default="geometric",
            choices=[
                ("geometric", "geometric"),
                ("classic-fpkm", "classic-fpkm"),
                ("quartile", "quartile"),
            ],
        )
        dispersion_method = StringField(
            label="Dispersion method",
            description=" Cuffdiff works by modeling the variance in "
            "fragment counts across replicates as a function of the "
            "mean fragment count across replicates. Strictly speaking, "
            "models a quantitity called dispersion - the variance "
            "present in a group of samples beyond what is expected "
            "from a simple Poisson model of RNA_Seq. You can control "
            "how Cuffdiff constructs its model of dispersion in locus "
            "fragment counts. Each condition that has replicates can "
            "receive its own model, or Cuffdiff can use a global model "
            "for all conditions. All of these policies are identical "
            "to those used by DESeq (Anders and Huber, Genome Biology, "
            "2010).",
            default="pooled",
            choices=[
                ("pooled", "pooled"),
                ("per-condition", "per-condition"),
                ("blind", "blind"),
                ("poisson", "poisson"),
            ],
        )

    class Output:
        """Output fields of the process Cuffdiff."""

        raw = FileField("Differential expression")
        de_json = JsonField(label="Results table (JSON)")
        de_file = FileField(label="Results table (file)")
        transcript_diff_exp = FileField(
            label="Differential expression (transcript level)"
        )
        tss_group_diff_exp = FileField(
            label="Differential expression (primary transcript)"
        )
        cds_diff_exp = FileField(label="Differential expression (coding sequence)")
        cuffdiff_output = FileField(label="Cuffdiff output")
        source = StringField(label="Gene ID database")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""

        cuffquants = inputs.case + inputs.control

        for c in cuffquants:
            if c.source != cuffquants[0].source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{c.source} and {cuffquants[0].source}."
                )
            if c.species != cuffquants[0].species:
                self.error(
                    "Input samples are of different Species: "
                    f"{c.species} and {cuffquants[0].species}."
                )
            if c.build != cuffquants[0].build:
                self.error(
                    "Input samples are of different Panel types: "
                    f"{c.build} and {cuffquants[0].build}."
                )

        for case in inputs.case:
            if case in inputs.control:
                self.error(
                    "Case and Control groups must contain unique "
                    f"samples. Sample {case.sample_name} is in both Case "
                    "and Control group."
                )

        case_paths = ",".join([case.cxb.path for case in inputs.case])
        control_paths = ",".join([control.cxb.path for control in inputs.control])
        labels = ",".join(inputs.labels)

        outputs.source = cuffquants[0].source
        outputs.species = cuffquants[0].species
        outputs.build = cuffquants[0].build
        outputs.feature_type = "gene"

        self.progress(0.1)

        params = [
            "-output-dir",
            "./",
            "-num-threads",
            self.requirements.resources.cores,
            "-labels",
            labels,
            "-FDR",
            inputs.fdr,
            "-library-type",
            inputs.library_type,
            "-library-norm-method",
            inputs.library_normalization,
            "-dispersion-method",
            inputs.dispersion_method,
            "-quiet",
        ]
        if inputs.genome:
            params.extend(["-frag-bias-correct", inputs.genome.fasta.path])
        if inputs.multi_read_correct:
            params.append("-multi-read-correct")

        return_code, _, _ = Cmd["cuffdiff"][params][
            inputs.annotation.annot.path, control_paths, case_paths
        ] & TEE(retcode=None)
        if return_code:
            self.error("Error while computing differential expression with Cuffdiff.")

        self.progress(0.90)

        exp_file = "cuffdiff.tab"
        os.rename("gene_exp.diff", exp_file)

        files_list = [
            "cds.*",
            "isoforms.*",
            "genes.*",
            "tss_groups.*",
            "read_groups.*",
            "promoters.diff",
            "splicing.diff",
            "cds_exp.diff",
            exp_file,
            "isoform_exp.diff",
            "tss_group_exp.diff",
        ]
        zip_file = "cuffdiff_output.zip"
        return_code, _, _ = Cmd["zip"][zip_file][files_list] & TEE(retcode=None)
        if return_code:
            self.error("Error while compressing Cuffdiff files.")

        args = [
            exp_file,
            "de_data.json",
            "de_file.tab.gz",
            "--gene_id",
            "gene_id",
            "--fdr",
            "q_value",
            "--pvalue",
            "p_value",
            "--logfc",
            "log2(fold_change)",
            "--stat",
            "test_stat",
        ]

        return_code, _, _ = Cmd["parse_diffexp.py"][args] & TEE(retcode=None)
        if return_code:
            self.error(f"Error while parsing DGE results.")

        (Cmd["gzip"][exp_file])()

        outputs.raw = f"{exp_file}.gz"
        outputs.de_json = "de_data.json"
        outputs.de_file = "de_file.tab.gz"
        outputs.transcript_diff_exp = "isoform_exp.diff"
        outputs.cds_diff_exp = "cds_exp.diff"
        outputs.tss_group_diff_exp = "tss_group_exp.diff"
        outputs.cuffdiff_output = "cuffdiff_output.zip"

        if inputs.create_sets:
            out_dir = "gene_sets"
            gene_set_args = [
                "--dge_file",
                "de_file.tab.gz",
                "--out_dir",
                out_dir,
                "--analysis_name",
                self.name,
                "--tool",
                "Cuffdiff",
                "--logfc",
                inputs.gene_logfc,
                "--fdr",
                inputs.gene_fdr,
            ]

            return_code, _, _ = Cmd["create_gene_sets.py"][gene_set_args] & TEE(
                retcode=None
            )
            if return_code:
                self.error(f"Error while creating gene sets.")

            for gene_file in sorted(Path(out_dir).glob("*.tab.gz")):
                gene_file.rename(Path() / gene_file.name)
                process_inputs = {
                    "src": str(gene_file.name),
                    "source": cuffquants[0].source,
                    "species": cuffquants[0].species,
                }
                self.run_process("upload-geneset", process_inputs)
