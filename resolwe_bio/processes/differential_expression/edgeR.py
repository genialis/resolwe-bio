"""Run EdgeR analysis."""

from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    IntegerField,
    JsonField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class EdgeR(Process):
    """Run EdgeR analysis.

    Empirical Analysis of Digital Gene Expression Data in R (edgeR).
    Differential expression analysis of RNA-seq expression profiles with
    biological replication. Implements a range of statistical methodology
    based on the negative binomial distributions, including empirical Bayes
    estimation, exact tests, generalized linear models and quasi-likelihood
    tests. As well as RNA-seq, it be applied to differential signal analysis
    of other types of genomic data that produce counts, including ChIP-seq,
    Bisulfite-seq, SAGE and CAGE. See
    [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
    for more information.
    """

    slug = "differentialexpression-edger"
    name = "edgeR"
    process_type = "data:differentialexpression:edger"
    version = "1.7.0"
    category = "Differential Expression"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {"cores": 1, "memory": 8192},
    }
    data_name = "Differential expression (case vs. control)"

    class Input:
        """Input fields to process EdgeR."""

        case = ListField(
            DataField("expression"),
            label="Case",
            description="Case samples (replicates)",
        )
        control = ListField(
            DataField("expression"),
            label="Control",
            description="Control samples (replicates)",
        )

        count_filter = IntegerField(
            label="Raw counts filtering threshold",
            default=10,
            description="Filter genes in the expression matrix input. "
            "Remove genes where the number of counts in all samples is "
            "below the threshold.",
        )
        create_sets = BooleanField(
            label="Create gene sets",
            description="After calculating differential gene "
            "expressions create gene sets for up-regulated genes, "
            "down-regulated genes and all genes.",
            default=False,
        )
        logfc = FloatField(
            label="Log2 fold change threshold for gene sets",
            description="Genes above Log2FC are considered as "
            "up-regulated and genes below -Log2FC as down-regulated.",
            default=1.0,
            hidden="!create_sets",
        )
        fdr = FloatField(
            label="FDR threshold for gene sets",
            default=0.05,
            hidden="!create_sets",
        )

    class Output:
        """Output fields of the process EdgeR."""

        raw = FileField("Differential expression")
        de_json = JsonField(label="Results table (JSON)")
        de_file = FileField(label="Results table (file)")
        source = StringField(label="Gene ID database")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if any(
            e.type == "data:expression:microarray:"
            for e in inputs.case + inputs.control
        ):
            self.error("Microarray expressions are not supported.")

        for t in inputs.case:
            if t in inputs.control:
                self.error(
                    "Case and Control groups must contain unique "
                    f"samples. Sample {t.sample_name} is in both Case "
                    "and Control group."
                )

        if len(inputs.case) < 2 or len(inputs.control) < 2:
            self.error(
                "Error in calculating edgeR dispersion, please provide more samples"
            )

        try:
            case_paths = [c.output.rc.path for c in inputs.case]
            control_paths = [c.output.rc.path for c in inputs.control]
        except AttributeError:
            self.error("Read counts are required when using edgeR")

        conditions = ["case"] * len(case_paths) + ["control"] * len(control_paths)

        expressions = inputs.case + inputs.control

        for exp in expressions:
            if exp.output.source != expressions[0].output.source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{exp.output.source} and {expressions[0].output.source}."
                )
            if exp.output.species != expressions[0].output.species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.output.species} and {expressions[0].output.species}."
                )
            if exp.output.build != expressions[0].output.build:
                self.error(
                    "Input samples are of different Panel types: "
                    f"{exp.output.build} and {expressions[0].output.build}."
                )
            if exp.output.feature_type != expressions[0].output.feature_type:
                self.error(
                    "Input samples are of different Feature type: "
                    f"{exp.output.feature_type} and {expressions[0].output.feature_type}."
                )

        for case in inputs.case:
            if case in inputs.control:
                self.error(
                    "Case and Control groups must contain unique "
                    f"samples. Sample {case.sample_name} is in both Case "
                    "and Control group."
                )

        self.progress(0.1)

        sample_files = case_paths + control_paths
        merge_args = [
            sample_files,
            "--experiments",
            sample_files,
            "--intersection",
            "--out",
            "counts.tab",
        ]
        return_code, _, _ = Cmd["expressionmerge.py"][merge_args] & TEE(retcode=None)
        if return_code:
            self.error("Error merging read counts.")

        filter_args = [
            "-counts",
            "counts.tab",
            "-filter",
            inputs.count_filter,
            "-out",
            "counts_filtered.tab",
        ]
        return_code, _, _ = Cmd["diffexp_filtering.R"][filter_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while filtering read counts.")

        args = [
            "counts_filtered.tab",
            "--sampleConditions",
            conditions,
        ]
        return_code, _, _ = Cmd["run_edger.R"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error computing differential expression (edgeR).")

        self.progress(0.95)

        edger_output = "diffexp_edgeR.tab"
        args = [
            edger_output,
            "de_data.json",
            "de_file.tab.gz",
            "--gene_id",
            "gene_id",
            "--fdr",
            "FDR",
            "--pvalue",
            "PValue",
            "--logfc",
            "logFC",
        ]

        return_code, _, _ = Cmd["parse_diffexp.py"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error while parsing DGE results.")

        (Cmd["gzip"][edger_output])()

        outputs.raw = f"{edger_output}.gz"
        outputs.de_json = "de_data.json"
        outputs.de_file = "de_file.tab.gz"
        outputs.source = expressions[0].output.source
        outputs.species = expressions[0].output.species
        outputs.build = expressions[0].output.build
        outputs.feature_type = expressions[0].output.feature_type

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
                "EdgeR",
                "--logfc",
                inputs.logfc,
                "--fdr",
                inputs.fdr,
            ]

            return_code, _, _ = Cmd["create_gene_sets.py"][gene_set_args] & TEE(
                retcode=None
            )
            if return_code:
                self.error("Error while creating gene sets.")

            for gene_file in sorted(Path(out_dir).glob("*.tab.gz")):
                gene_file.rename(Path() / gene_file.name)
                process_inputs = {
                    "src": str(gene_file.name),
                    "source": expressions[0].output.source,
                    "species": expressions[0].output.species,
                }
                self.run_process("upload-geneset", process_inputs)
