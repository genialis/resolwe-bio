"""Run DESeq analysis."""
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    JsonField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class Deseq(Process):
    """Run DESeq2 analysis.

    The DESeq2 package estimates variance-mean dependence in count data
    from high-throughput sequencing assays and tests for differential
    expression based on a model using the negative binomial
    distribution. See
    [here](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)
    and [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
    for more information.
    """

    slug = "differentialexpression-deseq2"
    name = "DESeq2"
    process_type = "data:differentialexpression:deseq2"
    version = "3.3.0"
    category = "Differential Expression"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.9.0"}
        },
        "resources": {"cores": 1, "memory": 8192},
    }
    data_name = "Differential expression (case vs. control)"

    class Input:
        """Input fields to process Deseq."""

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

        class Options:
            """Options."""

            beta_prior = BooleanField(
                label="Beta prior",
                default=False,
                description="Whether or not to put a zero-mean normal prior "
                "on the non-intercept coefficients.",
            )

        class FilterOptions:
            """Filtering options."""

            count = BooleanField(
                label="Filter genes based on expression count",
                default=True,
            )
            min_count_sum = IntegerField(
                label="Minimum gene expression count summed over all samples",
                default=10,
                description="Filter genes in the expression matrix input. "
                "Remove genes where the expression count sum over all samples "
                "is below the threshold.",
                hidden="!filter_options.count",
            )
            cook = BooleanField(
                label="Filter genes based on Cook's distance",
                default=False,
            )
            cooks_cutoff = FloatField(
                label="Threshold on Cook's distance",
                required=False,
                description="If one or more samples have Cook's distance "
                "larger than the threshold set here, the p-value for the row "
                "is set to NA. If left empty, the default threshold of 0.99 "
                "quantile of the F(p, m-p) distribution is used, where p is "
                "the number of coefficients being fitted and m is the number "
                "of samples. This test excludes Cook's distance of samples "
                "belonging to experimental groups with only two samples.",
                hidden="!filter_options.cook",
            )
            independent = BooleanField(
                label="Apply independent gene filtering",
                default=False,
            )
            alpha = FloatField(
                label="Significance cut-off used for optimizing independent "
                "gene filtering",
                default=0.1,
                description="The value should be set to adjusted p-value "
                "cut-off (FDR).",
                hidden="!filter_options.independent",
            )

        options = GroupField(Options, label="Gene filtering options")
        filter_options = GroupField(
            FilterOptions, label="Differential expression analysis options"
        )

    class Output:
        """Output fields of the process Deseq."""

        raw = FileField("Differential expression")
        de_json = JsonField(label="Results table (JSON)")
        de_file = FileField(label="Results table (file)")
        count_matrix = FileField(label="Count matrix")
        source = StringField(label="Gene ID database")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""

        expressions = inputs.case + inputs.control

        if any(e.type == "data:expression:microarray:" for e in expressions):
            self.error("Microarray expressions are not supported.")

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
                    f"{exp.build} and {expressions[0].build}."
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

        if all(e.type == "data:expression:nanostring:" for e in expressions):
            params = [
                "--cases",
                [e.output.exp.path for e in inputs.case],
                "--controls",
                [e.output.exp.path for e in inputs.control],
                "--format",
                "nanostring",
            ]

        elif all(e.type == "data:expression:rsem:" for e in expressions):
            params = [
                "--cases",
                [e.output.genes.path for e in inputs.case],
                "--controls",
                [e.output.genes.path for e in inputs.control],
                "--format",
                "rsem",
            ]

        elif all(e.type == "data:expression:salmon:" for e in expressions):
            params = [
                "--cases",
                [e.output.quant.path for e in inputs.case],
                "--controls",
                [e.output.quant.path for e in inputs.control],
                "--format",
                "salmon",
                "--tx2gene",
                inputs.case[0].output.txdb.path,
            ]

        else:
            if not all(hasattr(e.output.rc, "path") for e in expressions):
                self.error("Read counts are required when using DESeq2.")
            params = [
                "--cases",
                [e.output.rc.path for e in inputs.case],
                "--controls",
                [e.output.rc.path for e in inputs.control],
            ]

        if inputs.options.beta_prior:
            params.append("--beta-prior")
        if inputs.filter_options.count:
            params.extend(["--min-count-sum", inputs.filter_options.min_count_sum])
        if inputs.filter_options.cook:
            params.extend(["--cooks-cutoff", inputs.filter_options.cooks_cutoff])
        if inputs.filter_options.independent:
            params.extend(["--alpha", inputs.filter_options.alpha])

        return_code, _, _ = Cmd["deseq.R"][params] & TEE(retcode=None)
        self.progress(0.95)

        deseq_output = "diffexp_deseq2.tab"
        args = [
            deseq_output,
            "de_data.json",
            "de_file.tab.gz",
            "--gene_id",
            "gene_id",
            "--fdr",
            "padj",
            "--pvalue",
            "pvalue",
            "--logfc",
            "log2FoldChange",
            "--stat",
            "stat",
        ]

        return_code, _, _ = Cmd["parse_diffexp.py"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error while parsing DGE results.")

        (Cmd["gzip"][deseq_output])()
        (Cmd["gzip"]["count_matrix.tab"])()

        outputs.raw = f"{deseq_output}.gz"
        outputs.de_json = "de_data.json"
        outputs.de_file = "de_file.tab.gz"
        outputs.count_matrix = "count_matrix.tab.gz"
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
                "DESeq2",
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
