"""Run GATK GenotypeGVCFs tool."""
from plumbum import TEE

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


class GatkGenotypeGVCFs(Process):
    """Consolidate GVCFs and run joint calling using GenotypeGVCFs tool."""

    slug = "gatk-genotype-gvcfs"
    name = "GATK GenotypeGVCFs"
    category = "GATK"
    process_type = "data:variants:vcf:genotypegvcfs"
    version = "1.0.1"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
        },
    }
    data_name = "Cohort variants"

    class Input:
        """Input fields for GatkGenotypeGVCFs."""

        gvcfs = ListField(
            DataField("variants:gvcf"),
            label="Input data (GVCF)",
        )
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

        intervals = DataField(
            "bed",
            label="Intervals file (.bed)",
        )

        dbsnp = DataField("variants:vcf", label="dbSNP file")

        advanced = BooleanField(
            label="Show advanced options",
            description="Inspect and modify parameters.",
            default=False,
        )

        class AdvancedOptions:
            """Advanced options."""

            batch_size = IntegerField(
                label="Batch size",
                default=0,
                description="Batch size controls the number of samples "
                "for which readers are open at once and therefore provides "
                "a way to minimize memory consumption. However, it can "
                "take longer to complete. Use the consolidate flag if more "
                "than a hundred batches were used. This will improve feature "
                "read time. batchSize=0 means no batching "
                "(i.e. readers for all samples will be opened at once).",
            )

            consolidate = BooleanField(
                label="Consolidate",
                default=False,
                description="Boolean flag to enable consolidation. If "
                "importing data in batches, a new fragment is created for "
                "each batch. In case thousands of fragments are created, "
                "GenomicsDB feature readers will try to open ~20x as many "
                "files. Also, internally GenomicsDB would consume more "
                "memory to maintain bookkeeping data from all fragments. "
                "Use this flag to merge all fragments into one. Merging "
                "can potentially improve read performance, however overall "
                "benefit might not be noticeable as the top Java layers "
                "have significantly higher overheads. This flag has no "
                "effect if only one batch is used.",
            )

        advanced_options = GroupField(
            AdvancedOptions, label="Advanced options", hidden="!advanced"
        )

    class Output:
        """Output fields for GatkGenotypeGVCFs."""

        vcf = FileField(label="GVCF file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        variants = "cohort_variants.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"

        sample_map_file = "sample_map.txt"

        species = inputs.gvcfs[0].output.species
        if not all(gvcf.output.species == species for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples are of the same species.")

        build = inputs.gvcfs[0].output.build
        if not all(gvcf.output.build == build for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples have the same genome build.")

        with open(sample_map_file, "w") as sample_map:
            for gvcf in inputs.gvcfs:
                sample_map.write(f"{gvcf.entity_name}\t{gvcf.output.vcf.path}\n")

        db_import_args = [
            "--genomicsdb-workspace-path",
            "database",
            "-L",
            inputs.intervals.output.bed.path,
            "--sample-name-map",
            sample_map_file,
            "--batch-size",
            inputs.advanced_options.batch_size,
            "--reader-threads",
            min(self.requirements.resources.cores, 5),
        ]

        if inputs.advanced_options.consolidate:
            db_import_args.append("--seqBias")

        return_code, _, _ = Cmd["gatk"]["GenomicsDBImport"][db_import_args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK GenomicsDBImport tool failed.")

        genotype_gvcfs_inputs = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-V",
            "gendb://database",
            "-O",
            variants,
            "-L",
            inputs.intervals.output.bed.path,
            "-D",
            inputs.dbsnp.output.vcf.path,
            "-G",
            "StandardAnnotation",
            "-G",
            "AS_StandardAnnotation",
            "--only-output-calls-starting-in-intervals",
        ]

        return_code, _, _ = Cmd["gatk"]["GenotypeGVCFs"][genotype_gvcfs_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK GenotypeGVCFs tool failed.")

        # Compress and index the output variants file
        (Cmd["bgzip"]["-c", variants] > variants_gz)()
        Cmd["tabix"]["-p", "vcf", variants_gz]()

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = species
        outputs.build = build
