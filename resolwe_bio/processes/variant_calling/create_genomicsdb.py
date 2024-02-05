"""Run GATK GenomicsDBImport tool."""

import os
import shutil
from pathlib import Path

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
from resolwe.process.fields import DirField


class GenomicsDBImport(Process):
    """Import single-sample GVCFs into GenomicsDB before joint genotyping."""

    slug = "gatk-genomicsdb-import"
    name = "GATK GenomicsDBImport"
    category = "GATK"
    process_type = "data:genomicsdb"
    version = "1.3.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
            "storage": 200,
        },
    }
    data_name = '{{ "GATK GenomicsDB (%s %s)"|format(gvcfs|length, "samples added" if use_existing else "samples" ) }}'

    class Input:
        """Input fields for GenomicsDBImport."""

        gvcfs = ListField(
            DataField("variants:gvcf"),
            label="Input data (GVCF)",
        )

        intervals = DataField(
            "bed",
            label="Intervals file (.bed)",
            description="Intervals file is required if a new database will be "
            "created.",
            required=False,
        )

        use_existing = BooleanField(
            label="Add new samples to an existing GenomicsDB workspace",
            default=False,
        )

        existing_db = DataField(
            "genomicsdb",
            label="Select a GATK GenomicsDB object",
            description="Instead of creating a new database the GVCFs are "
            "added to this database and a new GenomicsDB object is created.",
            required=False,
            hidden="!use_existing",
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

            max_heap_size = IntegerField(
                label="Java maximum heap size in GB (Xmx)",
                default=28,
                description="Set the maximum Java heap size.",
            )

            use_cms_gc = BooleanField(
                label="Use CMS Garbage Collector in Java",
                default=True,
                description="The Concurrent Mark Sweep (CMS) implementation uses multiple garbage "
                "collector threads for garbage collection.",
            )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GenomicsDBImport."""

        database = DirField(label="GenomicsDB workspace")
        intervals = FileField(label="Intervals file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        database_folder = "database"
        sample_map_file = "sample_map.txt"

        species = inputs.gvcfs[0].output.species
        if any(gvcf.output.species != species for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples are of the same species.")

        build = inputs.gvcfs[0].output.build
        if any(gvcf.output.build != build for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples have the same genome build.")

        with open(sample_map_file, "w") as sample_map:
            for gvcf in inputs.gvcfs:
                sample_map.write(f"{gvcf.entity_name}\t{gvcf.output.vcf.path}\n")

        if inputs.use_existing and inputs.existing_db is None:
            self.error(
                "GATK GenomicsDB object has to be provided to add GVCFs to the existing "
                "database."
            )

        elif inputs.use_existing and inputs.existing_db:
            if species != inputs.existing_db.output.species:
                self.error("The existing database and GVCFs species differ.")
            if build != inputs.existing_db.output.build:
                self.error("The existing database and GVCFs build differ.")

            shutil.copytree(inputs.existing_db.output.database.path, database_folder)

            db_import_args = [
                "--genomicsdb-update-workspace-path",
                database_folder,
            ]
            intervals = Path(inputs.existing_db.output.intervals.path)

        elif inputs.intervals:
            db_import_args = [
                "--genomicsdb-workspace-path",
                database_folder,
                "-L",
                inputs.intervals.output.bed.path,
            ]
            intervals = Path(inputs.intervals.output.bed.path)

        else:
            self.error("Intervals file is required for creating a new database.")

        java_memory = min(
            int(self.requirements.resources.memory / 1024),
            inputs.advanced_options.max_heap_size,
        )
        java_options = f"-Xmx{java_memory}g"
        if inputs.advanced_options.use_cms_gc:
            java_options += " -XX:+UseConcMarkSweepGC"

        db_import_args.extend(
            [
                "--sample-name-map",
                sample_map_file,
                "--batch-size",
                inputs.advanced_options.batch_size,
                "--reader-threads",
                min(self.requirements.resources.cores, 5),
                "--verbosity",
                "DEBUG",
                "--tmp-dir",
                os.environ.get("TMPDIR"),
                "--java-options",
                java_options,
            ]
        )

        if inputs.advanced_options.consolidate:
            db_import_args.append("--consolidate")

        return_code, stdout, stderr = Cmd["gatk"]["GenomicsDBImport"][
            db_import_args
        ] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("GATK GenomicsDBImport tool failed.")

        output_bed = f"./{intervals.name}"
        Path(output_bed).symlink_to(str(intervals))

        outputs.intervals = output_bed
        outputs.database = database_folder
        outputs.species = species
        outputs.build = build
