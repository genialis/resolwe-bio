"""Run GATK GenotypeGVCFs tool."""

import os
from pathlib import Path

from joblib import Parallel, delayed, wrap_non_picklable_objects
from plumbum import TEE

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
from resolwe.process.fields import DirField

TMPDIR = os.environ.get("TMPDIR")


def create_vcf_path(interval_path):
    """Create a vcf path based on the interval name."""
    return f"cohort_variants/{interval_path.stem}.vcf"


@delayed
@wrap_non_picklable_objects
def run_genotype_gvcfs(interval_path, ref_seq_path, db_path, dbsnp_path, java_memory):
    """Run genotyping on a specifed interval."""
    variants_interval = create_vcf_path(interval_path)

    genotype_gvcfs_inputs = [
        "-R",
        ref_seq_path,
        "-V",
        f"gendb://{db_path}",
        "-O",
        variants_interval,
        "-L",
        interval_path,
        "-D",
        dbsnp_path,
        "-G",
        "StandardAnnotation",
        "-G",
        "AS_StandardAnnotation",
        "--only-output-calls-starting-in-intervals",
        "--java-options",
        f"-Xmx{java_memory}g",
        "--tmp-dir",
        TMPDIR,
    ]

    return_code, _, _ = Cmd["gatk"]["GenotypeGVCFs"][genotype_gvcfs_inputs] & TEE(
        retcode=None
    )
    return return_code


class GatkGenotypeGVCFs(Process):
    """Consolidate GVCFs and run joint calling using GenotypeGVCFs tool."""

    slug = "gatk-genotype-gvcfs"
    name = "GATK GenotypeGVCFs"
    category = "GATK"
    process_type = "data:variants:vcf:genotypegvcfs"
    version = "2.3.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 16,
            "memory": 32768,
            "storage": 200,
        },
    }
    data_name = "Cohort variants"

    class Input:
        """Input fields for GatkGenotypeGVCFs."""

        database = DataField("genomicsdb", label="GATK GenomicsDB")
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

        dbsnp = DataField("variants:vcf", label="dbSNP file")

        class AdvancedOptions:
            """Advanced options."""

            n_jobs = IntegerField(
                label="Number of concurent jobs",
                description="Use a fixed number of jobs for genotyping "
                "instead of determining it based on the number of available "
                "cores.",
                required=False,
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size in GB (Xmx)",
                default=28,
                description="Set the maximum Java heap size.",
            )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GatkGenotypeGVCFs."""

        vcf = FileField(label="GVCF file")
        vcf_dir = DirField(label="Folder with split GVCFs")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        variants = "cohort_variants.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"
        intervals_path = Path("intervals_folder")

        Path("cohort_variants").mkdir(exist_ok=True)
        intervals_path.mkdir(exist_ok=True)

        if inputs.advanced_options.n_jobs:
            n_jobs = max(inputs.advanced_options.n_jobs, 1)
        else:
            n_jobs = max(self.requirements.resources.cores, 1)

        split_intervals_inputs = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-L",
            inputs.database.output.intervals.path,
            "--scatter-count",
            n_jobs,
            "-O",
            str(intervals_path),
            "--tmp-dir",
            TMPDIR,
        ]
        return_code, _, _ = Cmd["gatk"]["SplitIntervals"][split_intervals_inputs] & TEE(
            retcode=None
        )

        if return_code is None:
            self.error("Could not create equally sized intervals from the bedfile.")

        intervals = [path for path in intervals_path.glob("*.interval_list")]

        java_memory = min(
            int(self.requirements.resources.memory / 1024),
            inputs.advanced_options.max_heap_size,
        )

        return_codes = Parallel(n_jobs=n_jobs)(
            run_genotype_gvcfs(
                interval_path,
                ref_seq_path=inputs.ref_seq.output.fasta.path,
                db_path=inputs.database.output.database.path,
                dbsnp_path=inputs.dbsnp.output.vcf.path,
                java_memory=java_memory,
            )
            for interval_path in intervals
        )

        outputs.vcf_dir = "cohort_variants"
        if any(return_codes):
            self.error("GATK GenotypeGVCFs tool failed.")

        merge_list_file = "merge_files.list"
        with open(merge_list_file, "w") as f:
            for interval_path in sorted(intervals):
                f.write(f"{create_vcf_path(interval_path)}\n")

        self.progress(0.8)

        merge_inputs = [
            "-I",
            merge_list_file,
            "-O",
            variants,
            "-R",
            inputs.ref_seq.output.fasta.path,
            "--CREATE_INDEX",
            "false",
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["GatherVcfs"][merge_inputs] & TEE(retcode=None)
        if return_code:
            self.error("GATK GatherVcfs tool failed.")

        # Compress and index the output variants file
        (Cmd["bgzip"]["-c", variants] > variants_gz)()
        Cmd["tabix"]["-p", "vcf", variants_gz]()

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = inputs.database.output.species
        outputs.build = inputs.database.output.build
