"""Run GATK GenotypeGVCFs tool."""
from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
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
    version = "1.0.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields for GatkGenotypeGVCFs."""

        gvcfs = ListField(
            DataField("variants:gvcf"),
            label="Input data (GVCF)",
        )
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

        intervals = DataField(
            "bed",
            label="Intervals file (.bed).",
        )

    class Output:
        """Output fields for GatkGenotypeGVCFs."""

        vcf = FileField(label="GVCF file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        variants = "variants.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"

        species = inputs.gvcfs[0].output.species
        if not all(gvcf.output.species == species for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples are of the same species.")

        build = inputs.gvcfs[0].output.build
        if not all(gvcf.output.build == build for gvcf in inputs.gvcfs):
            self.error("Not all of the input samples have the same genome build.")

        with open("sample_map.txt", "w") as sample_map:
            for gvcf in inputs.gvcfs:
                sample_map.write(f"{gvcf.entity_name}\t{gvcf.output.vcf.path}\n")

        db_import_args = [
            "--genomicsdb-workspace-path",
            "database",
            "-L",
            inputs.intervals.output.bed.path,
            "--sample-name-map",
            "sample_map.txt",
        ]

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
            "-G",
            "StandardAnnotation",
            "-G",
            "AS_StandardAnnotation",
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
