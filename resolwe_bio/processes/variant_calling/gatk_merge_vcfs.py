"""Run GATK MergeVcfs."""

import os

from plumbum import TEE

from resolwe.process import (
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


class GatkMergeVcfs(Process):
    """Combine multiple variant files into a single variant file using GATK MergeVcfs."""

    slug = "gatk-merge-vcfs"
    name = "GATK MergeVcfs"
    category = "GATK"
    process_type = "data:variants:vcf:mergevcfs"
    version = "1.2.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
            "storage": 200,
        },
    }
    data_name = "Combined variants"

    class Input:
        """Input fields for GatkMergeVcfs."""

        vcfs = ListField(DataField("variants:vcf"), label="Input data (VCFs)")

        class AdvancedOptions:
            """Advanced options."""

            ref_seq = DataField(
                "seq:nucleotide",
                label="Reference sequence",
                required=False,
                description="Optionally use a sequence dictionary file (.dict) "
                "if the input VCF does not contain a complete contig list.",
            )

            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of the garbage collectors.",
            )

            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GatkMergeVcfs."""

        vcf = FileField(label="Merged VCF")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        combined_variants = "combined_variants.vcf.gz"
        combined_variants_index = combined_variants + ".tbi"
        vcf_list = "input_variant_files.list"
        species = inputs.vcfs[0].output.species
        build = inputs.vcfs[0].output.build

        if not all(vcf.output.species == species for vcf in inputs.vcfs):
            self.error(
                "Species information must be the same for all of the input VCF files."
            )

        if not all(vcf.output.build == build for vcf in inputs.vcfs):
            self.error(
                "Genome build information must be the same for all of the input VCF files."
            )

        with open(vcf_list, "w") as input_vcfs:
            for vcf in inputs.vcfs:
                input_vcfs.write(f"{vcf.output.vcf.path}\n")

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced_options.java_gc_threads
        )
        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced_options.max_heap_size}g",
            f"I={vcf_list}",
            f"O={combined_variants}",
            f"TMP_DIR={TMPDIR}",
        ]

        if inputs.advanced_options.ref_seq:
            if inputs.advanced_options.ref_seq.output.species != species:
                self.error(
                    "The species information of the provided reference "
                    "sequence file does not match the species of the input VCFs."
                )
            if inputs.advanced_options.ref_seq.output.build != build:
                self.error(
                    "The genome build information of the provided reference "
                    "sequence file does not match the build of the input VCFs."
                )
            args.append(f"D={inputs.advanced_options.ref_seq.output.fasta_dict.path}")

        return_code, stdout, stderr = Cmd["gatk"]["MergeVcfs"][args] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("GATK MergeVcfs failed.")

        outputs.vcf = combined_variants
        outputs.tbi = combined_variants_index
        outputs.species = species
        outputs.build = build
