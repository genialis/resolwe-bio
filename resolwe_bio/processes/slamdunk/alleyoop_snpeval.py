"""Run Alleyoop snpeval tool on Slamdunk results."""
import os

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    IntegerField,
    Process,
    StringField,
)


class AlleyoopSnpEval(Process):
    """Run Alleyoop snpeval."""

    slug = "alleyoop-snpeval"
    process_type = "data:alleyoop:snpeval"
    name = "Alleyoop snpeval"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/slamdunk:2.0.0"},
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "Slamdunk"
    data_name = '{{ slamdunk|sample_name|default("?") }}'
    version = "1.2.1"

    class Input:
        """Input fields for AlleyoopSnpEval."""

        ref_seq = DataField(
            "seq:nucleotide", label="FASTA file containig sequences for aligning"
        )
        regions = DataField(
            "bed", label="BED file with coordinates of regions of interest"
        )
        slamdunk = DataField("alignment:bam:slamdunk", label="Slamdunk results")
        read_length = IntegerField(
            label="Maximum read length",
            description="Maximum length of reads in the input FASTQ file",
            default=150,
        )

    class Output:
        """Output fields to process AlleyoopSnpEval."""

        report = FileField(
            label="Tab-separated file with read counts, T>C read counts and SNP indication"
        )
        plot = FileField(label="SNP evaluation plot")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.slamdunk.output.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]

        args = [
            "-o",
            "snpeval",
            "-r",
            inputs.ref_seq.output.fasta.path,
            "-b",
            inputs.regions.output.bed.path,
            "-s",
            ".",
            "-l",
            inputs.read_length,
        ]

        (Cmd["ln"]["-s", inputs.slamdunk.output.variants.path, f"{name}_snp.vcf"])()

        return_code, _, _ = Cmd["alleyoop"]["snpeval"][args][
            inputs.slamdunk.output.bam.path
        ] & TEE(retcode=None)
        if return_code:
            self.error("Alleyoop snpeval analysis failed.")

        snp_file = os.path.join("snpeval", f"{name}_SNPeval.csv")
        snp_file_renamed = os.path.join("snpeval", f"{name}_SNPeval.txt")
        os.rename(snp_file, snp_file_renamed)

        outputs.report = snp_file_renamed
        outputs.plot = os.path.join("snpeval", f"{name}_SNPeval.pdf")
        outputs.species = inputs.slamdunk.output.species
        outputs.build = inputs.slamdunk.output.build
