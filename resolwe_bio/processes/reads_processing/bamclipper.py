"""Remove primer sequences from BAM alignments by soft-clipping."""
import os
from shutil import copy2

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    Process,
    SchedulingClass,
    StringField,
)


class Bamclipper(Process):
    """Remove primer sequence from BAM alignments by soft-clipping.

    This process is a wrapper for bamclipper which can be found at
    https://github.com/tommyau/bamclipper.
    """

    slug = "bamclipper"
    name = "Bamclipper"
    process_type = "data:alignment:bam:bamclipped:"
    version = "1.2.1"
    category = "Clipping"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:5.2.0"}
        },
    }
    data_name = '{{ alignment|sample_name|default("?") }}'

    class Input:
        """Input fields to process Bamclipper."""

        alignment = DataField("alignment:bam", label="Alignment BAM file")
        bedpe = DataField("bedpe", label="BEDPE file", required=False)
        skip = BooleanField(
            label="Skip Bamclipper step",
            description="Use this option to skip Bamclipper step.",
            default=False,
        )

    class Output:
        """Output fields to process Bamclipper."""

        bam = FileField(label="Clipped BAM file")
        bai = FileField(label="Index of clipped BAM file")
        stats = FileField(label="Alignment statistics")
        bigwig = FileField(label="BigWig file", required=False)
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        bam_build = inputs.alignment.output.build
        bam_species = inputs.alignment.output.species

        name = os.path.splitext(os.path.basename(inputs.alignment.output.bam.path))[0]

        # If so specified, skip bamclipper step. Prepare outputs to match those of as if bamclipping proceeded.
        if inputs.skip:
            bam = f"{name}.bamclipper_skipped.bam"
            bai = f"{bam}.bai"
            bigwig = f"{bam[:-4]}.bw"

            copy2(inputs.alignment.output.bam.path, bam)

            self.info("Skipping bamclipper step.")
        else:
            bedpe_build = inputs.bedpe.output.build
            bedpe_species = inputs.bedpe.output.species

            if bam_build != bedpe_build:
                self.error(
                    f"Builds of the genome {bam_build} and annotation "
                    f"{bedpe_build} do not match. Please provide genome and "
                    f"annotation with the same build."
                )

            if bam_species != bedpe_species:
                self.error(
                    f"Species of BAM ({bam_species}) and BEDPE ({bedpe_species}) "
                    "files do not match."
                )

            # Output of bamclipper.sh is a file that is appended a primerclipped just before file
            # extension, e.g. "file.bam" now becomes "file.primerclipped.bam".
            bc_inputs = [
                "-b",
                f"{inputs.alignment.output.bam.path}",
                "-p",
                f"{inputs.bedpe.output.bedpe.path}",
            ]
            Cmd["bamclipper.sh"](bc_inputs)

            bam = f"{name}.primerclipped.bam"
            bai = f"{bam}.bai"
            bigwig = f"{bam[:-4]}.bw"

            self.progress(0.5)

        # Calculate BAM statistics.
        stderr_file = "stderr.txt"
        (Cmd["samtools"]["index"][f"{bam}"] > stderr_file)()

        if not os.path.exists(f"{bai}"):
            self.error(f"Indexing of {bam} failed.")

        self.progress(0.7)

        # Print to console if errors have been generated.
        if os.path.exists(stderr_file):
            with open(stderr_file, "r") as f:
                all_lines = f.readlines()
                if len(all_lines) > 0:
                    for l in all_lines:
                        print(l)

        stats = f"{bam}_stats.txt"
        (Cmd["samtools"]["flagstat"][f"{bam}"] > stats)()

        self.progress(0.8)

        btb_inputs = [
            f"{bam}",
            f"{bam_species}",
            f"{self.requirements.resources.cores}",
        ]

        Cmd["bamtobigwig.sh"](btb_inputs)

        if not os.path.exists(bigwig):
            self.info("BigWig file not calculated.")
        else:
            outputs.bigwig = bigwig

        self.progress(0.9)

        outputs.bam = bam
        outputs.bai = bai
        outputs.stats = stats

        outputs.species = bam_species
        outputs.build = bam_build
