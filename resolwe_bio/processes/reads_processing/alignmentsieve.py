"""Filter BAM files using alignmentSieve from deeptools."""
import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class AlignmentSieve(Process):
    """Filter alignments of BAM files according to specified parameters.

    Program is bundled with deeptools. See [documentation](
    https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html)
    for more details.
    """

    slug = "alignmentsieve"
    name = "alignmentSieve"
    process_type = "data:alignment:bam:sieve"
    version = "1.2.1"
    category = "Alignment Filtering"
    data_name = 'Sieved BAM ({{ alignment|sample_name|default("?") }})'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.6.0"}
        },
        "resources": {
            "cores": 10,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process AlignmentSieve."""

        alignment = DataField(data_type="alignment:bam", label="Alignment BAM file")
        min_fragment_length = IntegerField(
            label="--minFragmentLength",
            description="The minimum fragment length needed for "
            "read/pair inclusion. This option is primarily useful in "
            "ATACseq experiments, for filtering mono- or di-nucleosome "
            "fragments. (Default: 0)",
            default=0,
        )
        max_fragment_length = IntegerField(
            label="--maxFragmentLength",
            description="The maximum fragment length needed for "
            "read/pair inclusion. A value of 0 indicates "
            "no limit. (Default: 0)",
            default=0,
        )

        class BigWigOptions:
            """Options for calculating BigWig."""

            bigwig_binsize = IntegerField(
                label="BigWig bin size",
                description="Size of the bins, in bases, for the output of the "
                "bigwig/bedgraph file. Default is 50.",
                default=50,
            )
            bigwig_timeout = IntegerField(
                label="BigWig timeout",
                description="Number of seconds before creation of BigWig timeouts. "
                "Default is after 480 seconds (8 minutes).",
                default=480,
            )

        bigwig_opts = GroupField(BigWigOptions, label="BigWig options")

    class Output:
        """Output fields to process AlignmentSieve."""

        bam = FileField(label="Sieved BAM file")
        bai = FileField(label="Index of sieved BAM file")
        stats = FileField(label="Alignment statistics")
        bigwig = FileField(label="BigWig file", required=False)
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run code of AlignmentSieve process."""

        basename = Path(inputs.alignment.output.bam.path).name
        assert basename.endswith(".bam")
        name = basename[:-4]

        bam_species = inputs.alignment.output.species
        output_bam_file = f"{name}_filtered.bam"
        filter_metrics = f"{name}_metrics.txt"
        bam_bigwig = f"{name}_filtered.bw"

        params = [
            "--verbose",
            "--bam",
            inputs.alignment.output.bam.path,
            "--outFile",
            output_bam_file,
            "--filterMetrics",
            filter_metrics,
            "--maxFragmentLength",
            inputs.max_fragment_length,
            "--minFragmentLength",
            inputs.min_fragment_length,
        ]

        return_code, _, _ = Cmd["alignmentSieve"][params] & TEE(retcode=None)

        if return_code:
            self.error("Error sieving the BAM file.")

        if not os.path.exists(output_bam_file):
            self.error(f"File {output_bam_file} not created.")

        outputs.bam = output_bam_file
        self.progress(0.3)

        return_code, stdout, stderr = Cmd["samtools"]["index"][output_bam_file] & TEE(
            retcode=None
        )

        if return_code:
            print(stderr)
            self.error(f"Failed to index {output_bam_file}.")

        outputs.bai = f"{output_bam_file}.bai"
        self.progress(0.5)

        stats = f"{name}_stats.txt"
        (Cmd["samtools"]["flagstat"][f"{output_bam_file}"] > stats)()

        outputs.stats = stats
        self.progress(0.75)

        btb_inputs = [
            output_bam_file,
            bam_species,
            self.requirements.resources.cores,
            "deeptools",
            inputs.bigwig_opts.bigwig_binsize,
            inputs.bigwig_opts.bigwig_timeout,
        ]

        check_bigwig, _, _ = Cmd["bamtobigwig.sh"][btb_inputs] & TEE(retcode=None)

        self.progress(0.99)

        if check_bigwig:
            self.info("BigWig file not calculated.")
        else:
            outputs.bigwig = bam_bigwig

        outputs.species = bam_species
        outputs.build = inputs.alignment.output.build
