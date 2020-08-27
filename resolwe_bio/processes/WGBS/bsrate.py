"""Calculate conversion rate for bisulfite sequencing."""
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    IntegerField,
    Process,
    SchedulingClass,
)


class BsConversionRate(Process):
    """Estimate bisulfite conversion rate in a control set.

    The program bsrate included in [Methpipe]
    (https://github.com/smithlabcode/methpipe) will estimate the bisulfite
    conversion rate.
    """

    slug = "bs-conversion-rate"
    name = "Bisulfite conversion rate"
    process_type = "data:wgbs:bsrate"
    version = "1.0.1"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/wgbs:1.2.0"}},
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    data_name = '{{ mr|sample_name|default("?") }}'

    class Input:
        """Input fields for BsConversionRate."""

        mr = DataField(
            "alignment:bam:walt",
            label="Aligned reads from bisulfite sequencing",
            description="Bisulfite specifc alignment such as WALT is required as .mr file type is used. Duplicates"
            "should be removed to reduce any bias introduced by incomplete conversion on PCR duplicate"
            "reads.",
        )
        skip = BooleanField(
            label="Skip Bisulfite conversion rate step",
            description="Bisulfite conversion rate step can be skipped.",
            default=False,
        )
        sequence = DataField(
            "seq:nucleotide",
            label="Unmethylated control sequence",
            description="Separate unmethylated control sequence FASTA file is required to estimate bisulfite"
            "conversion rate.",
            required=False,
        )
        count_all = BooleanField(
            label="Count all cytosines including CpGs", default=True
        )
        read_length = IntegerField(label="Average read length", default=150)
        max_mismatch = FloatField(
            label="Maximum fraction of mismatches", required=False
        )
        a_rich = BooleanField(label="Reads are A-rich", default=False)

    class Output:
        """Output fields."""

        report = FileField(label="Bisulfite conversion rate report")

    def run(self, inputs, outputs):
        """Run the analysis."""
        basename = os.path.basename(inputs.mr.mr.path)
        assert basename.endswith(".mr.gz")
        name = basename[:-6]
        report_file = f"{name}_spikein_bsrate.txt"

        skip_process = inputs.skip

        try:
            inputs.mr.spikein_mr.path
        except AttributeError:
            self.warning(
                "Selected sample lacks the alignment file for unmethylated control reads."
            )
            skip_process = True
        try:
            inputs.sequence.fasta.path
        except AttributeError:
            self.warning("Unmethylated control sequence was not provided.")
            skip_process = True

        if not skip_process:
            (Cmd["pigz"]["-cd", inputs.mr.spikein_mr.path] > f"{name}.mr")()

            args = [
                "-chrom",
                inputs.sequence.fasta.path,
                "-output",
                report_file,
            ]
            if inputs.count_all:
                args.append("-all")
            if inputs.max_mismatch:
                args.extend(["-max", inputs.max_mismatch])
            if inputs.a_rich:
                args.append("-a-rich")

            return_code, _, _ = Cmd["bsrate"][args][f"{name}.mr"] & TEE(retcode=None)
            if return_code:
                self.error("Bsrate analysis failed.")
        else:
            with open(report_file, "w") as f:
                f.write("Bisulfite conversion rate process skipped.")

        outputs.report = report_file
