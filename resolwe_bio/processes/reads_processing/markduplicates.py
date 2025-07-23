"""Remove duplicates using GATK's MarkDuplicates."""

import os
import shutil

from resolwe.process import (
    BooleanField,
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


class MarkDuplicates(Process):
    """Remove duplicate reads from BAM file.

    Tool from Picard, wrapped by GATK4. See GATK MarkDuplicates for more information.
    """

    slug = "markduplicates"
    name = "MarkDuplicates"
    process_type = "data:alignment:bam:markduplicate:"
    version = "1.7.1"
    category = "BAM processing"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
    }
    data_name = "{{ bam|name|default('?') }}"
    persistence = Persistence.CACHED

    class Input:
        """Input fields to process MarkDuplicates."""

        bam = DataField("alignment:bam", label="Alignment BAM file")
        skip = BooleanField(
            label="Skip MarkDuplicates step",
            description="MarkDuplicates step can be skipped.",
            default=False,
        )
        remove_duplicates = BooleanField(
            label="Remove duplicates",
            description="If true do not write duplicates to the output file "
            "instead of writing them with appropriate flags set.",
            default=False,
        )
        validation_stringency = StringField(
            label="Validation stringency",
            description="Validation stringency for all SAM files read by this "
            "program. Setting stringency to SILENT can improve "
            "performance when processing a BAM file in which "
            "variable-length data (read, qualities, tags) do not "
            "otherwise need to be decoded. Default is STRICT.",
            choices=[
                ("STRICT", "STRICT"),
                ("LENIENT", "LENIENT"),
                ("SILENT", "SILENT"),
            ],
            default="STRICT",
        )
        assume_sort_order = StringField(
            label="Assume sort order",
            description="If not null (default), assume that the input file "
            "has this order even if the header says otherwise."
            "Possible values are unsorted, queryname, coordinate "
            "and unknown.",
            choices=[
                ("", "as in BAM header (default)"),
                ("unsorted", "unsorted"),
                ("queryname", "queryname"),
                ("coordinate", "coordinate"),
                ("duplicate", "duplicate"),
                ("unknown", "unknown"),
            ],
            default="",
        )

        class Advanced:
            """Advanced options."""

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

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to process MarkDuplicates."""

        bam = FileField(label="Marked duplicates BAM file")
        bai = FileField(label="Index of marked duplicates BAM file")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")
        metrics_file = FileField(label="Metrics from MarkDuplicate process")

    def run(self, inputs, outputs):
        """Run analysis.

        Note that this process can have output as two different filenames. If the process
        is skipped, there will be no modification of the filename, whereas if
        markduplication goes through, it will append 'markduplicates'.
        """
        TMPDIR = os.environ.get("TMPDIR")

        # Prepare output file names.
        file_name = os.path.splitext(os.path.basename(inputs.bam.output.bam.path))[0]
        metrics_file = f"{file_name}_metrics.txt"
        # We do not append anything particular to this object (like we did with e.g.
        # _metrics.txt" for metrics file), because this step can be skipped and the
        # name (e.g. "file_markduplicated.bam") would only confuse the matter.

        species = inputs.bam.output.species
        build = inputs.bam.output.build

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        if not inputs.skip:
            if inputs.remove_duplicates:
                rmd = "true"
            else:
                rmd = "false"

            bam = file_name + ".markduplicates.bam"
            md_inputs = [
                "--java-options",
                f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                "--INPUT",
                f"{inputs.bam.output.bam.path}",
                "--VALIDATION_STRINGENCY",
                f"{inputs.validation_stringency}",
                "--OUTPUT",
                f"{bam}",
                "--METRICS_FILE",
                f"{metrics_file}",
                "--TMP_DIR",
                TMPDIR,
            ]

            if inputs.remove_duplicates:
                md_inputs.extend(["--REMOVE_DUPLICATES", f"{rmd}"])

            if inputs.assume_sort_order:
                md_inputs.extend(["--ASSUME_SORT_ORDER", f"{inputs.assume_sort_order}"])

            Cmd["gatk"]["MarkDuplicates"](md_inputs)
        else:
            # Process skipped, output filename matches input.
            bam = os.path.basename(inputs.bam.output.bam.path)
            shutil.copy2(inputs.bam.output.bam.path, bam)

            with open(metrics_file, "w") as f:
                f.write("MarkDuplicate process skipped.")

            if os.path.exists(metrics_file):
                print(f"{metrics_file} created.")

        self.progress(0.5)

        stderr_file = "stderr.txt"
        (Cmd["samtools"]["index"][f"{bam}"] > stderr_file)()

        if not os.path.exists(f"{bam}.bai"):
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

        self.progress(0.9)

        outputs.bam = bam
        outputs.bai = bam + ".bai"
        outputs.stats = stats

        outputs.species = species
        outputs.build = build
        outputs.metrics_file = metrics_file
