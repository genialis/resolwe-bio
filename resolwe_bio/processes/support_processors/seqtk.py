"""Subsample reads from FASTQ files with Seqtk."""

from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
)


class SeqtkSampleSingle(Process):
    """Subsample reads from FASTQ file (single-end).

    [Seqtk](https://github.com/lh3/seqtk) is a fast and lightweight tool for
    processing sequences in the FASTA or FASTQ format. The Seqtk "sample" command
    enables subsampling of the large FASTQ file(s).
    """

    slug = "seqtk-sample-single"
    name = "Subsample FASTQ (single-end)"
    process_type = "data:reads:fastq:single:seqtk"
    version = "1.5.2"
    category = "FASTQ processing"
    data_name = "{{ reads|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process SeqtkSampleSingle."""

        reads = DataField(
            "reads:fastq:single",
            label="Reads",
        )
        n_reads = IntegerField(
            label="Number of reads",
            default=1000000,
        )

        class Advanced:
            """Advanced options."""

            seed = IntegerField(
                label="Seed",
                default=11,
            )
            fraction = FloatField(
                label="Fraction",
                required=False,
                description="Use the fraction of reads [0 - 1.0] from the "
                "original input file instead of the absolute number of reads. "
                "If set, this will override the 'Number of reads' input parameter.",
            )
            two_pass = BooleanField(
                label="2-pass mode",
                default=False,
                description="Enable two-pass mode when down-sampling. "
                "Two-pass mode is twice as slow but with much reduced memory.",
            )

        advanced = GroupField(
            Advanced,
            label="Advanced options",
        )

    class Output:
        """Output fields to process SeqtkSampleSingle."""

        fastq = ListField(
            FileField(),
            label="Remaining reads",
        )
        fastqc_url = ListField(
            FileHtmlField(),
            label="Quality control with FastQC",
        )
        fastqc_archive = ListField(
            FileField(),
            label="Download FastQC archive",
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        if inputs.advanced.fraction and not 0 < inputs.advanced.fraction <= 1.0:
            self.error("Fraction of reads should be between 0 and 1.")

        basename = Path(inputs.reads.output.fastq[0].path).name
        assert basename.endswith(".fastq.gz")
        name = basename[:-9]

        input_reads = "input_reads.fastq.gz"
        final_reads = name + "_downsampled.fastq"

        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]]
            > input_reads
        )()

        args = [
            "-s",
            inputs.advanced.seed,
            input_reads,
        ]

        if inputs.advanced.two_pass:
            args.append("-2")

        if inputs.advanced.fraction:
            args.append(inputs.advanced.fraction)
        else:
            args.append(inputs.n_reads)

        (Cmd["seqtk"]["sample"][args] > final_reads)()

        Cmd["pigz"][final_reads]()

        args_fastqc = [
            f"{final_reads}.gz",
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]

        return_code, _, _ = Cmd["fastqc.sh"][args_fastqc] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing FASTQC report.")

        outputs.fastq = [f"{final_reads}.gz"]


class SeqtkSamplePaired(Process):
    """Subsample reads from FASTQ files (paired-end).

    [Seqtk](https://github.com/lh3/seqtk) is a fast and lightweight tool for
    processing sequences in the FASTA or FASTQ format. The Seqtk "sample" command
    enables subsampling of the large FASTQ file(s).
    """

    slug = "seqtk-sample-paired"
    name = "Subsample FASTQ (paired-end)"
    process_type = "data:reads:fastq:paired:seqtk"
    version = "1.5.2"
    category = "FASTQ processing"
    data_name = "{{ reads|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process SeqtkSamplePaired."""

        reads = DataField(
            "reads:fastq:paired",
            label="Reads",
        )
        n_reads = IntegerField(
            label="Number of reads",
            default=1000000,
        )

        class Advanced:
            """Advanced options."""

            seed = IntegerField(
                label="Seed",
                default=11,
            )
            fraction = FloatField(
                label="Fraction",
                required=False,
                description="Use the fraction of reads [0 - 1.0] from the "
                "orignal input file instead of the absolute number of reads. "
                "If set, this will override the 'Number of reads' input parameter.",
            )
            two_pass = BooleanField(
                label="2-pass mode",
                default=False,
                description="Enable two-pass mode when down-sampling. "
                "Two-pass mode is twice as slow but with much reduced memory.",
            )

        advanced = GroupField(
            Advanced,
            label="Advanced options",
        )

    class Output:
        """Output fields to process SeqtkSamplePaired."""

        fastq = ListField(
            FileField(),
            label="Remaining mate 1 reads",
        )
        fastq2 = ListField(
            FileField(),
            label="Remaining mate 2 reads",
        )
        fastqc_url = ListField(
            FileHtmlField(),
            label="Mate 1 quality control with FastQC",
        )
        fastqc_url2 = ListField(
            FileHtmlField(),
            label="Mate 2 quality control with FastQC",
        )
        fastqc_archive = ListField(
            FileField(),
            label="Download mate 1 FastQC archive",
        )
        fastqc_archive2 = ListField(
            FileField(),
            label="Download mate 2 FastQC archive",
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        if inputs.advanced.fraction and not 0 < inputs.advanced.fraction <= 1.0:
            self.error("Fraction of reads should be between 0 and 1.")

        basename1 = Path(inputs.reads.output.fastq[0].path).name
        basename2 = Path(inputs.reads.output.fastq2[0].path).name
        assert basename1.endswith(".fastq.gz")
        assert basename2.endswith(".fastq.gz")
        name_mate1 = basename1[:-9]
        name_mate2 = basename2[:-9]

        input_mate1 = "input_mate1.fastq.gz"
        input_mate2 = "input_mate2.fastq.gz"
        final_mate1 = name_mate1 + "_downsampled.fastq"
        final_mate2 = name_mate2 + "_downsampled.fastq"

        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]]
            > input_mate1
        )()
        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq2]]
            > input_mate2
        )()

        args1 = [
            "-s",
            inputs.advanced.seed,
            input_mate1,
        ]
        args2 = [
            "-s",
            inputs.advanced.seed,
            input_mate2,
        ]

        for arg in [args1, args2]:
            if inputs.advanced.two_pass:
                arg.append("-2")

            if inputs.advanced.fraction:
                arg.append(inputs.advanced.fraction)
            else:
                arg.append(inputs.n_reads)

        (Cmd["seqtk"]["sample"][args1] > final_mate1)()
        (Cmd["seqtk"]["sample"][args2] > final_mate2)()

        Cmd["pigz"][final_mate1]()
        Cmd["pigz"][final_mate2]()

        args_fastqc1 = [
            f"{final_mate1}.gz",
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]
        args_fastqc2 = [
            f"{final_mate2}.gz",
            "fastqc",
            "fastqc_archive2",
            "fastqc_url2",
        ]

        for arg in [args_fastqc1, args_fastqc2]:
            return_code, _, _ = Cmd["fastqc.sh"][arg] & TEE(retcode=None)
            if return_code:
                self.error("Error while preparing FASTQC report.")

        outputs.fastq = [f"{final_mate1}.gz"]
        outputs.fastq2 = [f"{final_mate2}.gz"]
