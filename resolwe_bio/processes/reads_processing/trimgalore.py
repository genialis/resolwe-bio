"""Run Trim Galore tool on paired-end sequencing data."""

import os

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
    Process,
    SchedulingClass,
    StringField,
)


class TrimGalorePaired(Process):
    """Process paired-end sequencing reads with Trim Galore.

    Trim Galore is a wrapper script that makes use of the publicly
    available adapter trimming tool Cutadapt and FastQC for quality
    control once the trimming process has completed.

    Low-quality ends are trimmed from reads in addition to adapter
    removal in a single pass. If no sequence was supplied, Trim Galore
    will attempt to auto-detect the adapter which has been used. For
    this it will analyse the first 1 million sequences of the first
    specified file and attempt to find the first 12 or 13bp of the
    following standard adapters: Illumina: AGATCGGAAGAGC, Small RNA:
    TGGAATTCTCGG, Nextera: CTGTCTCTTATA.

    If no adapter contamination can be detected within the first 1
    million sequences, or in case of a tie between several different
    adapters, Trim Galore defaults to illumina adapters.

    For additional information see official
    [user guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).
    """

    slug = "trimgalore-paired"
    name = "Trim Galore (paired-end)"
    process_type = "data:reads:fastq:paired:trimgalore"
    version = "1.3.2"
    category = "FASTQ processing"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"},
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    data_name = "{{ reads|name|default('?') }}"

    class Input:
        """Input fields of trimGalorePaired."""

        reads = DataField("reads:fastq:paired", label="Select paired-end reads")

        class QualityTrimming:
            """Quality trimming options."""

            quality = IntegerField(
                label="Quality cutoff",
                description="Trim low-quality ends from reads based on phred score.",
                default=20,
            )
            nextseq = IntegerField(
                label="NextSeq/NovaSeq trim cutoff",
                description="NextSeq/NovaSeq-specific quality "
                "trimming. Trims also dark cycles appearing as "
                "high-quality G bases. This will set a specific "
                "quality cutoff, but qualities of G bases are ignored. "
                "This can not be used with Quality cutoff and will "
                "override it.",
                required=False,
            )
            phred = StringField(
                label="Phred score encoding",
                description="Use either ASCII+33 quality scores as "
                "Phred scores (Sanger/Illumina 1.9+ encoding) or "
                "ASCII+64 quality scores (Illumina 1.5 encoding) for "
                "quality trimming",
                choices=[
                    ("--phred33", "ASCII+33"),
                    ("--phred64", "ASCII+64"),
                ],
                default="--phred33",
            )
            min_length = IntegerField(
                label="Minimum length after trimming",
                description="Discard reads that became shorter than "
                "selected length because of either quality or adapter "
                "trimming. Both reads of a read-pair need to be longer "
                "than specified length to be printed out to validated "
                "paired-end files. If only one read became too short "
                "there is the possibility of keeping such unpaired "
                "single-end reads with Retain unpaired. A value of 0 "
                "disables filtering based on length.",
                default=20,
            )
            max_n = IntegerField(
                label="Maximum number of Ns",
                description="Read exceeding this limit will result in "
                "the entire pair being removed from the trimmed output "
                "files.",
                required=False,
            )
            retain_unpaired = BooleanField(
                label="Retain unpaired reads after trimming",
                description="If only one of the two paired-end reads "
                "became too short, the longer read will be written.",
                default=False,
            )
            unpaired_len_1 = IntegerField(
                label="Unpaired read length cutoff for mate 1",
                default=35,
                hidden="!quality_trim.retain_unpaired",
            )
            unpaired_len_2 = IntegerField(
                label="Unpaired read length cutoff for mate 2",
                default=35,
                hidden="!quality_trim.retain_unpaired",
            )
            clip_r1 = IntegerField(
                label="Trim bases from 5' end of read 1",
                description="This may be useful if the qualities were "
                "very poor, or if there is some sort of unwanted bias "
                "at the 5' end.",
                required=False,
            )
            clip_r2 = IntegerField(
                label="Trim bases from 5' end of read 2",
                description="This may be useful if the qualities were "
                "very poor, or if there is some sort of unwanted bias "
                "at the 5' end. For paired-end bisulfite sequencing, "
                "it is recommended to remove the first few bp because "
                "the end-repair reaction may introduce a bias towards "
                "low methylation.",
                required=False,
            )
            three_prime_r1 = IntegerField(
                label="Trim bases from 3' end of read 1",
                description="Remove bases from the 3' end of read 1 "
                "after adapter/quality trimming has been performed. "
                "This may remove some unwanted bias from the 3' end "
                "that is not directly related to adapter sequence or "
                "basecall quality.",
                required=False,
            )
            three_prime_r2 = IntegerField(
                label="Trim bases from 3' end of read 2",
                description="Remove bases from the 3' end of read 2 "
                "after adapter/quality trimming has been performed. "
                "This may remove some unwanted bias from the 3' end "
                "that is not directly related to adapter sequence or "
                "basecall quality.",
                required=False,
            )

        class AdapterTrimming:
            """Adapter trimming options."""

            adapter = ListField(
                StringField(),
                label="Read 1 adapter sequence",
                description="Adapter sequences to be trimmed. "
                "Also see universal adapters field for predefined "
                "adapters. This is mutually exclusive with read 1 "
                "adapters file and universal adapters.",
                required=False,
                default=[],
            )
            adapter_2 = ListField(
                StringField(),
                label="Read 2 adapter sequence",
                description="Optional adapter sequence to be trimmed "
                "off read 2 of paired-end files. This is mutually "
                "exclusive with read 2 adapters file and universal "
                "adapters.",
                required=False,
                default=[],
            )
            adapter_file_1 = DataField(
                "seq:nucleotide",
                label="Read 1 adapters file",
                description="This is mutually exclusive with read 1 "
                "adapters and universal adapters.",
                required=False,
            )
            adapter_file_2 = DataField(
                "seq:nucleotide",
                label="Read 2 adapters file",
                description="This is mutually exclusive with read 2 "
                "adapters and universal adapters.",
                required=False,
            )
            universal_adapter = StringField(
                label="Universal adapters",
                description="Instead of default detection use specific "
                "adapters. Use 13bp of the Illumina universal adapter, "
                "12bp of the Nextera adapter or 12bp of the Illumina "
                "Small RNA 3' Adapter. Selecting to trim smallRNA "
                "adapters will also lower the length value to 18bp. "
                "If the smallRNA libraries are paired-end then read 2 "
                "adapter will be set to the Illumina small RNA 5' "
                "adapter automatically (GATCGTCGGACT) unless defined "
                "explicitly. This is mutually exclusive with manually "
                "defined adapters and adapter files.",
                choices=[
                    ("--illumina", "Illumina"),
                    ("--nextera", "Nextera"),
                    ("--small_rna", "Illumina small RNA"),
                ],
                required=False,
            )
            stringency = IntegerField(
                label="Overlap with adapter sequence required to trim",
                description="Defaults to a very stringent setting of "
                "1, i.e. even a single base pair of overlapping "
                "sequence will be trimmed of the 3' end of any read.",
                default=1,
            )
            error_rate = FloatField(
                label="Maximum allowed error rate",
                description="Number of errors divided by the length of "
                "the matching region",
                default=0.1,
            )

        class HardTrimming:
            """Hard trim options."""

            trim_5 = IntegerField(
                label="Hard trim sequences from 3' end",
                description="Instead of performing adapter-/quality "
                "trimming, this option will simply hard-trim sequences "
                "to bp from the 3' end. This is incompatible with "
                "other hard trimming options.",
                required=False,
            )
            trim_3 = IntegerField(
                label="Hard trim sequences from 5' end",
                description="Instead of performing adapter-/quality "
                "trimming, this option will simply hard-trim sequences "
                "to bp from the 5' end. This is incompatible with "
                "other hard trimming options.",
                required=False,
            )

        adapter_trim = GroupField(AdapterTrimming, label="Adapter trimming")
        quality_trim = GroupField(QualityTrimming, label="Quality trimming")
        hard_trim = GroupField(HardTrimming, label="Hard trimming")

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Remaining mate 1 reads")
        fastq2 = ListField(FileField(), label="Remaining mate 2 reads")
        report = FileField(label="Trim galore report", required=False)
        fastqc_url = ListField(
            FileHtmlField(), label="Mate 1 quality control with FastQC"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Mate 2 quality control with FastQC"
        )
        fastqc_archive = ListField(FileField(), label="Download mate 1 FastQC archive")
        fastqc_archive2 = ListField(FileField(), label="Download mate 2 FastQC archive")

    def run(self, inputs, outputs):
        """Run analysis."""

        if inputs.adapter_trim.adapter and inputs.adapter_trim.adapter_file_1:
            self.error(
                "Mate 1 adapters should be either a sequence or a file, but not both."
            )

        if inputs.adapter_trim.adapter_2 and inputs.adapter_trim.adapter_file_2:
            self.error(
                "Mate 2 adapters should be either a sequence or a file, but not both."
            )

        if inputs.hard_trim.trim_5 and inputs.hard_trim.trim_3:
            self.error("Only one type of hard trimming can be performed at once.")

        mate1_path = os.path.basename(inputs.reads.output.fastq[0].path)
        assert mate1_path.endswith(".fastq.gz")
        name_mate1 = mate1_path[:-9]
        mate2_path = os.path.basename(inputs.reads.output.fastq2[0].path)
        assert mate2_path.endswith(".fastq.gz")
        name_mate2 = mate2_path[:-9]

        merged_r1 = "input_reads_mate1.fastq.gz"
        merged_r1_name = merged_r1[:-9]
        merged_r2 = "input_reads_mate2.fastq.gz"
        merged_r2_name = merged_r2[:-9]
        (Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]] > merged_r1)()
        (Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq2]] > merged_r2)()

        params = [
            "--paired",
            "--cores",
            2,  # Actual core usage is 9, for more details see Trim Galore user guide
            inputs.quality_trim.phred,
            "--length",
            inputs.quality_trim.min_length,
            "--stringency",
            inputs.adapter_trim.stringency,
            "-e",
            inputs.adapter_trim.error_rate,
        ]

        if inputs.quality_trim.nextseq:
            params.extend(["--nextseq", inputs.quality_trim.nextseq])
        else:
            params.extend(["--quality", inputs.quality_trim.quality])

        if inputs.quality_trim.max_n:
            params.extend(["--max_n", inputs.quality_trim.max_n])
        if inputs.quality_trim.retain_unpaired:
            params.extend(
                [
                    "--retain_unpaired",
                    "--length_1",
                    inputs.quality_trim.unpaired_len_1,
                    "--length_2",
                    inputs.quality_trim.unpaired_len_2,
                ]
            )
        if inputs.quality_trim.clip_r1:
            params.extend(["--clip_R1", inputs.quality_trim.clip_r1])
        if inputs.quality_trim.clip_r2:
            params.extend(["--clip_R2", inputs.quality_trim.clip_r2])
        if inputs.quality_trim.three_prime_r1:
            params.extend(["--three_prime_clip_R1", inputs.quality_trim.three_prime_r1])
        if inputs.quality_trim.three_prime_r2:
            params.extend(["--three_prime_clip_R2", inputs.quality_trim.three_prime_r2])

        if inputs.adapter_trim.adapter:
            for adapter in inputs.adapter_trim.adapter:
                params.extend(["--adapter", adapter])
        if inputs.adapter_trim.adapter_2:
            for adapter in inputs.adapter_trim.adapter_2:
                params.extend(["--adapter2", adapter])
        if inputs.adapter_trim.adapter_file_1:
            params.extend(
                [
                    "--adapter",
                    f"file:{inputs.adapter_trim.adapter_file_1.output.fasta.path}",
                ]
            )
        if inputs.adapter_trim.adapter_file_2:
            params.extend(
                [
                    "--adapter2",
                    f"file:{inputs.adapter_trim.adapter_file_2.output.fasta.path}",
                ]
            )
        if inputs.adapter_trim.universal_adapter:
            params.append(inputs.adapter_trim.universal_adapter)
            if any(
                [
                    inputs.adapter_trim.adapter,
                    inputs.adapter_trim.adapter_2,
                    inputs.adapter_trim.adapter_file_1,
                    inputs.adapter_trim.adapter_file_2,
                ]
            ):
                self.error(
                    "You can not supply custom adapter sequence and use "
                    "the universal adapter sequence."
                )

        if inputs.hard_trim.trim_5 or inputs.hard_trim.trim_3:
            self.info(
                "Only hard trimming was performed. Skipped quality and adapter trimming."
            )
            if inputs.hard_trim.trim_5:
                params.extend(["--hardtrim5", inputs.hard_trim.trim_5])
            if inputs.hard_trim.trim_3:
                params.extend(["--hardtrim3", inputs.hard_trim.trim_3])

        return_code, _, _ = Cmd["trim_galore"][params][merged_r1, merged_r2] & TEE(
            retcode=None
        )
        if return_code:
            self.error("Error while trimming reads.")

        self.progress(0.7)

        trimmed_r1 = f"{name_mate1}_trim.fastq.gz"
        trimmed_r2 = f"{name_mate2}_trim.fastq.gz"
        if inputs.hard_trim.trim_5:
            os.rename(
                f"{merged_r1_name}.{inputs.hard_trim.trim_5}bp_5prime.fq.gz", trimmed_r1
            )
            os.rename(
                f"{merged_r2_name}.{inputs.hard_trim.trim_5}bp_5prime.fq.gz", trimmed_r2
            )

        elif inputs.hard_trim.trim_3:
            os.rename(
                f"{merged_r1_name}.{inputs.hard_trim.trim_3}bp_3prime.fq.gz", trimmed_r1
            )
            os.rename(
                f"{merged_r2_name}.{inputs.hard_trim.trim_3}bp_3prime.fq.gz", trimmed_r2
            )
        else:
            os.rename(f"{merged_r1_name}_val_1.fq.gz", trimmed_r1)
            os.rename(f"{merged_r2_name}_val_2.fq.gz", trimmed_r2)

            trim_report = "trim_galore_report.txt"
            (
                Cmd["cat"][
                    f"{merged_r1}_trimming_report.txt",
                    f"{merged_r2}_trimming_report.txt",
                ]
                > trim_report
            )()
            outputs.report = trim_report

        fastqc_args = [
            trimmed_r1,
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]

        self.progress(0.8)

        return_code, _, _ = Cmd["fastqc.sh"][fastqc_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing FASTQC report.")

        fastqc_args = [
            trimmed_r2,
            "fastqc",
            "fastqc_archive2",
            "fastqc_url2",
        ]
        return_code, _, _ = Cmd["fastqc.sh"][fastqc_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing FASTQC report.")

        self.progress(0.9)

        outputs.fastq = [trimmed_r1]
        outputs.fastq2 = [trimmed_r2]
