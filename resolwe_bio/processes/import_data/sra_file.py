"""Import-SRA."""
import glob
from pathlib import Path
from shutil import move

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    FileField,
    FileHtmlField,
    GroupField,
    IntegerField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def run_fastqc(fastqs, output_dir):
    """Run fastQC on given FASTQs.

    :param list fastqs: List of fastqs
    :param str output_dir: Output directory

    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    cmd = Cmd["fastqc"]
    for fastq in fastqs:
        cmd = cmd[fastq]
    cmd = cmd["--extract"]
    cmd = cmd[f"--outdir={str(output_path)}"]
    _, _, stderr = cmd & TEE

    return stderr


class ImportSra(Process):
    """Import reads from SRA.

    Import single or paired-end reads from Sequence Read Archive (SRA)
    via an SRA accession number. SRA stores raw sequencing data and
    alignment information from high-throughput sequencing platforms.

    """

    slug = "import-sra"
    name = "SRA data"
    process_type = "data:sra"
    version = "1.1.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.5.0"}},
        "resources": {"cores": 1, "memory": 1024, "network": True,},
    }
    data_name = "{{ sra_accession|first }}"

    class Input:
        """Input fields to process ImportSra."""

        sra_accession = ListField(StringField(), label="SRA accession(s)")
        show_advanced = BooleanField(label="Show advanced options", default=False)

        class Advanced:
            """Advanced options."""

            min_spot_id = IntegerField(label="Minimum spot ID", required=False)
            max_spot_id = IntegerField(label="Maximum spot ID", required=False)
            min_read_len = IntegerField(label="Minimum read length", required=False)
            clip = BooleanField(label="Clip adapter sequences", default=False)
            aligned = BooleanField(label="Dump only aligned sequences", default=False)
            unaligned = BooleanField(
                label="Dump only unaligned sequences", default=False
            )

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!show_advanced"
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        process_slugs = set()
        for srr in inputs.sra_accession:
            first_read = Cmd["fastq-dump"]("-X", 1, "-Z", "--split-spot", srr)
            num_lines = len(first_read.strip().split("\n"))
            if num_lines == 4:
                process_slug = "import-sra-single"
            else:
                process_slug = "import-sra-paired"
            process_slugs.add(process_slug)

        if len(process_slugs) > 1:
            self.error(
                "All reads must be either single-end or paired-end. Mixing SRA samples of "
                "different types is not allowed."
            )

        process_inputs = {
            "sra_accession": inputs.sra_accession,
            "advanced": {
                "clip": inputs.advanced.clip,
                "aligned": inputs.advanced.aligned,
                "unaligned": inputs.advanced.unaligned,
            },
        }
        if inputs.advanced.min_spot_id:
            process_inputs["advanced"]["min_spot_id"] = inputs.advanced.min_spot_id
        if inputs.advanced.max_spot_id:
            process_inputs["advanced"]["max_spot_id"] = inputs.advanced.max_spot_id
        if inputs.advanced.min_read_len:
            process_inputs["advanced"]["min_read_len"] = inputs.advanced.min_read_len
        if inputs.advanced.min_read_len:
            process_inputs["advanced"]["min_read_len"] = inputs.advanced.min_read_len

        self.run_process(process_slug, process_inputs)


class ImportSraSingle(Process):
    """Import single-end reads from SRA.

    Import single-end reads from Sequence Read Archive (SRA) via an SRA
    accession number. SRA stores raw sequencing data and alignment
    information from high-throughput sequencing platforms.

    """

    slug = "import-sra-single"
    name = "SRA data (single-end)"
    process_type = "data:reads:fastq:single"
    version = "1.1.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.5.0"}},
        "resources": {"cores": 1, "memory": 1024, "network": True,},
    }
    data_name = "{{ sra_accession|first }}"

    class Input:
        """Input fields to process ImportSraSingle."""

        sra_accession = ListField(StringField(), label="SRA accession(s)")
        show_advanced = BooleanField(label="Show advanced options", default=False)

        class Advanced:
            """Advanced options."""

            min_spot_id = IntegerField(label="Minimum spot ID", required=False)
            max_spot_id = IntegerField(label="Maximum spot ID", required=False)
            min_read_len = IntegerField(label="Minimum read length", required=False)
            clip = BooleanField(label="Clip adapter sequences", default=False)
            aligned = BooleanField(label="Dump only aligned sequences", default=False)
            unaligned = BooleanField(
                label="Dump only unaligned sequences", default=False
            )

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!show_advanced"
        )

    class Output:
        """Output fields to process ImportSraSingle."""

        fastq = ListField(FileField(), label="Reads file")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC",)
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run the analysis."""
        fastq_gz = []
        for srr in inputs.sra_accession:
            cmd = Cmd["fastq-dump"]["--gzip"]
            if inputs.advanced.min_spot_id:
                cmd = cmd["--minSpotId"][inputs.advanced.min_spot_id]
            if inputs.advanced.max_spot_id:
                cmd = cmd["--maxSpotId"][inputs.advanced.max_spot_id]
            if inputs.advanced.min_read_len:
                cmd = cmd["--minReadLen"][inputs.advanced.min_read_len]
            if inputs.advanced.clip:
                cmd = cmd["--clip"]
            if inputs.advanced.aligned:
                cmd = cmd["--aligned"]
            if inputs.advanced.unaligned:
                cmd = cmd["--unaligned"]
            return_code, _, _ = cmd[srr] & TEE(retcode=None)
            if return_code:
                self.error(f"Download of {srr} reads with fastq-dump failed.")

            srr_file = Path(f"{srr}.fastq.gz")
            if not srr_file.is_file():
                self.error(f"Download of {srr} reads with fastq-dump failed.")

            fastq_gz.append(str(srr_file))

        if not fastq_gz:
            self.error(
                f"Download of FASTQ files failed for all of the listed accession numbers: {inputs.sra_accession}."
            )
        outputs.fastq = fastq_gz

        print("Postprocessing FastQC...")
        stderr = run_fastqc([fastq_gz], "./fastqc")
        # FastQC writes both progress and errors to stderr and exits with code 0.
        # Catch if file is empty, wrong format... (Failed to process) or
        # if file path does not exist, file cannot be read (Skipping).
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            move(fastqc_zip, ".")

        fastqc = []
        fastqc_url = []
        for fq in fastq_gz:
            reads_name = Path(fq).name.replace(".fastq.gz", "")
            report_dir = Path("fastqc") / Path(f"{reads_name}_fastqc")
            if not report_dir.is_dir():
                continue

            fastqc_zip = Path(f"{reads_name}_fastqc.zip")
            if not fastqc_zip.is_file():
                self.error(f"FastQC failed to produce {fastqc_zip} file.")
            fastqc.append(str(fastqc_zip))

            fastqc_url.append(
                {
                    "file": str(report_dir / "fastqc_report.html"),
                    "refs": [str(report_dir)],
                }
            )
            encoding_file = report_dir / "fastqc_data.txt"
            encoding = Cmd["parse_encoding_type.py"](str(encoding_file))

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                print("Recoding input reads from Phred64 encoding to Phred33 encoding.")
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                Cmd[
                    "TrimmomaticSE",
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ]()
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")

            elif encoding != "Sanger / Illumina 1.9\n":
                self.error(
                    "Only Sanger / Illumina 1.9 / llumina 1.5 / Illumina 1.3 encoding is "
                    "supported."
                )

        outputs.fastqc_url = fastqc_url
        outputs.fastqc_archive = fastqc


class ImportSraPaired(Process):
    """Import paired-end reads from SRA.

    Import paired-end reads from Sequence Read Archive (SRA) via an SRA
    accession number. SRA stores raw sequencing data and alignment
    information from high-throughput sequencing platforms.

    """

    slug = "import-sra-paired"
    name = "SRA data (paired-end)"
    process_type = "data:reads:fastq:paired"
    version = "1.1.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.5.0"}},
        "resources": {"cores": 1, "memory": 1024, "network": True,},
    }
    data_name = "{{ sra_accession|first }}"

    class Input:
        """Input fields to process ImportSraPaired."""

        sra_accession = ListField(StringField(), label="SRA accession(s)")
        show_advanced = BooleanField(label="Show advanced options", default=False)

        class Advanced:
            """Advanced options."""

            min_spot_id = IntegerField(label="Minimum spot ID", required=False)
            max_spot_id = IntegerField(label="Maximum spot ID", required=False)
            min_read_len = IntegerField(label="Minimum read length", required=False)
            clip = BooleanField(label="Clip adapter sequences", default=False)
            aligned = BooleanField(label="Dump only aligned sequences", default=False)
            unaligned = BooleanField(
                label="Dump only unaligned sequences", default=False
            )

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!show_advanced"
        )

    class Output:
        """Output fields to process ImportSraSingle."""

        fastq = ListField(FileField(), label="Reads file (mate 1)")
        fastq2 = ListField(FileField(), label="Reads file (mate 2)")
        fastqc_url = ListField(
            FileHtmlField(), label="Quality control with FastQC (mate 1)",
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Quality control with FastQC (mate 2)",
        )
        fastqc_archive = ListField(
            FileField(), label="Download FastQC archive (mate 1)"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download FastQC archive (mate 2)"
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        fastq_gz_1 = []
        fastq_gz_2 = []
        fastqc_1 = []
        fastqc_2 = []
        fastqc_url_1 = []
        fastqc_url_2 = []
        fastqc_dir = Path("fastqc")
        for srr in inputs.sra_accession:
            cmd = Cmd["fastq-dump"]["--gzip"]["--split-files"]
            if inputs.advanced.min_spot_id:
                cmd = cmd["--minSpotId"][inputs.advanced.min_spot_id]
            if inputs.advanced.max_spot_id:
                cmd = cmd["--maxSpotId"][inputs.advanced.max_spot_id]
            if inputs.advanced.min_read_len:
                cmd = cmd["--minReadLen"][inputs.advanced.min_read_len]
            if inputs.advanced.clip:
                cmd = cmd["--clip"]
            if inputs.advanced.aligned:
                cmd = cmd["--aligned"]
            if inputs.advanced.unaligned:
                cmd = cmd["--unaligned"]
            return_code, _, _ = cmd[srr] & TEE(retcode=None)
            if return_code:
                self.error(f"Download of {srr} reads with fastq-dump failed.")

            mate1 = Path(f"{srr}_1.fastq.gz")
            mate2 = Path(f"{srr}_2.fastq.gz")
            if not mate1.exists() or not mate2.exists():
                self.error(f"Download of {srr} reads with fastq-dump failed.")
            fastq_gz_1.append(str(mate1))
            fastq_gz_2.append(str(mate2))
            fastqc_1.append(f"{srr}_1_fastqc.zip")
            fastqc_2.append(f"{srr}_2_fastqc.zip")
            fastqc_path_1 = fastqc_dir / f"{srr}_1_fastqc"
            fastqc_url_1.append(
                {
                    "file": str(fastqc_path_1 / "fastqc_report.html"),
                    "refs": [str(fastqc_path_1)],
                }
            )
            fastqc_path_2 = fastqc_dir / f"{srr}_2_fastqc"
            fastqc_url_2.append(
                {
                    "file": str(fastqc_path_2 / "fastqc_report.html"),
                    "refs": [str(fastqc_path_2)],
                }
            )

        outputs.fastq = fastq_gz_1
        outputs.fastq2 = fastq_gz_2

        print("Postprocessing FastQC...")
        stderr = run_fastqc(fastq_gz_1 + fastq_gz_2, "./fastqc")
        # FastQC writes both progress and errors to stderr and exits with code 0.
        # Catch if file is empty, wrong format... (Failed to process) or
        # if file path does not exist, file cannot be read (Skipping).
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            move(fastqc_zip, ".")

        for report_dir in fastqc_dir.iterdir():
            if not report_dir.is_dir():
                continue

            reads_name = report_dir.name.replace("_fastqc", "")
            encoding_file = report_dir / "fastqc_data.txt"
            encoding = Cmd["parse_encoding_type.py"](str(encoding_file))

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                print("Recoding input reads from Phred64 encoding to Phred33 encoding.")
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                Cmd[
                    "TrimmomaticSE",
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ]()
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")

            elif encoding != "Sanger / Illumina 1.9\n":
                self.error(
                    "Only Sanger / Illumina 1.9 / llumina 1.5 / Illumina 1.3 encoding is "
                    "supported."
                )

        outputs.fastqc_url = fastqc_url_1
        outputs.fastqc_url2 = fastqc_url_2
        outputs.fastqc_archive = fastqc_1
        outputs.fastqc_archive2 = fastqc_2
