"""Import reads (FASTQ)."""

import collections
import glob
import gzip
import shutil
from pathlib import Path

import dnaio
from dnaio.exceptions import FastqFormatError, FileFormatError
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
)

SUPPORTED_EXTENSIONS = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
)


def check_file(infile):
    """Check if the input file exists and has correct extensions."""
    fq_file = Path(infile)
    if not fq_file.is_file():
        message = "Input file {} does not exist".format(fq_file.name)
        return message

    if not fq_file.name.lower().endswith(SUPPORTED_EXTENSIONS):
        message = (
            "Unrecognized file name extension in file {}. "
            "Supported file name extensions are {}.".format(
                fq_file.name, SUPPORTED_EXTENSIONS
            )
        )
        return message

    message = "Correct input file."
    return message


def validate_fastq(fq, fq2=None):
    """Validate FASTQ files."""

    input_fastq = fq

    if fq2:
        input_fastq = fq + fq2

    # Reduce the probability of uploading the FASTQ files with the same
    # content multiple times (as multiple lanes or mates).
    if len(set(input_fastq)) != len(input_fastq):
        seen_files = [
            item
            for item, count in collections.Counter(input_fastq).items()
            if count > 1
        ]
        message = "Non-unique input file names detected: {}.".format(seen_files)
        return message

    if fq2 and len(fq) != len(fq2):
        message = (
            "The number of mate-pair files in split-lane samples must match. "
            "{} and {} input files were given for the -fq and -fq2 inputs, "
            "respectively.".format(len(fq), len(fq2))
        )
        return message

    if fq2:
        for mate1, mate2 in zip(fq, fq2):
            try:
                with gzip.open(mate1) as mate1, gzip.open(mate2) as mate2:
                    paired_reads = dnaio.open(mate1, file2=mate2, fileformat="fastq")
                    if not any(paired_reads):
                        message = "Mate-pair files {} and {} contain no read sequences.".format(
                            mate1.name, mate2.name
                        )
                        return message

                    else:
                        for read in paired_reads:
                            continue
                        message = (
                            "Successfully validated mate-pair files {} and {}.".format(
                                mate1.name, mate2.name
                            )
                        )
                        return message

            except (FastqFormatError, FileFormatError) as dnaio_error:
                message = "Format error in mate-pairs {} and {}. {}".format(
                    mate1.name, mate2.name, str(dnaio_error)
                )
                return message

    else:
        for fq in fq:
            fq = Path(fq)
            try:
                with gzip.open(fq) as read:
                    reads = dnaio.open(read, fileformat="fastq")
                    if not any(reads):
                        message = "Input file {} contains no read sequences.".format(
                            fq.name
                        )
                        return message

                    else:
                        for read in reads:
                            continue
                        message = "Successfully validated reads file {}.".format(
                            fq.name
                        )
                        return message

            except (FastqFormatError, FileFormatError) as dnaio_error:
                message = "Error in file {}. {}".format(fq.name, str(dnaio_error))
                return message


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


def parse_encoding_type(report_file):
    """Parse encoding type."""
    encoding = ""
    with open(report_file) as report:
        for line in report:
            if line.startswith("Encoding"):
                encoding = line.strip().split("\t")[1]
                break
        if encoding != "":
            return encoding
        else:
            return "Unknown"


def replace_extension(infile):
    """Replace extensions of file."""
    extensions = "".join(Path(str(infile)).suffixes[-2:])
    new_ext = ".fastq.gz"
    outfile = str(infile).replace(extensions, new_ext)
    return outfile


class UploadFastqSingle(Process):
    """Import single-end reads in FASTQ format.

    Import single-end reads in FASTQ format, which is a text-based format for
    storing both a biological sequence (usually nucleotide sequence) and its
    corresponding quality scores.
    """

    slug = "upload-fastq-single"
    name = "FASTQ file (single-end)"
    process_type = "data:reads:fastq:single"
    version = "2.6.0"
    category = "Import"
    data_name = '{{ src.0.file|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 1,
            "network": True,
        },
    }

    class Input:
        """Input fields to process UploadFastqSingle."""

        src = ListField(
            FileField(),
            label="Reads",
            description="Sequencing reads in FASTQ format. "
            "Supported extensions: .fastq.gz (preferred), .fq.* or .fastq.*",
        )
        merge_lanes = BooleanField(
            label="Merge lanes",
            default=False,
            description="Merge sample data split into multiple sequencing "
            "lanes into a single FASTQ file.",
        )

    class Output:
        """Output fields to process UploadFastqSingle."""

        fastq = ListField(FileField(), label="Reads file")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC")
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run upload."""

        fastqgz = []

        for read in inputs.src:
            read_imported = read.import_file(imported_format="compressed")
            stderr = check_file(infile=read_imported)
            if "Correct input file." not in stderr:
                self.error(stderr)
            renamed_reads = replace_extension(infile=read_imported)
            Path(read_imported).rename(renamed_reads)
            fastqgz.append(renamed_reads)

        stderr = validate_fastq(fq=fastqgz)
        if "Successfully validated reads" not in stderr:
            self.error(stderr)

        if inputs.merge_lanes:
            first_read = fastqgz[0][:-9]
            fastqz = f"{first_read}_merged.fastq.gz"
            with gzip.open(fastqz, "wb") as outfile:
                for read in fastqgz:
                    with gzip.open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            fastqgz = [fastqz]

        stderr = run_fastqc([fastqgz], "./fastqc")
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            shutil.move(fastqc_zip, ".")

        fastqc = []
        fastqc_url = []
        for fq in fastqgz:
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
            encoding = parse_encoding_type(report_file=encoding_file)

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                self.info(
                    "Recoding input reads from Phred64 encoding to Phred33 encoding."
                )
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                return_code, _, stderr = Cmd["TrimmomaticSE"][
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ] & TEE(retcode=None)
                if return_code:
                    print(stderr)
                    self.error("Error while running TrimmomaticSE.")
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")
            elif encoding != "Sanger / Illumina 1.9":
                self.error(
                    "Only Sanger / Illumina 1.9 / lllumina 1.5 / Illumina 1.3 encoding is supported."
                )

        outputs.fastq = fastqgz
        outputs.fastqc_url = fastqc_url
        outputs.fastqc_archive = fastqc


class UploadFastqPaired(Process):
    """Import paired-end reads in FASTQ format.

    Import paired-end reads in FASTQ format, which is a text-based format for
    storing both a biological sequence (usually nucleotide sequence) and its
    corresponding quality scores.
    """

    slug = "upload-fastq-paired"
    name = "FASTQ file (paired-end)"
    process_type = "data:reads:fastq:paired"
    version = "2.6.0"
    category = "Import"
    data_name = '{{ src1.0.file|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 1,
            "network": True,
        },
    }

    class Input:
        """Input fields to process UploadFastqPaired."""

        src1 = ListField(
            FileField(),
            label="Mate1",
            description="Sequencing reads in FASTQ format. "
            "Supported extensions: .fastq.gz (preferred), .fq.* or .fastq.*",
        )
        src2 = ListField(
            FileField(),
            label="Mate2",
            description="Sequencing reads in FASTQ format. "
            "Supported extensions: .fastq.gz (preferred), .fq.* or .fastq.*",
        )
        merge_lanes = BooleanField(
            label="Merge lanes",
            default=False,
            description="Merge sample data split into multiple sequencing "
            "lanes into a single FASTQ file.",
        )

    class Output:
        """Output fields to process UploadFastqPaired."""

        fastq = ListField(FileField(), label="Reads file (mate 1)")
        fastq2 = ListField(FileField(), label="Reads file (mate 2)")
        fastqc_url = ListField(
            FileHtmlField(), label="Quality control with FastQC (Upstream)"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Quality control with FastQC (Downstream)"
        )
        fastqc_archive = ListField(
            FileField(), label="Download FastQC archive (Upstream)"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download FastQC archive (Downstream)"
        )

    def run(self, inputs, outputs):
        """Run upload."""

        mate1_fastqgz = []
        mate1_fastqc = []
        mate1_fastqc_url = []

        mate2_fastqgz = []
        mate2_fastqc = []
        mate2_fastqc_url = []

        rep_dir = Path("fastqc")

        for read in inputs.src1:
            read_imported = read.import_file(imported_format="compressed")
            stderr = check_file(infile=read_imported)
            if "Correct input file." not in stderr:
                self.error(stderr)
            renamed_reads = replace_extension(infile=read_imported)
            Path(read_imported).rename(renamed_reads)
            name_mate1 = renamed_reads[:-9]
            mate1_fastqgz.append(renamed_reads)
            mate1_fastqc.append(f"{name_mate1}_fastqc.zip")
            mate1_fastqc_url.append(
                {
                    "file": str(
                        rep_dir / f"{name_mate1}_fastqc" / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{name_mate1}_fastqc")],
                }
            )
        for read in inputs.src2:
            read_imported = read.import_file(imported_format="compressed")
            stderr = check_file(infile=read_imported)
            if "Correct input file." not in stderr:
                self.error(stderr)
            renamed_reads = replace_extension(infile=read_imported)
            Path(read_imported).rename(renamed_reads)
            name_mate2 = renamed_reads[:-9]
            mate2_fastqgz.append(renamed_reads)
            mate2_fastqc.append(f"{name_mate2}_fastqc.zip")
            mate2_fastqc_url.append(
                {
                    "file": str(
                        rep_dir / f"{name_mate2}_fastqc" / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{name_mate2}_fastqc")],
                }
            )

        stderr = validate_fastq(fq=mate1_fastqgz, fq2=mate2_fastqgz)
        if "Successfully validated mate-pair" not in stderr:
            self.error(stderr)

        if inputs.merge_lanes:
            mate1_first_lane = mate1_fastqgz[0][:-9]
            fastqz_1 = f"{mate1_first_lane}_merged.fastq.gz"
            with gzip.open(fastqz_1, "wb") as outfile:
                for read in mate1_fastqgz:
                    with gzip.open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            mate1_fastqgz = [f"{mate1_first_lane}_merged.fastq.gz"]
            mate1_fastqc = [f"{mate1_first_lane}_merged_fastqc.zip"]
            mate1_fastqc_url = [
                {
                    "file": str(
                        rep_dir
                        / f"{mate1_first_lane}_merged_fastqc"
                        / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{mate1_first_lane}_merged_fastqc")],
                }
            ]

            mate2_first_lane = mate2_fastqgz[0][:-9]
            fastqz_2 = f"{mate2_first_lane}_merged.fastq.gz"
            with gzip.open(fastqz_2, "wb") as outfile:
                for read in mate2_fastqgz:
                    with gzip.open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            mate2_fastqgz = [f"{mate2_first_lane}_merged.fastq.gz"]
            mate2_fastqc = [f"{mate2_first_lane}_merged_fastqc.zip"]
            mate2_fastqc_url = [
                {
                    "file": str(
                        rep_dir
                        / f"{mate2_first_lane}_merged_fastqc"
                        / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{mate2_first_lane}_merged_fastqc")],
                }
            ]

        stderr = run_fastqc(mate1_fastqgz + mate2_fastqgz, "./fastqc")
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            shutil.move(fastqc_zip, ".")

        for report_dir in Path("fastqc").iterdir():
            if not report_dir.is_dir():
                continue

            reads_name = report_dir.name.replace("_fastqc", "")
            encoding_file = report_dir / "fastqc_data.txt"
            encoding = parse_encoding_type(report_file=encoding_file)

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                print("Recoding input reads from Phred64 encoding to Phred33 encoding.")
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                return_code, _, stderr = Cmd["TrimmomaticSE"][
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ] & TEE(retcode=None)
                if return_code:
                    print(stderr)
                    self.error("Recoding of input reads failed.")
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")

            elif encoding != "Sanger / Illumina 1.9":
                self.error(
                    "Only Sanger / Illumina 1.9 / lllumina 1.5 / Illumina 1.3 encoding is "
                    "supported."
                )

        outputs.fastq = mate1_fastqgz
        outputs.fastq2 = mate2_fastqgz
        outputs.fastqc_url = mate1_fastqc_url
        outputs.fastqc_url2 = mate2_fastqc_url
        outputs.fastqc_archive = mate1_fastqc
        outputs.fastqc_archive2 = mate2_fastqc


class FilesToFastqSingle(Process):
    """Convert FASTQ files to single-end reads."""

    slug = "files-to-fastq-single"
    name = "Convert files to reads (single-end)"
    process_type = "data:reads:fastq:single"
    version = "1.6.0"
    category = "Import"
    data_name = "Files to FASTQ single-end ({{ (src|first).file.file }})"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
    }

    class Input:
        """Input fields to process FilesToFastqSingle."""

        src = ListField(
            DataField("file"),
            label="Reads",
            description="Sequencing reads in FASTQ format",
        )
        merge_lanes = BooleanField(
            label="Merge lanes",
            default=False,
            description="Merge sample data split into multiple sequencing "
            "lanes into a single FASTQ file.",
        )

    class Output:
        """Output fields to process FilesToFastqSingle."""

        fastq = ListField(FileField(), label="Reads file")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC")
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run upload."""

        fastqgz = []

        for read in inputs.src:
            read_name = Path(read.output.file.path).name
            shutil.copy(read.output.file.path, ".")
            stderr = check_file(infile=read_name)
            if "Correct input file." not in stderr:
                self.error(stderr)
            if not read_name.endswith(".gz"):
                Cmd["pigz"][read_name]()
                renamed_reads = replace_extension(infile=f"{read_name}.gz")
                Path(f"{read_name}.gz").rename(renamed_reads)
                fastqgz.append(renamed_reads)
            else:
                renamed_reads = replace_extension(infile=read_name)
                Path(read_name).rename(renamed_reads)
                fastqgz.append(renamed_reads)

        stderr = validate_fastq(fq=fastqgz)
        if "Successfully validated reads" not in stderr:
            self.error(stderr)

        if inputs.merge_lanes:
            first_read = fastqgz[0][:-9]
            fastqz = f"{first_read}_merged.fastq.gz"
            with open(fastqz, "wb") as outfile:
                for read in fastqgz:
                    with open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            fastqgz = [fastqz]

        stderr = run_fastqc([fastqgz], "./fastqc")
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            shutil.move(fastqc_zip, ".")

        fastqc = []
        fastqc_url = []
        for fq in fastqgz:
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
            encoding = parse_encoding_type(report_file=encoding_file)
            self.info(encoding)

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                self.info(
                    "Recoding input reads from Phred64 encoding to Phred33 encoding."
                )
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                return_code, _, stderr = Cmd["TrimmomaticSE"][
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ] & TEE(retcode=None)
                if return_code:
                    print(stderr)
                    self.error("Error while running TrimmomaticSE.")
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")
            elif encoding != "Sanger / Illumina 1.9":
                self.error(
                    "Only Sanger / Illumina 1.9 / lllumina 1.5 / Illumina 1.3 encoding is supported."
                )

        outputs.fastq = fastqgz
        outputs.fastqc_url = fastqc_url
        outputs.fastqc_archive = fastqc


class FilesToFastqPaired(Process):
    """Convert FASTQ files to paired-end reads."""

    slug = "files-to-fastq-paired"
    name = "Convert files to reads (paired-end)"
    process_type = "data:reads:fastq:paired"
    version = "1.6.0"
    category = "Import"
    data_name = "Files to FASTQ paired-end ({{ (src1|first).file.file }}, {{(src2|first).file.file}})"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
    }

    class Input:
        """Input fields to process FilesToFastqPaired."""

        src1 = ListField(
            DataField("file"),
            label="Mate1",
        )
        src2 = ListField(
            DataField("file"),
            label="Mate2",
        )
        merge_lanes = BooleanField(
            label="Merge lanes",
            default=False,
            description="Merge sample data split into multiple sequencing "
            "lanes into a single FASTQ file.",
        )

    class Output:
        """Output fields to process FilesToFastqPaired."""

        fastq = ListField(FileField(), label="Reads file (mate 1)")
        fastq2 = ListField(FileField(), label="Reads file (mate 2)")
        fastqc_url = ListField(
            FileHtmlField(), label="Quality control with FastQC (Upstream)"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Quality control with FastQC (Downstream)"
        )
        fastqc_archive = ListField(
            FileField(), label="Download FasQC archive (Upstream)"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download FasQC archive (Downstream)"
        )

    def run(self, inputs, outputs):
        """Run upload."""

        mate1_fastqgz = []
        mate1_fastqc = []
        mate1_fastqc_url = []

        mate2_fastqgz = []
        mate2_fastqc = []
        mate2_fastqc_url = []

        rep_dir = Path("fastqc")

        for read in inputs.src1:
            read_name = Path(read.output.file.path).name
            shutil.copy(read.output.file.path, ".")
            stderr = check_file(infile=read_name)
            if "Correct input file." not in stderr:
                self.error(stderr)
            if not read_name.endswith(".gz"):
                Cmd["pigz"][read_name]()
                renamed_reads = replace_extension(infile=f"{read_name}.gz")
                Path(f"{read_name}.gz").rename(renamed_reads)
                name_mate1 = read_name[:-9]
            else:
                renamed_reads = replace_extension(infile=read_name)
                Path(read_name).rename(renamed_reads)
                name_mate1 = read_name[:-9]
            mate1_fastqgz.append(f"{name_mate1}.fastq.gz")
            mate1_fastqc.append(f"{name_mate1}_fastqc.zip")
            mate1_fastqc_url.append(
                {
                    "file": str(
                        rep_dir / f"{name_mate1}_fastqc" / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{name_mate1}_fastqc")],
                }
            )
        for read in inputs.src2:
            read_name = Path(read.output.file.path).name
            shutil.copy(read.output.file.path, ".")
            stderr = check_file(infile=read_name)
            if "Correct input file." not in stderr:
                self.error(stderr)
            if not read_name.endswith(".gz"):
                Cmd["pigz"][read_name]()
                renamed_reads = replace_extension(infile=f"{read_name}.gz")
                Path(f"{read_name}.gz").rename(renamed_reads)
                name_mate2 = read_name[:-9]
            else:
                renamed_reads = replace_extension(infile=read_name)
                Path(read_name).rename(renamed_reads)
                name_mate2 = read_name[:-9]
            mate2_fastqgz.append(f"{name_mate2}.fastq.gz")
            mate2_fastqc.append(f"{name_mate2}_fastqc.zip")
            mate2_fastqc_url.append(
                {
                    "file": str(
                        rep_dir / f"{name_mate2}_fastqc" / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{name_mate2}_fastqc")],
                }
            )

        stderr = validate_fastq(fq=mate1_fastqgz, fq2=mate2_fastqgz)
        if "Successfully validated mate-pair" not in stderr:
            self.error(stderr)

        if inputs.merge_lanes:
            mate1_first_lane = mate1_fastqgz[0][:-9]
            fastqz_1 = f"{mate1_first_lane}_merged.fastq.gz"
            with gzip.open(fastqz_1, "wb") as outfile:
                for read in mate1_fastqgz:
                    with gzip.open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            mate1_fastqgz = [f"{mate1_first_lane}_merged.fastq.gz"]
            mate1_fastqc = [f"{mate1_first_lane}_merged_fastqc.zip"]
            mate1_fastqc_url = [
                {
                    "file": str(
                        rep_dir
                        / f"{mate1_first_lane}_merged_fastqc"
                        / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{mate1_first_lane}_merged_fastqc")],
                }
            ]

            mate2_first_lane = mate2_fastqgz[0][:-9]
            fastqz_2 = f"{mate2_first_lane}_merged.fastq.gz"
            with gzip.open(fastqz_2, "wb") as outfile:
                for read in mate2_fastqgz:
                    with gzip.open(read, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            mate2_fastqgz = [f"{mate2_first_lane}_merged.fastq.gz"]
            mate2_fastqc = [f"{mate2_first_lane}_merged_fastqc.zip"]
            mate2_fastqc_url = [
                {
                    "file": str(
                        rep_dir
                        / f"{mate2_first_lane}_merged_fastqc"
                        / "fastqc_report.html"
                    ),
                    "refs": [str(rep_dir / f"{mate2_first_lane}_merged_fastqc")],
                }
            ]

        stderr = run_fastqc(mate1_fastqgz + mate2_fastqgz, "./fastqc")
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in glob.glob("fastqc/*_fastqc.zip"):
            shutil.move(fastqc_zip, ".")

        for report_dir in Path("fastqc").iterdir():
            if not report_dir.is_dir():
                continue

            reads_name = report_dir.name.replace("_fastqc", "")
            encoding_file = report_dir / "fastqc_data.txt"
            encoding = parse_encoding_type(report_file=encoding_file)

            if encoding == "Illumina 1.5" or encoding == "Illumina 1.3":
                print("Recoding input reads from Phred64 encoding to Phred33 encoding.")
                Path(f"{reads_name}.fastq.gz").rename("input_reads.fastq.gz")
                return_code, _, stderr = Cmd["TrimmomaticSE"][
                    "-phred64",
                    "input_reads.fastq.gz",
                    "reformated.fastq.gz",
                    "TOPHRED33",
                ] & TEE(retcode=None)
                if return_code:
                    print(stderr)
                    self.error("Recoding of input reads failed.")
                Path("reformated.fastq.gz").rename(f"{reads_name}.fastq.gz")

            elif encoding != "Sanger / Illumina 1.9":
                self.error(
                    "Only Sanger / Illumina 1.9 / lllumina 1.5 / Illumina 1.3 encoding is "
                    "supported."
                )

        outputs.fastq = mate1_fastqgz
        outputs.fastq2 = mate2_fastqgz
        outputs.fastqc_url = mate1_fastqc_url
        outputs.fastqc_url2 = mate2_fastqc_url
        outputs.fastqc_archive = mate1_fastqc
        outputs.fastqc_archive2 = mate2_fastqc
