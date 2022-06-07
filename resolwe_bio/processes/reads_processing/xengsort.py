"""Xenograft read sorting with Xengsort."""

import gzip
import re
import shutil
from math import ceil
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)

HASH_FUNCTION = "linear945:linear9123641:linear349341847"
READ_TYPES = ["host", "graft", "both", "neither", "ambiguous"]


def concatenate_files(filenames, out_fname, decompress=False):
    """Concatenate files and optionally decompress them."""
    handler = gzip.open if decompress else open
    with open(out_fname, "wb") as out_handle:
        for fname in filenames:
            with handler(fname, "rb") as handle:
                shutil.copyfileobj(
                    fsrc=handle, fdst=out_handle, length=1024 * 1024 * 10
                )


def create_filename(basename, suffix):
    """Create a secure version of a filename."""
    # Keep only alphanumeric characters and underscores.
    basename = re.sub(pattern=r"[\W]", repl="_", string=basename)
    # Replace consecutive underscores with a single one.
    basename = re.sub(pattern=r"(_)\1+", repl=r"\1", string=basename)
    return f"{basename}.{suffix}"


def concatenate_reads(filenames, out_fasta, error):
    """Concatenate reads in FASTQ files."""
    try:
        concatenate_files(filenames=filenames, out_fname=out_fasta, decompress=True)
    except Exception as e:
        error(
            f"Failed to concatenate FASTQ files for {out_fasta.name}. The error was: {repr(e)}"
        )


class XengsortIndex(Process):
    """Build an index for sorting xenograft reads with Xengsort.

    Xengsort is an alignment free method for sorting reads from
    xenograft experiments. Description of the method and evaluation on
    several datasets is provided in the
    [article](https://doi.org/10.1186/s13015-021-00181-w).
    """

    slug = "xengsort-index"
    name = "Xengsort index"
    process_type = "data:xengsort:index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/xengsort:1.0.0"},
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
        },
    }
    category = "Xenograft processing"
    data_name = "Xengsort index"
    version = "1.0.1"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    class Input:
        """Input fields to XengsortIndex process."""

        graft_refs = ListField(
            DataField("seq:nucleotide"),
            label="Graft reference sequences (nucleotide FASTA)",
        )

        host_refs = ListField(
            DataField("seq:nucleotide"),
            label="Host reference sequences (nucleotide FASTA)",
        )

        n_kmer = IntegerField(
            label="Number of distinct k-mers [--nobjects]",
            required=False,
            description="The number of k-mers that will be stored in the hash table. This depends "
            "on the used reference genomes and must be estimated beforehand. If the number of "
            "distinct k-mers is known beforehand it should be specified. For all 25-mers in the "
            "human and mouse genome and transcriptome, this number is roughly 4,500,000,000. "
            "If this is not set, the number is estimated with ntCard tool and increased by two "
            "percent to account for errors.",
        )

        class Advanced:
            """Advanced options."""

            kmer_size = IntegerField(label="k-mer size [--kmersize]", default=25)

            aligned_cache = BooleanField(
                label="Use power-of-two aligned pages [--aligned]",
                default=False,
                description="Indicates whether each bucket should consume a number of bits that "
                "is a power of 2. Using --aligned ensures that each bucket stays within the same "
                "cache line, but may waste space (padding bits), yielding faster speed but "
                "larger space requirements. By default no bits are used for padding and buckets "
                "may cross cache line boundaries [--unaligned]. This is slightly slower, but may "
                "save a little or a lot of space.",
            )

            fixed_hashing = BooleanField(
                label="Use fixed hash function [--hashfunctions]",
                default=True,
                description="Hash function used to store the key-value pairs is defined by "
                "--hashfunction parameter. With this option selected a fixed hash function "
                "(linear945:linear9123641:linear349341847) is used. When this is not selected a "
                "different random functions are chosen each time. It is recommended to have them "
                "chosen randomly unless you need strictly reproducible behavior.",
            )

            page_size = IntegerField(
                label="Number of elements stored in one bucket (or page) [--pagesize]",
                default=4,
            )

            fill = FloatField(
                label="Fill rate of the hash table [--fill]",
                range=[0.0, 1.0],
                default=0.88,
                description="This determines the desired fill rate or load factor of the hash "
                "table. It should be set between 0.0 and 1.0. It is beneficial to leave part of "
                "the hash table empty for faster lookups. Together with the number of distinct "
                "k-mers [--nobjects], the number of slots in the table is calculated as "
                "ceil(nobjects/fill).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to process XengsortIndex."""

        index = FileField(label="Xengsort index")
        stats = FileField(label="Xengsort statistics")
        graft_species = StringField(label="Graft species")
        graft_build = StringField(label="Graft build")
        host_species = StringField(label="Host species")
        host_build = StringField(label="Host build")

    def run(self, inputs, outputs):
        """Run analysis."""

        graft_species = inputs.graft_refs[0].output.species
        graft_build = inputs.graft_refs[0].output.build
        if any(ref.output.species != graft_species for ref in inputs.graft_refs):
            self.error(
                "All input graft reference sequences must be from the same species."
            )

        if any(ref.output.build != graft_build for ref in inputs.graft_refs):
            self.error("All input graft reference sequences must share the same build.")

        host_species = inputs.host_refs[0].output.species
        host_build = inputs.host_refs[0].output.build
        if any(ref.output.species != host_species for ref in inputs.host_refs):
            self.error(
                "All input host reference sequences must be from the same species."
            )

        if any(ref.output.build != host_build for ref in inputs.host_refs):
            self.error("All input host reference sequences must share the same build.")

        outputs.graft_species = graft_species
        outputs.graft_build = graft_build
        outputs.host_species = host_species
        outputs.host_build = host_build

        graft_fasta = Path("graft.fasta")
        concatenate_files(
            filenames=[ref.output.fasta.path for ref in inputs.graft_refs],
            out_fname=graft_fasta,
        )

        host_fasta = Path("host.fasta")
        concatenate_files(
            filenames=[ref.output.fasta.path for ref in inputs.host_refs],
            out_fname=host_fasta,
        )

        if inputs.n_kmer:
            n_kmer = inputs.n_kmer
        else:
            nthll_params = [
                f"--threads={self.requirements.resources.cores}",
                f"--kmer={inputs.advanced.kmer_size}",
                graft_fasta,
                host_fasta,
            ]

            return_code, stdout, stderr = Cmd["nthll"][nthll_params] & TEE(retcode=None)

            if return_code:
                print(stderr)
                self.error("Failed to estimate the number of distinct k-mers.")

            # Example standard output is "F0, Exp# of distnt kmers(k=25): 47270"
            stdout_search = re.search(
                pattern=r"^F0, Exp# of distnt kmers\(k=\d+\): (\d+)$",
                string=stdout,
                flags=re.MULTILINE,
            )

            try:
                n_kmer = int(stdout_search.group(1))
            except (AttributeError, ValueError, IndexError) as e:
                print(e.message)
                self.error("Failed to extract the number of distinct k-mers.")

            # Account for possible errors in estimation by increasing
            # the number of distinct k-mers by 2 percent.
            n_kmer = ceil(n_kmer * 1.02)

        index_file = create_filename(
            basename="_".join([graft_species, graft_build, host_species, host_build]),
            suffix="h5",
        )

        index_params = [
            index_file,
            "--host",
            host_fasta,
            "--graft",
            graft_fasta,
            "--kmersize",
            inputs.advanced.kmer_size,
            "--shortcutbits",
            0,
            "--nobjects",
            n_kmer,
            "--type",
            "3FCVbb",
            "--aligned" if inputs.advanced.aligned_cache else "--unaligned",
            "--hashfunctions",
            HASH_FUNCTION if inputs.advanced.fixed_hashing else "random",
            "--pagesize",
            inputs.advanced.page_size,
            "--fill",
            inputs.advanced.fill,
            "--threads",
            self.requirements.resources.cores,
        ]

        return_code, stdout, stderr = Cmd["xengsort"]["index"][index_params] & TEE(
            retcode=None
        )

        if return_code:
            print(stderr)
            self.error("Failed to calculate Xengsort index.")

        stats_file = "index_stats.txt"
        with open(stats_file, "w") as stats_handle:
            stats_handle.writelines(stdout)

        outputs.index = index_file
        outputs.stats = stats_file


class XengsortClassify(Process):
    """Classify xenograft reads with Xengsort.

    Xengsort is an alignment free method for sorting reads from
    xenograft experiments. It classifies sequencing reads into five
    categories based on their origin: host, graft, both, neither, and
    ambiguous. Categories “host” and “graft” are for reads that can be
    clearly assigned to one of the species. Category “both” is for reads
    that match equally well to both references. Category “neither” is
    for reads that contain many k-mers that cannot be found in the
    key-value store; these could point to technical problems (primer
    dimers) or contamination of the sample with other species. Finally,
    category “ambiguous” is for reads that provide conflicting
    information. Such reads should not usually be seen; they could
    result from PCR hybrids between host and graft during library
    preparation.

    Description of the method and evaluation on several
    datasets is provided in the
    [article](https://doi.org/10.1186/s13015-021-00181-w).
    """

    slug = "xengsort-classify"
    name = "Xengsort classify"
    process_type = "data:xengsort:classification"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/xengsort:1.0.0"},
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "Xenograft processing"
    data_name = "{{ reads|name|default('?') }}"
    version = "1.0.0"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    class Input:
        """Input fields to XengsortClassify process."""

        reads = DataField("reads:fastq", label="Reads")
        index = DataField("xengsort:index", label="Xengsort genome index")

        upload_reads = StringField(
            label="Select reads to upload",
            description="All read categories are returned in this process but only the ones "
            "selected are uploaded as separate FASTQ files. This should be used for categories "
            "of reads that will be used in further analyses.",
            choices=[
                ("none", "none"),
                ("all", "all"),
                ("graft", "graft"),
                ("graft, both", "graft, both"),
                ("graft, host", "graft, host"),
                ("graft, host, both", "graft, host, both"),
            ],
            default="none",
        )

        merge_both = BooleanField(
            label="Upload merged graft and both reads",
            description="Merge graft reads with the reads that can originate from both genomes "
            "and upload it as graft reads. In any workflow, the latter reads, classified as both "
            "may pose problems, because one may not be able to decide on the species of origin "
            "due to ultra-conserved regions between species.",
            hidden="upload_reads == 'none'",
            default=False,
        )

        class Advanced:
            """Advanced options."""

            chunksize = FloatField(
                label="Chunk size in MB [--chunksize]",
                default=16.0,
                description="Controll the memory usage by setting chunk size per thread.",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to XengsortClassify process."""

        stats = FileField(label="Xengsort classification statistics")
        host1 = FileField(label="Host reads (mate 1)")
        host2 = FileField(label="Host reads (mate 2)", required=False)
        graft1 = FileField(label="Graft reads (mate 1)")
        graft2 = FileField(label="Graft reads (mate 2)", required=False)
        both1 = FileField(label="Both reads (mate 1)")
        both2 = FileField(label="Both reads (mate 2)", required=False)
        neither1 = FileField(label="Neither reads (mate 1)")
        neither2 = FileField(label="Neither reads (mate 2)", required=False)
        ambiguous1 = FileField(label="Ambiguous reads (mate 1)")
        ambiguous2 = FileField(label="Ambiguous reads (mate 2)", required=False)
        graft_species = StringField(label="Graft species")
        graft_build = StringField(label="Graft build")
        host_species = StringField(label="Host species")
        host_build = StringField(label="Host build")

    def run(self, inputs, outputs):
        """Run analysis."""

        concatenated_r1 = "mate_1.fastq"
        concatenate_reads(
            filenames=[fastq.path for fastq in inputs.reads.output.fastq],
            out_fasta=concatenated_r1,
            error=self.error,
        )

        is_paired = inputs.reads.type.startswith("data:reads:fastq:paired:")
        if is_paired:
            concatenated_r2 = "mate_2.fastq"
            concatenate_reads(
                filenames=[fastq.path for fastq in inputs.reads.output.fastq2],
                out_fasta=concatenated_r2,
                error=self.error,
            )

        if is_paired:
            classify_params = ["--fastq", concatenated_r1, "--pairs", concatenated_r2]
        else:
            classify_params = ["--fastq", concatenated_r1]

        fastq_file = Path(inputs.reads.output.fastq[0].path).name
        assert fastq_file.endswith(".fastq.gz")
        name = fastq_file[:-9]

        classify_params.extend(
            [
                "--index",
                inputs.index.output.index.path,
                "--threads",
                self.requirements.resources.cores,
                "--prefix",
                name,
                "--prefetchlevel",
                0,
                "--chunksize",
                inputs.advanced.chunksize,
            ]
        )

        return_code, stdout, stderr = Cmd["xengsort"]["classify"][
            classify_params
        ] & TEE(retcode=None)

        if return_code:
            print(stderr)
            self.error("Failed to classify reads with Xengsort.")

        stats_file = "classification_stats.txt"
        with open(stats_file, "w") as stats_handle:
            stats_handle.writelines(stdout)

        outputs.stats = stats_file

        output_files = {}
        for read_type in READ_TYPES:
            mates = [1, 2] if is_paired else [1]
            for mate in mates:
                if is_paired:
                    output_file = Path(f"{name}-{read_type}.{mate}.fq")
                else:
                    output_file = Path(f"{name}-{read_type}.fq")

                if output_file.is_file():
                    output_file = output_file.rename(output_file.with_suffix(".fastq"))
                    return_code, _, _ = Cmd["pigz"][output_file] & TEE(retcode=None)
                    if return_code:
                        self.error(f"Compression of {read_type} reads failed.")

                    output_files[f"{read_type}{mate}"] = f"{output_file}.gz"

        for output_key, output_file in output_files.items():
            setattr(outputs, output_key, output_file)

        outputs.graft_species = inputs.index.output.graft_species
        outputs.graft_build = inputs.index.output.graft_build
        outputs.host_species = inputs.index.output.host_species
        outputs.host_build = inputs.index.output.host_build

        if inputs.merge_both:
            merged1 = f"{name}-graft-both{'.1' if is_paired else ''}.fastq.gz"
            concatenate_files(
                filenames=[output_files["graft1"], output_files["both1"]],
                out_fname=merged1,
            )
            output_files["graft1"] = merged1

            if is_paired:
                merged2 = f"{name}-graft-both.2.fastq.gz"
                concatenate_files(
                    filenames=[output_files["graft2"], output_files["both2"]],
                    out_fname=merged2,
                )
                output_files["graft2"] = merged2

        if inputs.upload_reads == "all":
            upload_list = [
                (f"{read_type}1", f"{read_type}2") if is_paired else (f"{read_type}1",)
                for read_type in READ_TYPES
            ]
        elif inputs.upload_reads != "none":
            upload_list = [
                (f"{read_type}1", f"{read_type}2") if is_paired else (f"{read_type}1",)
                for read_type in inputs.upload_reads.split(", ")
            ]
        else:
            upload_list = []

        for output_type in upload_list:
            if is_paired:
                upload_slug = "upload-fastq-paired"
                upload_inputs = {
                    "src1": [str(output_files[output_type[0]])],
                    "src2": [str(output_files[output_type[1]])],
                }
            else:
                upload_slug = "upload-fastq-single"
                upload_inputs = {"src": [str(output_files[output_type[0]])]}

            self.run_process(
                slug=upload_slug,
                inputs=upload_inputs,
            )
