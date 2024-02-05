"""GATK joint sample analysis."""

import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DateField,
    FileField,
    GroupField,
    IntegerField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class CheMut(Process):
    """CheMut varint calling using multiple BAM input files."""

    slug = "vc-chemut"
    name = "Variant calling (CheMut)"
    category = "WGS"
    process_type = "data:variants:vcf:chemut"
    version = "3.0.1"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    data_name = "Called variants (CheMut)"

    class Input:
        """Input fields for CheMut."""

        genome = DataField(data_type="seq:nucleotide", label="Reference genome")
        parental_strains = ListField(
            DataField(data_type="alignment:bam"), label="Parental strains"
        )
        mutant_strains = ListField(
            DataField(data_type="alignment:bam"), label="Mutant strains"
        )
        base_recalibration = BooleanField(
            label="Do variant base recalibration", default=False
        )
        known_sites = DataField(
            data_type="variants:vcf",
            label="dbSNP file",
            description="Database of known polymorphic sites.",
            required=False,
        )
        known_indels = ListField(
            DataField(data_type="variants:vcf"),
            label="Known indels",
            required=False,
            hidden="!base_recalibration",
        )

        class ReadsInfo:
            """Reads information."""

            PL = StringField(
                label="Platform/technology",
                description="Platform/technology used to produce the reads.",
                choices=[
                    ("Capillary", "Capillary"),
                    ("Ls454", "Ls454"),
                    ("Illumina", "Illumina"),
                    ("SOLiD", "SOLiD"),
                    ("Helicos", "Helicos"),
                    ("IonTorrent", "IonTorrent"),
                    ("Pacbio", "Pacbio"),
                ],
                default="Illumina",
            )
            LB = StringField(label="Library", default="x")
            PU = StringField(
                label="Platform unit",
                default="x",
                description="Platform unit (e.g. flowcell-barcode.lane for "
                "Illumina or slide for SOLiD). Unique identifier.",
            )
            CN = StringField(
                label="Sequencing center",
                default="x",
                description="Name of sequencing center producing the read.",
            )
            DT = DateField(
                label="Date",
                default="2017-01-01",
                description="Date the run was produced.",
            )

        class HaplotypeCaller:
            """Options for GATK HaplotypeCaller."""

            intervals = DataField(
                data_type="bed",
                label="Intervals (from BED file)",
                description="Use this option to perform the analysis over only part of the genome.",
                required=False,
            )
            ploidy = IntegerField(
                label="Sample ploidy",
                description="Ploidy (number of chromosomes) per sample. For pooled data, set "
                "to (Number of samples in each pool * Sample Ploidy).",
                default=2,
            )
            stand_call_conf = IntegerField(
                label="Min call confidence threshold",
                default=30,
                description="The minimum phred-scaled confidence threshold at which "
                "variants should be called.",
            )
            mbq = IntegerField(
                label="Min Base Quality",
                default=10,
                description="Minimum base quality required to consider a base for calling.",
            )
            max_reads = IntegerField(
                label="Max reads per alignment start site",
                default=50,
                description="Maximum number of reads to retain per alignment start position. "
                "Reads above this threshold will be downsampled. Set to 0 to disable.",
            )

        class Advanced:
            """Advanced options."""

            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of "
                "the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        reads_info = GroupField(ReadsInfo, label="Reads information")
        hc = GroupField(HaplotypeCaller, label="HaplotypeCaller options")
        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields for CheMut."""

        vcf = FileField(label="Called variants file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        if (inputs.base_recalibration and not inputs.known_sites) or (
            inputs.base_recalibration and not inputs.known_indels
        ):
            self.error(
                "Variant base recalibration requires known sites/indels "
                "information in the form of user provided VCF files."
            )

        samples = [
            Path(bam.output.bam.path).stem
            for bam in inputs.parental_strains + inputs.mutant_strains
        ]

        if len(samples) > len(set(samples)):
            self.error("Sample names must be unique.")

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )
        TMPDIR = os.environ.get("TMPDIR")
        samples_list = "samples.list"
        for counter, bam in list(
            enumerate(inputs.parental_strains + inputs.mutant_strains, start=1)
        ):
            bam_file = Path(bam.output.bam.path).stem

            args_markduplicates = [
                "--java-options",
                f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                f"INPUT={bam.output.bam.path}",
                f"OUTPUT={bam_file}_inds.bam",
                "METRICS_FILE=junk.txt",
                "VALIDATION_STRINGENCY=LENIENT",
                f"TMP_DIR={TMPDIR}",
            ]
            return_code, stdout, stderr = Cmd["gatk"]["MarkDuplicates"][
                args_markduplicates
            ] & TEE(retcode=None)
            if return_code:
                print(stdout, stderr)
                self.error("GATK MarkDuplicates tool failed.")

            args_groups = [
                "--java-options",
                f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                f"INPUT={bam_file}_inds.bam",
                f"OUTPUT={bam_file}_indh.bam",
                f"RGID=ReadGroup_{counter}",
                f"RGLB={inputs.reads_info.LB}",
                f"RGPL={inputs.reads_info.PL}",
                f"RGPU={inputs.reads_info.PU}",
                f"RGCN={inputs.reads_info.CN}",
                f"RGDT={inputs.reads_info.DT}",
                f"TMP_DIR={TMPDIR}",
                "CREATE_INDEX=TRUE",
            ]

            if bam in inputs.parental_strains:
                args_groups.append(f"RGSM=parental_{bam_file}")
            else:
                args_groups.append(f"RGSM=mut_{bam_file}")

            return_code, stdout, stderr = Cmd["gatk"]["AddOrReplaceReadGroups"][
                args_groups
            ] & TEE(retcode=None)
            if return_code:
                print(stdout, stderr)
                self.error("GATK AddOrReplaceReadGroups tool failed.")

            counter += 1

            if inputs.base_recalibration:
                args_br = [
                    "--java-options",
                    f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                    "--input",
                    f"{bam_file}_indh.bam",
                    "--reference",
                    inputs.genome.output.fasta.path,
                    "--output",
                    "recal_data.table",
                    "--tmp-dir",
                    TMPDIR,
                ]

                if inputs.known_sites:
                    args_br.extend(
                        ["--known-sites", inputs.known_sites.output.vcf.path]
                    )
                if inputs.known_indels:
                    for indel in inputs.known_indels:
                        args_br.extend(["--known-sites", indel.output.vcf.path])

                return_code, stdout, stderr = Cmd["gatk"]["BaseRecalibrator"][
                    args_br
                ] & TEE(retcode=None)
                if return_code:
                    print(stdout, stderr)
                    self.error("GATK BaseRecalibrator tool failed.")

                args_ab = [
                    "--java-options",
                    f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
                    "--input",
                    f"{bam_file}_indh.bam",
                    "--output",
                    f"{bam_file}_final.bam",
                    "--reference",
                    inputs.reference.output.fasta.path,
                    "--bqsr-recal-file",
                    "recal_data.table",
                    "--tmp-dir",
                    TMPDIR,
                ]

                return_code, stdout, stderr = Cmd["gatk"]["ApplyBQSR"][args_ab] & TEE(
                    retcode=None
                )
                if return_code:
                    print(stdout, stderr)
                    self.error("GATK ApplyBQSR tool failed.")

                with open(samples_list, "a") as f:
                    f.write(f"{bam_file}_final.bam\n")

            else:
                with open(samples_list, "a") as f:
                    f.write(f"{bam_file}_indh.bam\n")

        variants_gz = "GATKvariants_raw.vcf.gz"
        variants_index = variants_gz + ".tbi"
        args_hc = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--input",
            samples_list,
            "--reference",
            inputs.genome.output.fasta.path,
            "--output",
            variants_gz,
            "--standard-min-confidence-threshold-for-calling",
            inputs.hc.stand_call_conf,
            "--min-base-quality-score",
            inputs.hc.mbq,
            "--max-reads-per-alignment-start",
            inputs.hc.max_reads,
            "--sample-ploidy",
            inputs.hc.ploidy,
        ]

        if inputs.known_sites:
            args_hc.extend(["--known-sites", inputs.known_sites.output.vcf.path])
        if inputs.hc.intervals:
            args_hc.extend(["-L", inputs.hc.intervals.output.bed.path])

        return_code, stdout, stderr = Cmd["gatk"]["HaplotypeCaller"][args_hc] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("GATK HaplotypeCaller tool failed.")

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = inputs.parental_strains[0].output.species
        outputs.build = inputs.parental_strains[0].output.build
