import os

from resolwe.flow.models import Process
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class SlamdunkProcessorTestCase(KBBioProcessTestCase):
    @tag_process("slamdunk-all-paired")
    def test_slamdunk_paired(self):
        with self.preparation_stage():
            paired_reads = self.prepare_paired_reads(
                ["hs_slamseq_R1_complemented.fastq.gz"], ["hs_slamseq_R2.fastq.gz"]
            )
            transcripts = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.fasta"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            bedfile = self.run_process(
                "upload-bed",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.bed"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
        inputs = {
            "reads": paired_reads.id,
            "ref_seq": transcripts.id,
            "regions": bedfile.id,
            "filter_multimappers": True,
            "max_alignments": 1,
            "read_length": 75,
        }
        slamdunk = self.run_process("slamdunk-all-paired", inputs)
        self.assertFile(
            slamdunk,
            "tcount",
            os.path.join("slamseq", "output", "hs_slamseq_tcount.tsv"),
        )

    @tag_process("alleyoop-rates")
    def test_alleyoop_rates(self):
        with self.preparation_stage():
            reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.fasta"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            process = Process.objects.create(
                name="Upload Slamdunk data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alignment:bam:slamdunk:",
                input_schema=[
                    {"name": "src", "type": "basic:file:",},
                    {"name": "index", "type": "basic:file:",},
                ],
                output_schema=[
                    {"name": "bam", "type": "basic:file:",},
                    {"name": "bai", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai
re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )
            slamdunk = self.run_process(
                process.slug,
                {
                    "src": os.path.join(
                        "slamseq", "output", "hs_slamseq_slamdunk_mapped_filtered.bam"
                    ),
                    "index": os.path.join(
                        "slamseq",
                        "output",
                        "hs_slamseq_slamdunk_mapped_filtered.bam.bai",
                    ),
                },
            )

        alleyoop_rates = self.run_process(
            "alleyoop-rates", {"ref_seq": reference.id, "slamdunk": slamdunk.id}
        )
        self.assertFile(
            alleyoop_rates,
            "report",
            os.path.join("slamseq", "output", "hs_alleyoop_overallrates.txt"),
        )

    @tag_process("alleyoop-utr-rates")
    def test_alleyoop_utrrates(self):
        with self.preparation_stage():
            reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.fasta"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            bedfile = self.run_process(
                "upload-bed",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.bed"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            process = Process.objects.create(
                name="Upload Slamdunk data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alignment:bam:slamdunk:",
                input_schema=[
                    {"name": "src", "type": "basic:file:",},
                    {"name": "index", "type": "basic:file:",},
                ],
                output_schema=[
                    {"name": "bam", "type": "basic:file:",},
                    {"name": "bai", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai
re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )
            slamdunk = self.run_process(
                process.slug,
                {
                    "src": os.path.join(
                        "slamseq", "output", "hs_slamseq_slamdunk_mapped_filtered.bam"
                    ),
                    "index": os.path.join(
                        "slamseq",
                        "output",
                        "hs_slamseq_slamdunk_mapped_filtered.bam.bai",
                    ),
                },
            )
        rates = self.run_process(
            "alleyoop-utr-rates",
            {"ref_seq": reference.id, "regions": bedfile.id, "slamdunk": slamdunk.id},
        )
        self.assertFile(
            rates,
            "report",
            os.path.join("slamseq", "output", "hs_alleyoop_mutationrates.txt"),
        )

    @tag_process("alleyoop-summary")
    def test_alleyoop_summary(self):
        with self.preparation_stage():
            process = Process.objects.create(
                name="Upload Slamdunk data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alignment:bam:slamdunk:",
                input_schema=[
                    {"name": "src", "type": "basic:file:",},
                    {"name": "index", "type": "basic:file:",},
                    {"name": "tsv", "type": "basic:file:",},
                ],
                output_schema=[
                    {"name": "bam", "type": "basic:file:",},
                    {"name": "bai", "type": "basic:file:",},
                    {"name": "tcount", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai

re-import {{ tsv.file_temp|default(tsv.file) }} {{ tsv.file }} "tsv" "tsv" 0.1 extract
re-save-file tcount "${NAME}".tsv

re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )
            slamdunk = self.run_process(
                process.slug,
                {
                    "src": os.path.join(
                        "slamseq", "output", "hs_slamseq_slamdunk_mapped_filtered.bam"
                    ),
                    "index": os.path.join(
                        "slamseq",
                        "output",
                        "hs_slamseq_slamdunk_mapped_filtered.bam.bai",
                    ),
                    "tsv": os.path.join("slamseq", "output", "hs_slamseq_tcount.tsv"),
                },
            )
        summary = self.run_process("alleyoop-summary", {"slamdunk": [slamdunk.id]})
        self.assertFile(
            summary,
            "report",
            os.path.join("slamseq", "output", "hs_alleyoop_summary.txt"),
        )

    @tag_process("alleyoop-snpeval")
    def test_alleyoop_snpeval(self):
        with self.preparation_stage():
            reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.fasta"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            bedfile = self.run_process(
                "upload-bed",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.bed"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            process = Process.objects.create(
                name="Upload Slamdunk data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alignment:bam:slamdunk:",
                input_schema=[
                    {"name": "src", "type": "basic:file:",},
                    {"name": "index", "type": "basic:file:",},
                    {"name": "snp", "type": "basic:file:",},
                ],
                output_schema=[
                    {"name": "bam", "type": "basic:file:",},
                    {"name": "bai", "type": "basic:file:",},
                    {"name": "variants", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai

re-import {{ snp.file_temp|default(snp.file) }} {{ snp.file }} "vcf" "vcf" 0.1 extract
re-save-file variants "${NAME}".vcf

re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )
            slamdunk = self.run_process(
                process.slug,
                {
                    "src": os.path.join(
                        "slamseq", "output", "hs_slamseq_slamdunk_mapped_filtered.bam"
                    ),
                    "index": os.path.join(
                        "slamseq",
                        "output",
                        "hs_slamseq_slamdunk_mapped_filtered.bam.bai",
                    ),
                    "snp": os.path.join("slamseq", "output", "hs_slamseq_snp.vcf"),
                },
            )
        summary = self.run_process(
            "alleyoop-snpeval",
            {"ref_seq": reference.id, "regions": bedfile.id, "slamdunk": slamdunk.id},
        )
        self.assertFile(
            summary,
            "report",
            os.path.join("slamseq", "output", "hs_alleyoop_SNPeval.txt"),
        )

    @with_resolwe_host
    @tag_process("alleyoop-collapse")
    def test_alleyoop_collapse(self):
        with self.preparation_stage():
            process = Process.objects.create(
                name="Upload Slamdunk data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alignment:bam:slamdunk:",
                input_schema=[{"name": "src", "type": "basic:file:",},],
                output_schema=[
                    {"name": "tcount", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "tsv" "tsv" 0.1 extract
re-save-file tcount "${NAME}".tsv
re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )

            slamdunk = self.run_process(
                process.slug,
                {"src": os.path.join("slamseq", "output", "hs_slamseq_tcount.tsv")},
            )

        alleyoop_collapse = self.run_process(
            "alleyoop-collapse", {"slamdunk": slamdunk.id,}
        )
        self.assertFile(
            alleyoop_collapse,
            "tcount",
            os.path.join("slamseq", "output", "collapsed_tcount.txt"),
        )

    @with_resolwe_host
    @tag_process("slam-count")
    def test_slam_count(self):
        with self.preparation_stage():
            process = Process.objects.create(
                name="Upload Slamdunk collapse data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                contributor=self.contributor,
                type="data:alleyoop:collapse:",
                input_schema=[{"name": "src", "type": "basic:file:",},],
                output_schema=[
                    {"name": "tcount", "type": "basic:file:",},
                    {"name": "species", "type": "basic:string:",},
                    {"name": "build", "type": "basic:string:",},
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 extract
re-save-file tcount "${NAME}".txt
re-save species "Homo sapiens"
re-save build "Gencode 32"
""",
                },
            )

            collapsed_input = self.run_process(
                process.slug,
                {"src": os.path.join("slamseq", "output", "collapsed_tcount.txt")},
            )

        slam_count = self.run_process("slam-count", {"tcount": collapsed_input.id,})
        self.assertFile(
            slam_count,
            "exp",
            os.path.join("slamseq", "output", "tcount_exp_tpm.txt.gz"),
            compression="gzip",
        )
        self.assertFile(
            slam_count,
            "exp_set",
            os.path.join("slamseq", "output", "tcount_exp_set.txt.gz"),
            compression="gzip",
        )
        self.assertFields(slam_count, "exp_type", "TPM")
        self.assertFields(slam_count, "source", "ENSEMBL")
        self.assertFields(slam_count, "species", "Homo sapiens")
        self.assertFields(slam_count, "build", "Gencode 32")
        self.assertFields(slam_count, "feature_type", "gene")
