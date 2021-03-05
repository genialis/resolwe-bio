import os

from resolwe.flow.models import Process
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WgbsProcessorTestCase(BioProcessTestCase):
    @tag_process("walt-index")
    def test_walt_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn="./wgbs/input/g.en ome.fasta.gz",
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )

        walt_index = self.run_process("walt-index", {"ref_seq": ref_seq.id})
        self.assertDir(
            walt_index, "index", os.path.join("wgbs", "output", "walt_index.tar.gz")
        )
        self.assertFile(
            walt_index, "fasta", os.path.join("wgbs", "output", "genome.fasta")
        )
        self.assertFile(
            walt_index,
            "fastagz",
            os.path.join("wgbs", "output", "genome.fasta.gz"),
            compression="gzip",
        )
        self.assertFile(
            walt_index, "fai", os.path.join("wgbs", "output", "genome.fasta.fai")
        )
        self.assertFields(walt_index, "species", "Dictyostelium discoideum")
        self.assertFields(walt_index, "build", "dd-05-2009")

    @tag_process("walt")
    def test_walt(self):
        with self.preparation_stage():
            inputs = {
                "src": os.path.join("wgbs", "input", "hg19_chr2 17k.fasta.gz"),
                "species": "Homo sapiens",
                "build": "hg19",
            }
            ref_seq = self.run_process("upload-fasta-nucl", inputs)
            walt_index = self.run_process("walt-index", {"ref_seq": ref_seq.id})
            reads_paired = self.prepare_paired_reads(
                mate1=[
                    os.path.join("wgbs", "input", "3A_WT_WGBS_chr2_17k_R1.fastq.gz")
                ],
                mate2=[
                    os.path.join("wgbs", "input", "3A_WT_WGBS_chr2_17k_R2.fastq.gz")
                ],
            )

        inputs = {
            "genome": walt_index.id,
            "reads": reads_paired.id,
            "spikein_options": {"spikein_name": "chr2"},
        }
        walt = self.run_process("walt", inputs)
        self.assertFile(
            walt, "stats", os.path.join("wgbs", "output", "walt_report.txt")
        )
        self.assertFileExists(walt, "bam")
        self.assertFileExists(walt, "bai")
        self.assertFileExists(walt, "mr")
        self.assertFileExists(walt, "unmapped")
        self.assertFileExists(walt, "spikein_mr")
        self.assertFields(walt, "species", "Homo sapiens")
        self.assertFields(walt, "build", "hg19")

    @tag_process("methcounts")
    def test_methcounts(self):
        with self.preparation_stage():
            inputs = {
                "src": os.path.join("wgbs", "input", "hg19_chr2 17k.fasta.gz"),
                "species": "Homo sapiens",
                "build": "hg19",
            }
            genome = self.run_process("upload-fasta-nucl", inputs)

            # Mock upload WALT alignment process
            process = Process.objects.create(
                name="Upload WALT alignment file (.mr) mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/s4q6j6e8/resolwebio/base:ubuntu-20.04-03042021",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:alignment:bam:walt:",
                input_schema=[
                    {
                        "name": "fc",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "mr",
                        "type": "basic:file:",
                    },
                    {
                        "name": "species",
                        "type": "basic:string:",
                    },
                    {
                        "name": "build",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "mr" "mr" 0.1 compress
re-save-file mr "${NAME}".mr.gz
re-save species 'Homo sapiens'
re-save build 'hg19'
""",
                },
            )

            inputs = {"fc": os.path.join("wgbs", "output", "3A_WT_WGBS_chr2_17k.mr.gz")}
            walt = self.run_process(process.slug, inputs)

        inputs = {
            "genome": genome.id,
            "alignment": walt.id,
        }
        methcounts = self.run_process("methcounts", inputs)
        self.assertFile(
            methcounts, "stats", os.path.join("wgbs", "output", "methconts_report.txt")
        )
        self.assertFileExists(methcounts, "meth")
        self.assertFileExists(methcounts, "bigwig")
        self.assertFields(methcounts, "species", "Homo sapiens")
        self.assertFields(methcounts, "build", "hg19")

    @tag_process("hmr")
    def test_hmr(self):
        with self.preparation_stage():
            # Mock upload methcounts process
            process = Process.objects.create(
                name="Upload methcounts file (.meth) mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/s4q6j6e8/resolwebio/base:ubuntu-20.04-03042021",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:wgbs:methcounts:",
                input_schema=[
                    {
                        "name": "fc",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "meth",
                        "type": "basic:file:",
                    },
                    {
                        "name": "species",
                        "type": "basic:string:",
                    },
                    {
                        "name": "build",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "meth" "meth" 0.1 compress
re-save-file meth "${NAME}".meth.gz
re-save species 'Homo sapiens'
re-save build 'hg19'
""",
                },
            )

            inputs = {
                "fc": os.path.join("wgbs", "output", "3A_WT_WGBS_chr2_17k.meth.gz")
            }
            methcounts = self.run_process(process.slug, inputs)

        hmr = self.run_process("hmr", {"methcounts": methcounts.id})
        self.assertFile(
            hmr,
            "hmr",
            os.path.join("wgbs", "output", "3A_WT_WGBS_chr2_17k.hmr.gz"),
            compression="gzip",
        )
        self.assertFields(hmr, "species", "Homo sapiens")
        self.assertFields(hmr, "build", "hg19")

    @tag_process("bs-conversion-rate")
    def test_bsrate(self):
        with self.preparation_stage():
            sequence = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join(
                        "wgbs", "input", "unmethylated_lambda_J02459.fasta.gz"
                    ),
                    "species": "Homo sapiens",
                    "build": "J02459",
                },
            )
            # Mock upload WALT alignment process
            process = Process.objects.create(
                name="Upload WALT alignment file (.mr) mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/s4q6j6e8/resolwebio/base:ubuntu-20.04-03042021",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:alignment:bam:walt:",
                input_schema=[
                    {
                        "name": "fc",
                        "type": "basic:file:",
                    },
                    {
                        "name": "f",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "mr",
                        "type": "basic:file:",
                    },
                    {
                        "name": "spikein_mr",
                        "type": "basic:file:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "mr" "mr" 0.1 compress
re-save-file mr "${NAME}".mr.gz
re-import {{ f.file_temp|default(f.file) }} {{ f.file }} "mr" "mr" 0.1 compress
re-save-file spikein_mr "${NAME}".mr.gz
""",
                },
            )

            inputs = {
                "fc": os.path.join("wgbs", "input", "spike_in_lambda.mr.gz"),
                "f": os.path.join("wgbs", "input", "spike_in_lambda.mr.gz"),
            }
            walt = self.run_process(process.slug, inputs)

        bsrate = self.run_process(
            "bs-conversion-rate", {"mr": walt.id, "sequence": sequence.id}
        )
        self.assertFile(
            bsrate,
            "report",
            os.path.join("wgbs", "output", "Escherichia_phage_bsrate.txt"),
        )
