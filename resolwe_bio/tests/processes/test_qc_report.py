# import os
# import pandas as pd
# import plotly.graph_objects as go
# from pathlib import Path

from resolwe.flow.models import Data, Process
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import KBBioProcessTestCase


class QcReportTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("qc-repor")
    def test_qc_report_process_fc(self):
        with self.preparation_stage():
            inputs = {
                "src": [
                    "qc_report/input/hs_single bbduk_star_htseq_reads_single.fastq.gz",
                    "qc_report/input/hs sim_reads_single.fastq.gz",
                    "qc_report/input/merged_single_end_reads.fastq.gz",
                ]
            }
            reads = self.run_processor("upload-fastq-single", inputs)

            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/hs genome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            adapters = self.prepare_ref_seq()

            inputs = {
                "src": "qc_report/input/hs annotation.gtf.gz",
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "build": "ens_90",
            }
            annotation = self.run_process("upload-gtf", inputs)

            inputs = {"annotation": annotation.id, "ref_seq": star_index_fasta.id}
            star_index = self.run_process("alignment-star-index", inputs)

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genomeSAindexNbases": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genomeSAindexNbases": 2,
                    },
                },
            )

            fc = self.run_process(
                "workflow-bbduk-star-fc-quant-single",
                {
                    "reads": reads.id,
                    "star_index": star_index.id,
                    "adapters": [adapters.id],
                    "annotation": annotation.id,
                    "stranded": "forward",
                    "qc": {
                        "rrna_reference": rrna_star_index.id,
                        "globin_reference": globin_star_index.id,
                    },
                },
            )

            multiqc = Data.objects.filter(process__slug="multiqc").last()

        inputs = {"multiqc": multiqc}
        qc_report = self.run_process("qc_report", inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFileExists(qc_report, "general_stats_table.tsv")
        self.assertFileExists(qc_report, "qc_report.html")

    @with_resolwe_host
    @tag_process("qc-repor")
    def test_qc_report_process_salmon(self):
        with self.preparation_stage():
            inputs = {
                "src": [
                    "qc_report/input/hs_single bbduk_star_htseq_reads_single.fastq.gz",
                    "qc_report/input/hs sim_reads_single.fastq.gz",
                    "qc_report/input/merged_single_end_reads.fastq.gz",
                ]
            }
            reads = self.run_processor("upload-fastq-single", inputs)

            annotation = self.prepare_annotation(
                fn="qc_report/input/hs annotation.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens92",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/hs genome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            star_index = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation.id,
                    "ref_seq": star_index_fasta.id,
                },
            )
            cdna = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/hs cdna.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": cdna.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            adapters = self.prepare_ref_seq()
            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genomeSAindexNbases": 2,
                    },
                },
            )
            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qc_report/input/Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genomeSAindexNbases": 2,
                    },
                },
            )

            inputs = {
                "reads": reads.id,
                "genome": star_index.id,
                "salmon_index": salmon_index.id,
                "annotation": annotation.id,
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
                "preprocessing": {
                    "adapters": [adapters.id],
                },
                "quantification": {
                    "min_assigned_frag": 1,
                },
            }

            self.run_process("workflow-bbduk-salmon-qc-single", inputs)

            multiqc = Data.objects.filter(process__slug="multiqc").last()

        inputs = {"multiqc": multiqc}
        qc_report = self.run_process("qc_report", inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFileExists(qc_report, "general_stats_table.tsv")
        self.assertFileExists(qc_report, "qc_report.html")
