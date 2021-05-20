from pathlib import Path

from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WgsProcessorTestCase(BioProcessTestCase):
    @tag_process("wgs-preprocess")
    def test_wgs_preprocess(self):
        def filter_startedon(line):
            """Filter stared on header line."""
            return line.startswith(b"# Started on:") or line.startswith(
                b"# MarkDuplicates"
            )

        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": inputs / "hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            reads = self.prepare_paired_reads(
                mate1=[inputs / "TP53_1.fastq.gz"],
                mate2=[inputs / "TP53_2.fastq.gz"],
            )

            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "dbsnp_TP53.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        analysis_ready_bam = self.run_process(
            "wgs-preprocess",
            {
                "reads": reads.id,
                "ref_seq": ref_seq.id,
                "bwa_index": bwa_index.id,
                "known_sites": [dbsnp.id],
            },
        )

        self.assertFile(
            analysis_ready_bam,
            "stats",
            outputs / "wgs_preprocess_bam_stats.txt",
        )

        self.assertFile(
            analysis_ready_bam,
            "metrics_file",
            outputs / "wgs_preprocess_markdups_metrics.txt",
            file_filter=filter_startedon,
        )
