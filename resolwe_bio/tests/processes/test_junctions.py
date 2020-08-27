from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class JunctionsProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("regtools-junctions-annotate")
    def test_regtools_annotate(self):
        with self.preparation_stage():
            reads = self.prepare_reads(["SRR2141558_chr1_20k.fq.gz"])
            annotation = self.prepare_annotation(
                fn="Homo_sapiens.GRCh38.92.chr1_20k.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="GRCh38_ens92",
            )

            genome = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens.GRCh38.92.chr1_20k.fa.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens92",
                },
            )

            star_index = self.run_process(
                "alignment-star-index",
                {"annotation": annotation.id, "ref_seq": genome.id},
            )

            star_sj = self.run_process(
                "alignment-star",
                {
                    "genome": star_index.id,
                    "reads": reads.id,
                },
            )

            bed_upload = self.run_process(
                "upload-bed",
                {
                    "src": "SRR2141558_chr1_20k_SJ.bed",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens92",
                },
            )

        inputs = {
            "alignment_star": star_sj.id,
            "genome": genome.id,
            "annotation": annotation.id,
        }
        junctions_star = self.run_process("regtools-junctions-annotate", inputs)
        self.assertFile(junctions_star, "bed", "regtools.bed.gz", compression="gzip")
        self.assertFile(
            junctions_star, "splice_junctions", "regtools_SJ.txt.gz", compression="gzip"
        )
        self.assertFile(
            junctions_star, "novel_splice_junctions", "regtools_novel_SJ.txt"
        )

        inputs = {
            "alignment": star_sj.id,
            "genome": genome.id,
            "annotation": annotation.id,
        }
        junctions_bam = self.run_process("regtools-junctions-annotate", inputs)
        self.assertFile(
            junctions_bam, "bed", "regtools_from_bam.bed.gz", compression="gzip"
        )
        self.assertFile(
            junctions_bam,
            "splice_junctions",
            "regtools_SJ_from_bam.txt.gz",
            compression="gzip",
        )
        self.assertFile(
            junctions_bam, "novel_splice_junctions", "regtools_novel_SJ_from_bam.txt"
        )

        inputs = {
            "input_bed_junctions": bed_upload.id,
            "genome": genome.id,
            "annotation": annotation.id,
        }
        junctions_bed = self.run_process("regtools-junctions-annotate", inputs)
        self.assertFile(
            junctions_bed, "splice_junctions", "regtools_SJ.txt.gz", compression="gzip"
        )
        self.assertFile(
            junctions_bed, "novel_splice_junctions", "regtools_novel_SJ.txt"
        )

    @with_resolwe_host
    @tag_process("regtools-junctions-annotate")
    def test_no_sj_inputs(self):
        with self.preparation_stage():
            reads_no_junctions = self.prepare_reads(["SRR2124780_1 1k.fastq.gz"])
            annotation_no_junctions = self.prepare_annotation(
                fn="HS chr21_short.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="GRCh38_ens90",
            )

            genome_no_junctions = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "HS chr21_ensembl.fa.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

            star_index_no_junctions = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation_no_junctions.id,
                    "ref_seq": genome_no_junctions.id,
                },
            )

            star_no_junctions = self.run_process(
                "alignment-star",
                {
                    "genome": star_index_no_junctions.id,
                    "reads": reads_no_junctions.id,
                },
            )

            empty_bam_upload = self.run_process(
                "upload-bam",
                {
                    "src": "empty.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

        inputs = {
            "alignment_star": star_no_junctions.id,
            "genome": genome_no_junctions.id,
            "annotation": annotation_no_junctions.id,
        }
        annotate_star_no_junctions = self.run_process(
            "regtools-junctions-annotate", inputs
        )
        warning_msg = [
            "STAR SJ.out.tab file has no entries. There will be no splice junctions detected.",
            "Bed file has no entries.",
            "BigBed index can not be created.",
            "Bed file with novel splice junctions has no entries. BigBed index can not be created.",
        ]
        self.assertEqual(annotate_star_no_junctions.process_warning, warning_msg)

        inputs = {
            "alignment": star_no_junctions.id,
            "genome": genome_no_junctions.id,
            "annotation": annotation_no_junctions.id,
        }
        annotate_bam_no_junctions = self.run_process(
            "regtools-junctions-annotate", inputs
        )
        warning_msg = [
            "Bed file has no entries.",
            "BigBed index can not be created.",
            "Bed file with novel splice junctions has no entries. BigBed index can not be created.",
        ]
        self.assertEqual(annotate_bam_no_junctions.process_warning, warning_msg)

        inputs["alignment"] = empty_bam_upload.id
        annotate_empty_bam = self.run_process("regtools-junctions-annotate", inputs)
        warning_msg = [
            "Bam file has no entries. There will be no splice junctions detected.",
            "Bed file has no entries.",
            "BigBed index can not be created.",
            "Bed file with novel splice junctions has no entries. BigBed index can not be created.",
        ]
        self.assertEqual(annotate_empty_bam.process_warning, warning_msg)
