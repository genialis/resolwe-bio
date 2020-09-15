from os.path import join
from pathlib import Path

from resolwe.flow.models import Collection, Data, Process, Relation
from resolwe.flow.models.entity import RelationPartition, RelationType
from resolwe.test import tag_process

from resolwe_bio.expression_filters.relation import background_pairs
from resolwe_bio.models import Sample
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class ChipSeqProcessorTestCase(BioProcessTestCase):

    fixtures = ["relationtypes.yaml"]

    @tag_process("chipseq-peakscore", "chipseq-genescore")
    def test_chipseq(self):
        with self.preparation_stage():
            inputs = {
                "src": "chip_seq_control.bam",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            control_bam = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "chip_seq_case.bam",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            case_bam = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "chip_seq.bed",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            bed = self.run_process("upload-bed", inputs)

            inputs = {
                "case": case_bam.pk,
                "control": control_bam.pk,
                "settings": {
                    "nomodel": True,
                    "pvalue": 0.001,
                    "slocal": 2000,
                    "extsize": 100,
                    "call_summits": True,
                },
            }
            macs2 = self.run_process("macs2-callpeak", inputs)

        inputs = {"peaks": macs2.pk, "bed": bed.pk}
        peak_score = self.run_process("chipseq-peakscore", inputs)
        self.assertFile(peak_score, "peak_score", "chip_seq_peakscore_genomicContext")

        inputs = {"peakscore": peak_score.id}
        gene_score = self.run_process("chipseq-genescore", inputs)
        self.assertFile(gene_score, "genescore", "chip_seq_geneScore.xls")

    @tag_process("macs14")
    def test_macs14(self):
        with self.preparation_stage():
            inputs = {
                "src": "macs14_control.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            control_bam = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs14_case.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_bam = self.run_process("upload-bam", inputs)

        inputs = {"treatment": case_bam.id, "control": control_bam.id}
        macs14 = self.run_process("macs14", inputs)

        self.assertFields(macs14, "species", "Homo sapiens")
        self.assertFields(macs14, "build", "hg19")
        self.assertFile(macs14, "peaks_bed", "macs14_peaks.bed.gz")
        self.assertFile(macs14, "peaks_bigbed_igv_ucsc", "macs14_peaks.bb")
        self.assertFile(macs14, "peaks_tbi_jbrowse", "macs14_peaks.gz.tbi")
        self.assertFile(macs14, "summits_tbi_jbrowse", "macs14_summits.gz.tbi")
        self.assertFile(macs14, "treat_bigwig", "macs14_treat.bw")

    @tag_process("macs2-callpeak", "archive-samples")
    def test_macs2(self):
        with self.preparation_stage():
            inputs = {
                "src": "macs2/input/SRR5675973_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_1 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675974_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            control_1 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675975_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_2 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675976_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            control_2 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/paired_end_discodeum.bam",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            paired_end = self.run_process("upload-bam", inputs)

            promoters = self.run_process(
                "upload-bed",
                {
                    "src": "macs2/input/promoter_regions.bed",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

        inputs = {
            "case": case_1.id,
            "control": control_1.id,
            "promoter": promoters.id,
            "settings": {
                "extsize": 298,
                "nomodel": True,
                "bedgraph": True,
            },
        }
        macs_sample1 = self.run_process("macs2-callpeak", inputs)

        self.assertFields(macs_sample1, "species", "Homo sapiens")
        self.assertFields(macs_sample1, "build", "hg19")
        self.assertFile(
            macs_sample1, "case_prepeak_qc", "macs2/output/case_prepeak_qc.txt"
        )
        self.assertFile(
            macs_sample1, "control_prepeak_qc", "macs2/output/control_prepeak_qc.txt"
        )
        self.assertFile(macs_sample1, "chip_qc", "macs2/output/postpeak_qc_report.txt")
        self.assertFileExists(macs_sample1, "called_peaks")
        self.assertFileExists(macs_sample1, "narrow_peaks")
        self.assertFileExists(macs_sample1, "narrow_peaks_bigbed_igv_ucsc")
        self.assertFileExists(macs_sample1, "summits")
        self.assertFileExists(macs_sample1, "summits_tbi_jbrowse")
        self.assertFileExists(macs_sample1, "summits_bigbed_igv_ucsc")
        self.assertFileExists(macs_sample1, "treat_pileup")
        self.assertFileExists(macs_sample1, "treat_pileup_bigwig")
        self.assertFileExists(macs_sample1, "control_lambda")
        self.assertFileExists(macs_sample1, "control_lambda_bigwig")
        self.assertFileExists(macs_sample1, "case_bam")
        self.assertFileExists(macs_sample1, "case_bai")
        self.assertFileExists(macs_sample1, "control_bam")
        self.assertFileExists(macs_sample1, "control_bai")

        # Use filtered tagAlign file for MACS2
        inputs = {
            "case": case_2.id,
            "control": control_2.id,
            "tagalign": True,
            "settings": {
                "bedgraph": True,
            },
        }
        macs_sample2 = self.run_process("macs2-callpeak", inputs)

        self.assertFields(macs_sample2, "species", "Homo sapiens")
        self.assertFields(macs_sample2, "build", "hg19")
        self.assertFileExists(macs_sample2, "chip_qc")
        self.assertFileExists(macs_sample2, "called_peaks")
        self.assertFileExists(macs_sample2, "narrow_peaks")
        self.assertFileExists(macs_sample2, "narrow_peaks_bigbed_igv_ucsc")
        self.assertFileExists(macs_sample2, "summits")
        self.assertFileExists(macs_sample2, "summits_tbi_jbrowse")
        self.assertFileExists(macs_sample2, "summits_bigbed_igv_ucsc")

        # Test "archive-samples"' QC report merge
        inputs = {
            "data": [macs_sample1.id, macs_sample2.id],
            "fields": ["chip_qc", "case_prepeak_qc"],
        }
        self.run_process("archive-samples", inputs)

        # Run "archive-samples" without QC reports
        inputs["fields"] = ["called_peaks"]
        self.run_process("archive-samples", inputs)

        # Test BAMPE reads
        inputs = {
            "case": paired_end.id,
            "settings": {
                "format": "BAMPE",
                "pvalue": 0.05,
                "duplicates": "all",
                "bedgraph": True,
            },
        }
        macs_paired = self.run_process("macs2-callpeak", inputs)
        self.assertFields(macs_paired, "species", "Dictyostelium discoideum")
        self.assertFields(macs_paired, "build", "dd-05-2009")
        self.assertFileExists(macs_paired, "chip_qc")
        self.assertFileExists(macs_paired, "called_peaks")
        self.assertFileExists(macs_paired, "narrow_peaks")
        self.assertFileExists(macs_paired, "narrow_peaks_bigbed_igv_ucsc")
        self.assertFileExists(macs_paired, "summits")
        self.assertFileExists(macs_paired, "summits_tbi_jbrowse")
        self.assertFileExists(macs_paired, "summits_bigbed_igv_ucsc")

        inputs["settings"]["shift"] = 10
        macs_paired_err = self.run_process("macs2-callpeak", inputs, Data.STATUS_ERROR)
        err_msg = [
            "Shift values other than 0 are not supported when the format is BAMPE."
        ]
        self.assertEqual(macs_paired_err.process_error, err_msg)

        inputs = {
            "case": case_1.id,
            "settings": {
                "format": "BAMPE",
            },
        }
        macs_paired_err = self.run_process("macs2-callpeak", inputs, Data.STATUS_ERROR)
        err_msg = ["No paired-end reads were detected but BAMPE format was selected."]
        self.assertEqual(macs_paired_err.process_error, err_msg)

    @skipUnlessLargeFiles("rose2_case.bam", "rose2_control.bam")
    @tag_process("rose2")
    def test_rose2(self):
        with self.preparation_stage():
            inputs = {
                "src": join("large", "rose2_case.bam"),
                "species": "Homo sapiens",
                "build": "hg19",
            }
            bam = self.run_process("upload-bam", inputs)

            inputs = {
                "src": join("large", "rose2_control.bam"),
                "species": "Homo sapiens",
                "build": "hg19",
            }
            control = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs14_chr22.bed",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            macs_peaks = self.run_process("upload-bed", inputs)

            inputs = {
                "src": "hg19_encode_blacklist_chr22.bed",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            mask = self.run_process("upload-bed", inputs)

        inputs = {
            "input_upload": macs_peaks.id,
            "rankby": bam.id,
            "control": control.id,
            "stitch": 5000,
            "tss": 2500,
            "mask": mask.id,
        }
        rose2 = self.run_process("rose2", inputs)

        # remove changing lines from the rose2 output
        def filter_created(line):
            return line.startswith(b"#Created")

        self.assertFile(
            rose2,
            "all_enhancers",
            "rose2_enhancer_table.txt",
            file_filter=filter_created,
        )

    @tag_process("qc-prepeak")
    def test_qc_prepeak(self):
        with self.preparation_stage():
            inputs = {
                "src": "prepeak_se.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            bam = self.run_process("upload-bam", inputs)

        inputs = {
            "alignment": bam.id,
            "n_sub": 7000,
            "q_treshold": 25,
        }
        prepeak = self.run_process("qc-prepeak", inputs)

        self.assertFields(prepeak, "species", "Homo sapiens")
        self.assertFields(prepeak, "build", "hg19")
        self.assertFields(prepeak, "fraglen", 215)
        self.assertFile(prepeak, "chip_qc", "prepeak_se_qc_report.txt")
        self.assertFileExists(prepeak, "tagalign")

        with self.preparation_stage():
            inputs = {
                "src": "prepeak_pe.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            bam = self.run_process("upload-bam", inputs)

        inputs = {
            "alignment": bam.id,
            "n_sub": 7000,
            "q_treshold": 25,
        }
        prepeak = self.run_process("qc-prepeak", inputs)

        self.assertFields(prepeak, "species", "Homo sapiens")
        self.assertFields(prepeak, "build", "hg19")
        self.assertFields(prepeak, "fraglen", 235)
        self.assertFile(prepeak, "chip_qc", "prepeak_pe_qc_report.txt")
        self.assertFileExists(prepeak, "tagalign")

        # Test Tn5 shifting (ATAC-seq) and peak shifting
        inputs["tn5"] = True
        inputs["shift"] = 0
        prepeak = self.run_process("qc-prepeak", inputs)

        self.assertFields(prepeak, "species", "Homo sapiens")
        self.assertFields(prepeak, "build", "hg19")
        self.assertFields(prepeak, "fraglen", 0)
        self.assertFile(prepeak, "chip_qc", "prepeak_pe_qc_report_tn5.txt")
        self.assertFile(
            prepeak, "tagalign", "prepeak_pe_tn5.tagAlign.gz", compression="gzip"
        )

    @tag_process("macs2-batch")
    def test_macs2_batch(self):
        with self.preparation_stage():
            inputs = {
                "src": "macs2/input/SRR5675973_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_1 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675974_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            background_1 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675975_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_2 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675976_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            background_2 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675973_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_wo_background = self.run_process("upload-bam", inputs)

            promoters = self.run_process(
                "upload-bed",
                {
                    "src": "macs2/input/promoter_regions.bed",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

            collection = Collection.objects.create(
                name="Test collection", contributor=self.contributor
            )

            rel_type_background = RelationType.objects.get(name="background")

            background1 = Relation.objects.create(
                contributor=self.contributor,
                collection=collection,
                type=rel_type_background,
                category="Background",
            )

            background2 = Relation.objects.create(
                contributor=self.contributor,
                collection=collection,
                type=rel_type_background,
                category="Background2",
            )

            RelationPartition.objects.create(
                relation=background1, entity=case_1.entity, label="case"
            )
            RelationPartition.objects.create(
                relation=background1, entity=background_1.entity, label="background"
            )
            RelationPartition.objects.create(
                relation=background2, entity=case_2.entity, label="case"
            )
            RelationPartition.objects.create(
                relation=background2, entity=background_2.entity, label="background"
            )

            self.assertEqual(
                background_pairs(
                    [
                        {"__id": case_1.id, "__type": case_1.process.type},
                        {"__id": background_1.id, "__type": background_1.process.type},
                        {"__id": case_2.id, "__type": case_2.process.type},
                        {"__id": background_2.id, "__type": background_2.process.type},
                        {
                            "__id": case_wo_background.id,
                            "__type": case_wo_background.process.type,
                        },
                    ]
                ),
                [
                    (case_1.id, background_1.id),
                    (case_2.id, background_2.id),
                    (case_wo_background.id, None),
                ],
            )

        self.run_process(
            "macs2-batch",
            {
                "alignments": [
                    case_1.id,
                    background_1.id,
                    case_2.id,
                    background_2.id,
                    case_wo_background.id,
                ],
                "promoter": promoters.id,
                "tagalign": True,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        macs2 = Data.objects.filter(process__slug="macs2-callpeak").last()

        self.assertFields(macs2, "species", "Homo sapiens")
        self.assertFields(macs2, "build", "hg19")
        self.assertFile(
            macs2, "case_prepeak_qc", "macs2/output/batch_case_prepeak_qc.txt"
        )
        self.assertFile(
            macs2, "chip_qc", "macs2/output/batch_case_postpeak_qc_report.txt"
        )
        self.assertFileExists(macs2, "called_peaks")
        self.assertFileExists(macs2, "narrow_peaks")
        self.assertFileExists(macs2, "narrow_peaks_bigbed_igv_ucsc")
        self.assertFileExists(macs2, "summits")
        self.assertFileExists(macs2, "summits_tbi_jbrowse")
        self.assertFileExists(macs2, "summits_bigbed_igv_ucsc")
        self.assertFileExists(macs2, "treat_pileup")
        self.assertFileExists(macs2, "treat_pileup_bigwig")
        self.assertFileExists(macs2, "control_lambda")
        self.assertFileExists(macs2, "control_lambda_bigwig")

    @tag_process("macs2-rose2-batch")
    def test_macs2_rose2_batch(self):
        with self.preparation_stage():
            inputs = {
                "src": "macs2/input/SRR5675973_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            case_1 = self.run_process("upload-bam", inputs)

            inputs = {
                "src": "macs2/input/SRR5675974_chr17.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            background_1 = self.run_process("upload-bam", inputs)

            promoters = self.run_process(
                "upload-bed",
                {
                    "src": "macs2/input/promoter_regions.bed",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

            collection = Collection.objects.create(
                name="Test collection", contributor=self.contributor
            )

            rel_type_background = RelationType.objects.get(name="background")

            background = Relation.objects.create(
                contributor=self.contributor,
                collection=collection,
                type=rel_type_background,
                category="Background",
            )

            RelationPartition.objects.create(
                relation=background, entity=case_1.entity, label="case"
            )
            RelationPartition.objects.create(
                relation=background, entity=background_1.entity, label="background"
            )

            self.assertEqual(
                background_pairs(
                    [
                        {"__id": case_1.id, "__type": case_1.process.type},
                        {"__id": background_1.id, "__type": background_1.process.type},
                    ]
                ),
                [(case_1.id, background_1.id)],
            )

        self.run_process(
            "macs2-rose2-batch",
            {"alignments": [case_1.id, background_1.id], "promoter": promoters.id},
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        rose2 = Data.objects.filter(process__slug="rose2").last()

        # remove changing lines from the rose2 output
        def filter_created(line):
            return line.startswith(b"#Created")

        self.assertFile(
            rose2,
            "all_enhancers",
            "macs2/output/rose2_enhancer_table.txt",
            file_filter=filter_created,
        )

    @tag_process("chipqc")
    def test_chipqc(self):
        with self.preparation_stage():

            def set_sample_name(data, sample_name):
                """Set sample name."""
                sample = Sample.objects.get(data=data)
                sample.name = sample_name
                sample.save()

            alignment = self.run_process(
                "upload-bam",
                {
                    "src": join("chipqc", "input", "SRR5675976_chr17.bam"),
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            macs2_mock = Process.objects.create(
                name="Upload macs2 data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "resolwebio/base:ubuntu-18.04",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:chipseq:callpeak:macs2:",
                entity_type="sample",
                entity_descriptor_schema="sample",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "called_peaks",
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
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "xls" "xls" 0.1 extract
re-save-file called_peaks "${NAME}".xls

re-save species "Homo sapiens"
re-save build "hg19"
""",
                },
            )
            macs2 = self.run_process(
                macs2_mock.slug,
                {"src": join("chipqc", "input", "SRR5675973_chr17_peaks.xls")},
            )
            set_sample_name(macs2, "SRR5675976_chr17.bam")

            macs14_mock = Process.objects.create(
                name="Upload macs14 data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "resolwebio/base:ubuntu-18.04",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:chipseq:callpeak:macs14:",
                entity_type="sample",
                entity_descriptor_schema="sample",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "peaks_bed",
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
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bed" "bed" 0.1 compress
re-save-file peaks_bed "${NAME}".bed.gz

re-save species "Homo sapiens"
re-save build "hg19"
""",
                },
            )
            macs14 = self.run_process(
                macs14_mock.slug,
                {"src": join("chipqc", "input", "macs14_peaks.bed.gz")},
            )
            set_sample_name(macs14, "SRR5675976_chr17.bam")

        chipqc = self.run_process(
            "chipqc",
            {
                "alignment": alignment.id,
                "peaks": macs2.id,
            },
        )
        self.assertFileExists(chipqc, "ccplot")
        self.assertFileExists(chipqc, "coverage_histogram")
        self.assertFileExists(chipqc, "peak_profile")
        self.assertFileExists(chipqc, "peaks_barplot")
        self.assertFileExists(chipqc, "peaks_density_plot")

        # Test annotation
        chipqc_annotation = self.run_process(
            "chipqc",
            {
                "alignment": alignment.id,
                "peaks": macs2.id,
                "calculate_enrichment": True,
            },
        )
        self.assertFileExists(chipqc_annotation, "enrichment_heatmap")

        chipqc_macs14 = self.run_process(
            "chipqc",
            {
                "alignment": alignment.id,
                "peaks": macs14.id,
            },
        )
        self.assertFileExists(chipqc_macs14, "ccplot")
        self.assertFileExists(chipqc_macs14, "peak_profile")

    @tag_process(
        "workflow-subsample-bwa-aln-single", "workflow-subsample-bwa-aln-paired"
    )
    def test_subsample_bwa_align(self):
        with self.preparation_stage():
            base_input = Path("test_subsample_align")
            input_folder = base_input / "input"
            output_folder = base_input / "output"
            with self.preparation_stage():
                ref_seq = self.prepare_ref_seq(
                    fn=str(input_folder / "g.en ome.fasta.gz"),
                    species="Homo sapiens",
                    build="hg38",
                )
                reads = self.prepare_reads()
                reads_paired = self.prepare_paired_reads(
                    mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                    mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
                )

            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            self.run_process(
                "workflow-subsample-bwa-aln-single",
                {
                    "reads": reads.id,
                    "genome": bwa_index.id,
                    "downsampling": {"n_reads": 10},
                },
            )

            alignment_single = Data.objects.filter(process__slug="alignment-bwa-aln")[0]
            self.assertFile(
                alignment_single,
                "stats",
                output_folder / "downsampled_aligned_stats.txt",
            )

            self.run_process(
                "workflow-subsample-bwa-aln-paired",
                {
                    "reads": reads_paired.id,
                    "genome": bwa_index.id,
                    "downsampling": {"n_reads": 15},
                },
            )

            for data in Data.objects.all():
                self.assertStatus(data, Data.STATUS_DONE)

            alignment_paired = Data.objects.filter(
                process__slug="alignment-bwa-aln"
            ).last()
            self.assertFile(
                alignment_paired,
                "stats",
                output_folder / "downsampled_aligned_paired_stats.txt",
            )
