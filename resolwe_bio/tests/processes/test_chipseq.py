# pylint: disable=missing-docstring
from os.path import join

from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class ChipSeqProcessorTestCase(BioProcessTestCase):

    @tag_process('chipseq-peakscore', 'chipseq-genescore')
    def test_chipseq(self):
        with self.preparation_stage():
            inputs = {
                'src': 'chip_seq_control.bam',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            control_bam = self.run_process("upload-bam", inputs)

            inputs = {
                'src': 'chip_seq_case.bam',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            case_bam = self.run_process("upload-bam", inputs)

            inputs = {
                'src': 'chip_seq.bed',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            bed = self.run_process('upload-bed', inputs)

            inputs = {
                'case': case_bam.pk,
                'control': control_bam.pk,
                'settings': {'nomodel': True,
                             'pvalue': 0.001,
                             'slocal': 2000,
                             'extsize': 100,
                             'call_summits': True}}
            macs2 = self.run_process("macs2-callpeak", inputs)

        inputs = {
            'peaks': macs2.pk,
            'bed': bed.pk}
        peak_score = self.run_process('chipseq-peakscore', inputs)
        self.assertFile(peak_score, 'peak_score', 'chip_seq_peakscore_genomicContext')

        inputs = {'peakscore': peak_score.id}
        gene_score = self.run_process('chipseq-genescore', inputs)
        self.assertFile(gene_score, 'genescore', 'chip_seq_geneScore.xls')

    @tag_process('macs14')
    def test_macs14(self):
        with self.preparation_stage():
            inputs = {
                'src': 'macs14_control.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            control_bam = self.run_process("upload-bam", inputs)

            inputs = {
                'src': 'macs14_case.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            case_bam = self.run_process("upload-bam", inputs)

        inputs = {"treatment": case_bam.id,
                  "control": control_bam.id}
        macs14 = self.run_process("macs14", inputs)

        self.assertFields(macs14, 'species', 'Homo sapiens')
        self.assertFields(macs14, 'build', 'hg19')
        self.assertFile(macs14, 'peaks_bed', 'macs14_peaks.bed.gz')
        self.assertFile(macs14, 'peaks_bigbed_igv_ucsc', 'macs14_peaks.bb')
        self.assertFile(macs14, 'peaks_tbi_jbrowse', 'macs14_peaks.gz.tbi')
        self.assertFile(macs14, 'summits_tbi_jbrowse', 'macs14_summits.gz.tbi')
        self.assertFile(macs14, 'treat_bigwig', 'macs14_treat.bw')

    @tag_process('macs2-callpeak', 'archive-samples')
    def test_macs2(self):
        with self.preparation_stage():
            inputs = {
                'src': 'macs2_case.bam',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            case_bam = self.run_process("upload-bam", inputs)

            inputs = {
                'src': 'promoters_chr21_subregion.bed',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            promoters = self.run_process('upload-bed', inputs)

            inputs = {
                'alignment': case_bam.id,
            }
            prepeak = self.run_process("qc-prepeak", inputs)

        inputs = {
            'case': case_bam.id,
            'promoter': promoters.id,
            'settings': {
                'extsize': 298,
                'nomodel': True,
                'bedgraph': True,
            }
        }
        macs2 = self.run_process("macs2-callpeak", inputs)

        self.assertFields(macs2, 'species', 'Homo sapiens')
        self.assertFields(macs2, 'build', 'hg19')
        self.assertFile(macs2, 'chip_qc', 'postpeak_qc_report.txt')
        self.assertFileExists(macs2, 'called_peaks')
        self.assertFileExists(macs2, 'narrow_peaks')
        self.assertFileExists(macs2, 'narrow_peaks_bigbed_igv_ucsc')
        self.assertFileExists(macs2, 'summits')
        self.assertFileExists(macs2, 'summits_tbi_jbrowse')
        self.assertFileExists(macs2, 'summits_bigbed_igv_ucsc')
        self.assertFileExists(macs2, 'treat_pileup')
        self.assertFileExists(macs2, 'treat_pileup_bigwig')
        self.assertFileExists(macs2, 'control_lambda')
        self.assertFileExists(macs2, 'control_lambda_bigwig')

        # With "qc-prepeak" data object as input
        inputs = {
            'case': case_bam.id,
            'case_prepeak': prepeak.id,
            'promoter': promoters.id,
            'settings': {
                'bedgraph': True,
            }
        }
        macs2 = self.run_process("macs2-callpeak", inputs)

        self.assertFields(macs2, 'species', 'Homo sapiens')
        self.assertFields(macs2, 'build', 'hg19')
        self.assertFileExists(macs2, 'chip_qc')
        self.assertFileExists(macs2, 'called_peaks')
        self.assertFileExists(macs2, 'narrow_peaks')
        self.assertFileExists(macs2, 'narrow_peaks_bigbed_igv_ucsc')
        self.assertFileExists(macs2, 'summits')
        self.assertFileExists(macs2, 'summits_tbi_jbrowse')
        self.assertFileExists(macs2, 'summits_bigbed_igv_ucsc')

        # Test "archive-samples"' QC report merge
        sample_1 = case_bam.entity_set.first()
        sample_1.data.add(macs2)
        sample_1.name = 'Sample_1'
        sample_1.save()

        inputs = {
            'data': [prepeak.id, macs2.id],
            'fields': ['chip_qc'],
        }
        self.run_process('archive-samples', inputs)

    @skipUnlessLargeFiles('rose2_case.bam', 'rose2_control.bam')
    @tag_process('rose2')
    def test_rose2(self):
        with self.preparation_stage():
            inputs = {
                'src': join('large', 'rose2_case.bam'),
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

            inputs = {
                'src': join('large', 'rose2_control.bam'),
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            control = self.run_process("upload-bam", inputs)

            inputs = {
                'src': 'macs14_chr22.bed',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            macs_peaks = self.run_process('upload-bed', inputs)

            inputs = {
                'src': 'hg19_encode_blacklist_chr22.bed',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            mask = self.run_process('upload-bed', inputs)

        inputs = {
            "input_upload": macs_peaks.id,
            "rankby": bam.id,
            "control": control.id,
            "stitch": 5000,
            "tss": 2500,
            "mask": mask.id
        }
        rose2 = self.run_process("rose2", inputs)

        # remove changing lines from the rose2 output
        def filter_created(line):
            return line.startswith(b'#Created')

        self.assertFile(rose2, 'all_enhancers', 'rose2_enhancer_table.txt', file_filter=filter_created)

    @tag_process('qc-prepeak')
    def test_qc_prepeak(self):
        with self.preparation_stage():
            inputs = {
                'src': 'prepeak_se.bam',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            bam = self.run_process("upload-bam", inputs)

        inputs = {
            "alignment": bam.id,
            'n_sub': 7000,
            'q_treshold': 25,
            'tn5': True,
        }
        prepeak = self.run_process("qc-prepeak", inputs)

        self.assertFields(prepeak, 'species', 'Homo sapiens')
        self.assertFields(prepeak, 'build', 'hg19')
        self.assertFields(prepeak, 'fraglen', 215)
        self.assertFile(prepeak, 'chip_qc', 'prepeak_se_qc_report.txt')
        self.assertFileExists(prepeak, 'tagalign')

        with self.preparation_stage():
            inputs = {
                'src': 'prepeak_pe.bam',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            bam = self.run_process("upload-bam", inputs)

        inputs = {
            "alignment": bam.id,
            'n_sub': 7000,
            'q_treshold': 25,
            'tn5': True,
        }
        prepeak = self.run_process("qc-prepeak", inputs)

        self.assertFields(prepeak, 'species', 'Homo sapiens')
        self.assertFields(prepeak, 'build', 'hg19')
        self.assertFields(prepeak, 'fraglen', 225)
        self.assertFile(prepeak, 'chip_qc', 'prepeak_pe_qc_report.txt')
        self.assertFileExists(prepeak, 'tagalign')
