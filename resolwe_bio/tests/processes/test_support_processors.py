# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class SupportProcessorTestCase(BioProcessTestCase):

    @tag_process('reference_compatibility')
    def test_reference_compatibility(self):
        with self.preparation_stage():
            inputs = {
                'src': 'sp_test.fasta',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            genome = self.run_process('upload-genome', inputs)

            mapping = self.prepare_bam()
            annotation = self.prepare_annotation()

        inputs = {'reference': genome.pk, 'bam': mapping.pk, 'annot': annotation.pk}
        compatibility_test = self.run_process('reference_compatibility', inputs)
        self.assertFile(compatibility_test, 'report_file', 'sp_test_compatibility_report.txt')

    @tag_process('bam-split')
    def test_bam_split(self):
        with self.preparation_stage():
            bam = self.prepare_bam(fn='hybrid.bam', species='Mus musculus',
                                   build='mm10_dm6')

            header = self.run_process('upload-header-sam', {'src': 'mm10_header.sam'})
            header2 = self.run_process('upload-header-sam', {'src': 'dm6_header.sam'})

        inputs = {
            'bam': bam.id,
        }
        bam1 = self.run_process('bam-split', inputs)
        bam2 = Data.objects.last()

        self.assertFile(bam1, 'bam', 'hybrid_mm10.bam')
        self.assertFile(bam1, 'bai', 'hybrid_mm10.bam.bai')
        self.assertFields(bam1, 'species', 'Mus musculus')
        self.assertFields(bam1, 'build', 'mm10')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')
        self.assertFields(bam2, 'species', 'Drosophila melanogaster')
        self.assertFields(bam2, 'build', 'dm6')

        inputs['header'] = header.id
        inputs['header2'] = header2.id
        bam1 = self.run_process('bam-split', inputs)
        bam2 = Data.objects.last()

        self.assertFile(bam1, 'bam', 'hybrid_mm10.bam')
        self.assertFile(bam1, 'bai', 'hybrid_mm10.bam.bai')
        self.assertFields(bam1, 'species', 'Mus musculus')
        self.assertFields(bam1, 'build', 'mm10')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')
        self.assertFields(bam2, 'species', 'Drosophila melanogaster')
        self.assertFields(bam2, 'build', 'dm6')

    @tag_process('feature_location')
    def test_feature_location(self):
        with self.preparation_stage():
            inputs = {
                'src': 'mm10_small.gtf.gz',
                'source': 'UCSC',
                'species': 'Mus musculus',
                'build': 'mm10'
            }
            annotation = self.run_process('upload-gtf', inputs)

        inputs = {'annotation': annotation.pk,
                  'feature_type': 'exon',
                  'id_type': 'transcript_id',
                  'summarize_exons': True}
        features = self.run_process('feature_location', inputs)
        self.assertJSON(features, features.output['feature_location'], '', 'feature_locations.json.gz')

    @tag_process('gff-to-gtf')
    def test_gff_to_gtf(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation_gff()

        gff_to_gtf = self.run_process('gff-to-gtf', {'annotation': annotation.id})
        self.assertFile(gff_to_gtf, 'annot', 'gff_to_gtf_annotation.gtf')

    @tag_process('archive-samples')
    def test_ars(self):
        with self.preparation_stage():
            txt_file = self.run_process('upload-file', {'src': '56G_masterfile_test.txt'})
            bam_input = {
                'src': 'bamplot_alignment.bam',
                'species': 'Mus musculus',
                'build': 'GRCh38 _ens90',
            }
            bam = self.run_process('upload-bam', bam_input)

            read_inputs = {'src': ['rRNA forw.fastq.gz', 'rRNA_rew.fastq.gz']}
            reads = self.run_process('upload-fastq-single', read_inputs)

            vcf_input = {
                'src': 'igv_human.lf.vcf',
                'species': 'Homo sapiens',
                'build': 'b37',
            }
            vcf = self.run_process('upload-variants-vcf', vcf_input)

            expression_1 = self.prepare_expression(f_exp='exp_1_rc.tab.gz', f_type="RC")
            expression_2 = self.prepare_expression(f_exp='exp_2_rc.tab.gz', f_type="RC")
            expression_5 = self.prepare_expression(f_exp='exp_5_rc.tab.gz', f_type="RC")
            expression_3 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")
            expression_4 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        self.run_process('archive-samples', {
            'data': [txt_file.id, bam.id, reads.id, vcf.id, expression_1.id, expression_2.id, expression_5.id,
                     expression_4.id, expression_3.id],
            'fields': ['file', 'bam', 'bai', 'fastq', 'fastqc_url', 'fastqc_archive', 'vcf', 'exp']})

    @tag_process('prepare-geo-chipseq')
    def test_prepare_geo_chipseq(self):
        with self.preparation_stage():
            reads_1 = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])
            reads_2 = self.prepare_reads()
            reads_3 = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])

            case_bam = self.prepare_bam(fn='macs14_case.bam', species='Homo sapiens',
                                        build='hg19')
            control_bam = self.prepare_bam(fn='macs14_control.bam', species='Homo sapiens',
                                           build='hg19')

            inputs = {"treatment": case_bam.id,
                      "control": control_bam.id}
            macs14_1 = self.run_process("macs14", inputs)
            reads_1.entity_set.first().data.add(macs14_1)

            # Run macs14 without control/background sample
            del inputs['control']
            macs14_2 = self.run_process("macs14", inputs)
            reads_3.entity_set.first().data.add(macs14_2)

            sample_1 = reads_1.entity_set.last().name
            sample_2 = reads_2.entity_set.last().name

        inputs = {
            "reads": [reads_1.id, reads_2.id, reads_3.id],
            "macs14": [macs14_1.id, macs14_2.id],
            "relations": ["{}:{}".format(sample_1, sample_2)],
            "name": "prepare_geo"
        }
        prepare_geo_chipseq = self.run_process("prepare-geo-chipseq", inputs)

        self.assertFile(prepare_geo_chipseq, 'table', 'prepare_geo_ChIP-Seq.txt')

    @tag_process('prepare-geo-rnaseq')
    def test_prepare_geo_rnaseq(self):
        with self.preparation_stage():
            reads_1 = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])
            reads_2 = self.prepare_reads()

            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

            # Delete expression samples
            expression_1.entity_set.all().delete()
            expression_2.entity_set.all().delete()

            # Add expressions to reads samples
            reads_1.entity_set.first().data.add(expression_1)
            reads_2.entity_set.first().data.add(expression_2)

        inputs = {
            "reads": [reads_1.id, reads_2.id],
            "expressions": [expression_1.pk, expression_2.pk],
            "name": "prepare_geo"
        }
        prepare_geo_rnaseq = self.run_process("prepare-geo-rnaseq", inputs)

        self.assertFile(prepare_geo_rnaseq, 'table', 'prepare_geo_RNA-Seq.txt')

    @tag_process('library-strandedness')
    def test_library_strandedness(self):
        with self.preparation_stage():
            cds = self.run_process('upload-fasta-nucl', {'src': 'salmon_cds.fa.gz'})

            inputs = {
                'nucl': cds.id,
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            }
            salmon_index = self.run_process('salmon-index', inputs)

            single_reads = self.prepare_reads(['reads rsem.fq.gz'])
            paired_reads = self.prepare_paired_reads(mate1=['reads rsem.fq.gz'], mate2=['reads rsem2.fq.gz'])

        single_input = {
            'reads': single_reads.id,
            'salmon_index': salmon_index.id,
        }

        lib_strandedness_single = self.run_process('library-strandedness', single_input)
        self.assertFields(lib_strandedness_single, 'strandedness', 'U')
        self.assertFields(lib_strandedness_single, 'fragment_ratio', 1.0)

        paired_input = {
            'reads': paired_reads.id,
            'salmon_index': salmon_index.id,
        }

        lib_strandedness_paired = self.run_process('library-strandedness', paired_input)
        self.assertFields(lib_strandedness_paired, 'strandedness', 'IU')
        self.assertFields(lib_strandedness_paired, 'fragment_ratio', 1.0)

    @tag_process('prepare-annotation-dexseq')
    def test_prepare_dexseq_gtf(self):
        with self.preparation_stage():
            inputs = {
                'src': 'annotation_ensembl_mus_musculus_short.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Mus musculus',
                'build': 'GRCm38.p5',
            }
            annotation = self.run_process('upload-gtf', inputs)

        inputs = {'annotation': annotation.pk}
        prepare_dexseq = self.run_process('prepare-annotation-dexseq', inputs)
        self.assertFileExists(prepare_dexseq, 'annot')
