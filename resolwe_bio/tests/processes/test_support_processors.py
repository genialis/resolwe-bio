# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase

from resolwe.flow.models import Data
from resolwe.test import tag_process


class SupportProcessorTestCase(BioProcessTestCase):

    @tag_process('reference_compatibility')
    def test_reference_compatibility(self):
        with self.preparation_stage():
            inputs = {"src": "sp_test.fasta"}
            genome = self.run_process('upload-genome', inputs)

            mapping = self.prepare_bam()
            annotation = self.prepare_annotation()

        inputs = {'reference': genome.pk, 'bam': mapping.pk, 'annot': annotation.pk}
        compatibility_test = self.run_process('reference_compatibility', inputs)
        self.assertFile(compatibility_test, 'report_file', 'sp_test_compatibility_report.txt')

    @tag_process('bam-split')
    def test_bam_split(self):
        with self.preparation_stage():
            bam = self.prepare_bam('hybrid.bam')

            header = self.run_process('upload-header-sam', {'src': 'mm10_header.sam'})
            header2 = self.run_process('upload-header-sam', {'src': 'dm6_header.sam'})

        inputs = {
            'bam': bam.id,
            'organism': 'mm10',
            'organism2': 'dm6'
        }
        self.run_process('bam-split', inputs)

        bam2, bam1 = Data.objects.order_by('-id')[0:2]

        self.assertFile(bam1, 'bam', 'hybrid_mm10.bam')
        self.assertFile(bam1, 'bai', 'hybrid_mm10.bam.bai')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')

        inputs['header'] = header.id
        inputs['header2'] = header2.id
        self.run_process('bam-split', inputs)

        bam2, bam1 = Data.objects.order_by('-id')[0:2]

        self.assertFile(bam1, 'bam', 'hybrid_mm10.bam')
        self.assertFile(bam1, 'bai', 'hybrid_mm10.bam.bai')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')

    @tag_process('feature_location')
    def test_feature_location(self):
        with self.preparation_stage():
            inputs = {'src': 'mm10_small.gtf.gz', 'source': 'UCSC'}
            annotation = self.run_process('upload-gtf', inputs)

        inputs = {'annotation': annotation.pk,
                  'feature_type': 'exon',
                  'id_type': 'transcript_id',
                  'summarize_exons': True}
        features = self.run_process('feature_location', inputs)
        self.assertJSON(features, features.output['feature_location'], '', 'feature_locations.json.gz')

    @tag_process('go-genesets')
    def test_generate_go_genesets(self):
        with self.preparation_stage():
            inputs = {'src': 'go_genesets.mgi.gz', 'source': 'MGI', 'species': 'Mus musculus'}
            gaf = self.run_process('upload-gaf', inputs)

        inputs = {'gaf': gaf.id, 'source': 'MGI_ID'}
        go_genesets = self.run_process('go-genesets', inputs)
        self.assertFields(go_genesets, 'num_genesets', 6)

        last_geneset = Data.objects.last()
        self.assertFile(last_geneset, 'geneset', 'go_geneset.tab.gz', compression='gzip')

    @tag_process('gff-to-gtf')
    def test_gff_to_gtf(self):
        with self.preparation_stage():
            inputs = {'src': 'annotation.gff.gz', 'source': 'DICTYBASE'}
            annotation = self.run_process('upload-gff3', inputs)

        inputs = {'annotation': annotation.pk}
        gff_to_gtf = self.run_process('gff-to-gtf', inputs)
        self.assertFile(gff_to_gtf, 'gtf', 'gff_to_gtf_annotation.gtf')

    @tag_process('archive-samples')
    def test_archive_samples(self):
        with self.preparation_stage():
            txt_file = self.run_process('upload-file', {'src': '56G_masterfile_test.txt'})
            bam = self.run_process('upload-bam', {'src': 'bamplot_alignment.bam'})

            read_inputs = {'src': ['rRNA forw.fastq.gz', 'rRNA_rew.fastq.gz']}
            reads = self.run_process('upload-fastq-single', read_inputs)

            vcf = self.run_process('upload-variants-vcf', {'src': 'igv_human.lf.vcf'})

        self.run_process('archive-samples', {
            'data': [txt_file.id, bam.id, reads.id, vcf.id],
            'fields': ['file', 'bam', 'bai', 'fastq', 'fastqc_url', 'fastqc_archive', 'vcf']})
