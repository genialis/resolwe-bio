# pylint: disable=missing-docstring
from os.path import join

from resolwe.test import tag_process
from resolwe_bio.utils.test import skipUnlessLargeFiles, BioProcessTestCase


class CoverageProcessorTestCase(BioProcessTestCase):

    @tag_process('coverage-garvan')
    def test_gencover(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

            # GTF import
            inputs = {
                'src': 'annotation_ok.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            annotation_gtf = self.run_process('upload-gtf', inputs)

            # redundant GTF import
            inputs = {
                'src': 'annotation_red.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            annotation_gtf_red = self.run_process('upload-gtf', inputs)

            inputs = {
                'genome': genome.pk,
                'reads': reads.pk,
                'annotation': annotation.pk,
                'PE_options': {
                    'library_type': "fr-unstranded"}}
            aligned_reads = self.run_process('alignment-tophat2', inputs)

            # samtools mapping
            inputs = {
                'mapping': aligned_reads.pk,
                'genome': genome.pk}
            variants = self.run_process('vc-samtools', inputs)

        # Coverage report
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        self.run_process('coverage-garvan', inputs)

        # Missing gene in BAM file test
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf_red.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        exon_cov = self.run_process('coverage-garvan', inputs)
        self.assertFile(exon_cov, 'exon_coverage', 'exons_coverage.txt.gz', compression='gzip')

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    @tag_process('coveragebed')
    def test_amplicon_coverage(self):
        with self.preparation_stage():
            bam_input = {
                'src': join('large', '56GSID_10k_mate1_RG.bam'),
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            bam = self.run_process('upload-bam', bam_input)
            master_file = self.prepare_amplicon_master_file()

        coverage = self.run_process('coveragebed', {
            'alignment': bam.id,
            'master_file': master_file.id,
        })
        self.assertFile(coverage, 'cov_metrics', '56GSID_10k_covMetrics.txt')
        self.assertFile(coverage, 'mean_cov', '56GSID_10k_ampmeancov.covd')
        self.assertFileExists(coverage, 'amplicon_cov')
        self.assertFileExists(coverage, 'covplot_html')
