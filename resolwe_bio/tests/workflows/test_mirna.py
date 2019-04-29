# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe.flow.models import Data
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class MicroRNATestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process('workflow-mirna')
    def test_mirna_workflow(self):
        # Prepare data for aligning the reads with bowtie2 and annotation file for featureCounts.
        with self.preparation_stage():
            inputs = {
                'src': 'genome_rsem.fa.gz',
                'species': 'Homo sapiens',
                'build': 'fake_genome_RSEM',
                'advanced': {
                    'bowtie_index': 'genome_rsem_bt_index.tar.gz',
                    'bowtie2_index': 'genome_rsem_bt2_index.tar.gz',
                    'bwa_index': 'genome_rsem_bwa_index.tar.gz',
                    'hisat2_index': 'genome_rsem_hisat2_index.tar.gz',
                    'subread_index': 'genome_rsem_subread_index.tar.gz',
                },
            }
            genome = self.run_process('upload-genome', inputs)
            single_reads = self.prepare_reads(['reads rsem.fq.gz'])
            annotation = self.prepare_annotation('annotation_rsem.gtf.gz', species='Homo sapiens',
                                                 build='fake_genome_RSEM')

        inputs = {
            'genome': genome.pk,
            'reads': single_reads.pk,
            'annotation': annotation.pk,
            'feature_class': 'exon',
        }

        # Run process and assert.
        self.run_process('workflow-mirna', inputs)
        workflow = Data.objects.last()

        # check featureCount summary
        self.assertFile(workflow, 'rc', 'mirna_featurecounts_rc.tab.gz', compression='gzip')
        self.assertFile(workflow, 'fpkm', 'mirna_featurecounts_fpkm.tab.gz', compression='gzip')
        self.assertFile(workflow, 'exp', 'mirna_featurecounts_tpm.tab.gz', compression='gzip')
