# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase
from resolwe_bio.models import Sample


class AlignmentProcessorTestCase(BioProcessTestCase):

    @tag_process('alignment-bowtie')
    def test_bowtie(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads_single = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {
            'genome': genome.id,
            'reads': reads_single.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'r': "-a -m 1 --best --strata"}
        }
        alignment = self.run_process('alignment-bowtie', inputs)
        self.assertFile(alignment, 'stats', 'bowtie_single_reads_report.tab.gz', compression='gzip')
        self.assertFields(alignment, 'species', 'Dictyostelium discoideum')
        self.assertFields(alignment, 'build', 'dd-05-2009')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'r': "-a -m 1 --best --strata"},
            'use_se': True
        }
        alignment = self.run_process('alignment-bowtie', inputs)
        self.assertFile(alignment, 'stats', 'bowtie_use_SE_report.tab.gz', compression='gzip')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'r': "-a -m 1 --best --strata"}
        }
        alignment = self.run_process('alignment-bowtie', inputs)
        self.assertFile(alignment, 'stats', 'bowtie_paired_reads_report.tab.gz', compression='gzip')

    @tag_process('alignment-bowtie2')
    def test_bowtie2(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_paired_end_report.txt')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"},
            'PE_options': {'use_se': True}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_use_SE_report.txt')

    @tag_process('alignment-tophat2')
    def test_tophat(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

            annotation = self.prepare_annotation_gff()

        inputs = {
            'genome': genome.id,
            'reads': reads.id,
            'annotation': annotation.id,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_process('alignment-tophat2', inputs)
        self.assertFile(aligned_reads, 'stats', 'tophat_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'annotation': annotation.id,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_process('alignment-tophat2', inputs)
        self.assertFile(aligned_reads, 'stats', 'tophat_paired_reads_report.txt')

    @tag_process('alignment-star-index')
    def test_star_index(self):
        with self.preparation_stage():
            annotation_gtf = self.prepare_annotation(fn='HS chr21_short.gtf.gz', source='UCSC',
                                                     species='Homo sapiens', build='hg19')
            star_index_fasta = self.prepare_adapters(fn='HS chr21_ensembl.fa.gz')
            genome = self.prepare_genome()
            annotation_gff3 = self.prepare_annotation_gff()

        inputs_gtf = {'annotation': annotation_gtf.id, 'genome2': star_index_fasta.id}
        self.run_process('alignment-star-index', inputs_gtf)

        inputs_gff3 = {'annotation': annotation_gff3.id, 'genome': genome.id}
        self.run_process('alignment-star-index', inputs_gff3)

    @tag_process('alignment-star-index', 'alignment-star')
    def test_star(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])
            paired_reads = self.prepare_paired_reads(mate1=['SRR2124780_1 1k.fastq.gz'],
                                                     mate2=['SRR2124780_2 1k.fastq.gz'])
            annotation = self.prepare_annotation(fn='HS chr21_short.gtf.gz', source='UCSC',
                                                 species='Homo sapiens', build='hg19')
            inputs = {
                'src': 'HS chr21_ensembl.fa.gz',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            star_index_fasta = self.run_process('upload-fasta-nucl', inputs)
            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}

            star_index = self.run_process('alignment-star-index', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        inputs = {
            'genome': star_index.id,
            'reads': reads.id,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True}}
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star_single.tab.gz', compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'hg19')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'star_expression_single.tab.gz', compression='gzip')
        self.assertFields(exp, 'source', 'UCSC')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'hg19')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs = {
            'genome': star_index.id,
            'reads': paired_reads.id,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True}}
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star_paired.tab.gz', compression='gzip')
        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'star_expression_paired.tab.gz', compression='gzip')
        self.assertFields(exp, 'source', 'UCSC')

    @tag_process('alignment-bwa-aln')
    def test_bwa_bt(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {'genome': genome.id, 'reads': reads.id}
        aligned_reads = self.run_process('alignment-bwa-aln', inputs)
        self.assertFile(aligned_reads, 'stats', 'bwa_bt_reads_report.txt')

        inputs = {'genome': genome.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-bwa-aln', inputs)
        self.assertFile(aligned_reads, 'stats', 'bwa_bt_paired_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

    @tag_process('alignment-bwa-sw')
    def test_bwa_sw(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {'genome': genome.id, 'reads': reads.id}
        aligned_reads = self.run_process('alignment-bwa-sw', inputs)
        self.assertFile(aligned_reads, 'bam', 'bwa_sw_reads_mapped.bam')
        self.assertFile(aligned_reads, 'stats', 'bwa_sw_reads_report.txt')

        inputs = {'genome': genome.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-bwa-sw', inputs)
        self.assertFile(aligned_reads, 'bam', 'bwa_sw_paired_reads_mapped.bam')
        self.assertFile(aligned_reads, 'stats', 'bwa_sw_paired_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

    @tag_process('alignment-bwa-mem')
    def test_bwa_mem(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {'genome': genome.id, 'reads': reads.id}
        aligned_reads = self.run_process('alignment-bwa-mem', inputs)
        self.assertFile(aligned_reads, 'stats', 'bwa_mem_reads_report.txt')

        inputs = {'genome': genome.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-bwa-mem', inputs)
        self.assertFile(aligned_reads, 'stats', 'bwa_mem_paired_reads_report.txt')
        self.assertFile(aligned_reads, 'unmapped', 'bwa_mem_unmapped_reads.fastq.gz', compression='gzip', sort=True)
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

    @tag_process('alignment-hisat2')
    def test_hisat2(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            sample = Sample.objects.get(data=reads)
            sample.name = 'Single reads'
            sample.save()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])
            sample_paired = Sample.objects.get(data=reads_paired)
            sample_paired.name = 'Paired-end reads'
            sample_paired.save()

        inputs = {
            'genome': genome.id,
            'reads': reads.id}
        aligned_reads = self.run_process('alignment-hisat2', inputs)
        self.assertFile(aligned_reads, 'stats', 'hisat2_report.txt')
        self.assertFile(aligned_reads, 'unmapped_f', 'hisat2_unmapped.fastq.gz', compression='gzip', sort=True)
        self.assertFileExists(aligned_reads, 'splice_junctions')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-hisat2', inputs)
        self.assertFile(aligned_reads, 'stats', 'hisat2_paired_report.txt')
        self.assertFile(aligned_reads, 'unmapped_f', 'hisat2_unmapped_1.fastq.gz', compression='gzip', sort=True)
        self.assertFile(aligned_reads, 'unmapped_r', 'hisat2_unmapped_2.fastq.gz', compression='gzip', sort=True)
        self.assertFileExists(aligned_reads, 'splice_junctions')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

    @tag_process('alignment-subread')
    def test_subread(self):
        with self.preparation_stage():
            genome = self.prepare_genome()

            inputs = {
                'src': 'my.strange.genome name$.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            genome_2 = self.run_process('upload-genome', inputs)

            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {'genome': genome.id, 'reads': reads.id}
        aligned_reads = self.run_process('alignment-subread', inputs)
        self.assertFile(aligned_reads, 'stats', 'subread_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {'genome': genome_2.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-subread', inputs)
        self.assertFile(aligned_reads, 'stats', 'subread_paired_reads_report.txt')
