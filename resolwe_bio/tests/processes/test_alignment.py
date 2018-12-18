# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class AlignmentProcessorTestCase(KBBioProcessTestCase):

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

        # Values for alignment options are default according to the documentation. However, L may not be set
        # correctly as there is some incongruency. See https://github.com/BenLangmead/bowtie2/issues/215
        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"},
            'alignment_options': {'N': 0, 'gbar': 4, 'L': 22, 'mp': '6', 'rdg': '5,3', 'rfg': '5,3',
                                  'score_min': 'L,-0.6,-0.6'}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_reads_report.txt')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"},
            'alignment_options': {'N': 0, 'gbar': 4, 'L': 22, 'mp': '6,2', 'rdg': '5,3', 'rfg': '5,3',
                                  'score_min': 'L,-0.6,-0.6'}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_paired_end_report.txt')

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
            'trimming': {'trim_iter': 2, 'trim_nucl': 4},
            'reporting': {'rep_mode': "def"},
            'PE_options': {'use_se': True},
            'alignment_options': {'N': 0, 'gbar': 4, 'L': 22, 'mp': '6,2', 'rdg': '5,3', 'rfg': '5,3',
                                  'score_min': 'L,-0.6,-0.6'}
        }
        aligned_reads = self.run_process('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_use_SE_report.txt')

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

        # Test genome indexing without the input annotation file
        inputs_gtf = {'genome2': star_index_fasta.id}
        self.run_process('alignment-star-index', inputs_gtf)

        inputs_gff3 = {'annotation': annotation_gff3.id, 'genome': genome.id}
        star_index = self.run_process('alignment-star-index', inputs_gff3)
        self.assertAlmostEqual(star_index.output['index']['size'], 1566163859, delta=5)

    @with_resolwe_host
    @tag_process('alignment-star')
    def test_star(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])
            paired_reads = self.prepare_paired_reads(mate1=['SRR2124780_1 1k.fastq.gz'],
                                                     mate2=['SRR2124780_2 1k.fastq.gz'])
            annotation = self.prepare_annotation(fn='HS chr21_short.gtf.gz', source='ENSEMBL',
                                                 species='Homo sapiens', build='GRCh38_ens90')
            inputs = {
                'src': 'HS chr21_ensembl.fa.gz',
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
                'source': 'ENSEMBL',
            }
            star_index_fasta = self.run_process('upload-fasta-nucl', inputs)
            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}

            star_index = self.run_process('alignment-star-index', inputs)

            star_index_wo_annot = self.run_process('alignment-star-index', {'genome2': star_index_fasta.id})

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        inputs = {
            'genome': star_index.id,
            'reads': reads.id,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True,
            },
            'two_pass_mapping': {
                'two_pass_mode': True,
            },
        }
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star_single.tab.gz', compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'star_expression_single.tab.gz', compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs = {
            'genome': star_index_wo_annot.id,
            'reads': reads.id,
            'annotation': annotation.id,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True,
            },
            'two_pass_mapping': {
                'two_pass_mode': True,
            },
        }
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star_single.tab.gz', compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'star_expression_single.tab.gz', compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs = {
            'genome': star_index.id,
            'reads': paired_reads.id,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True,
            },
            'two_pass_mapping': {
                'two_pass_mode': True,
            },
        }
        aligned_reads = self.run_process('alignment-star', inputs)
        self.assertFile(aligned_reads, 'bigwig', 'bigwig_star_ens_paired.bw')
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star_paired.tab.gz', compression='gzip')
        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'star_expression_paired.tab.gz', compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFile(exp, 'exp_set', 'star_out_exp_set.txt.gz', compression='gzip')
        self.assertJSON(exp, exp.output['exp_set_json'], '', 'star_exp_set.json.gz')

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
        self.assertFile(aligned_reads, 'bigwig', 'hisat2_paired_bigwig.bw')
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
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            genome_2 = self.run_process('upload-genome', inputs)

            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                     mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])

        inputs = {'genome': genome.id, 'reads': reads.id}
        aligned_reads = self.run_process('alignment-subread', inputs)
        self.assertFile(aligned_reads, 'stats', 'subread_reads_report.txt')
        self.assertFile(aligned_reads, 'bam', 'subread_single.bam')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {'genome': genome_2.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-subread', inputs)
        self.assertFile(aligned_reads, 'stats', 'subread_paired_reads_report.txt')
        self.assertFile(aligned_reads, 'bam', 'subread_paired.bam')

    @tag_process('alignment-hisat2')
    def test_hisat2_bigwig(self):
        with self.preparation_stage():
            paired_reads_ucsc = self.prepare_paired_reads(mate1=['SRR2124780_1 1k.fastq.gz'],
                                                          mate2=['SRR2124780_2 1k.fastq.gz'])
            paired_reads_ncbi = self.prepare_paired_reads(mate1=['GRCh38.p12_NCBIchr21_R1.fq.gz'],
                                                          mate2=['GRCh38.p12_NCBIchr21_R2.fq.gz'])

            inputs = {
                'src': 'hg38_chr21_9M.fa.gz',
                'species': 'Homo sapiens',
                'build': 'hg38',
            }
            genome_ucsc = self.run_process('upload-genome', inputs)

            inputs = {
                'src': 'GRCh38.p12_NCBIchr21_9M.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'GRCh38.p12',
            }
            genome_ncbi = self.run_process('upload-genome', inputs)

        inputs = {
            'genome': genome_ucsc.id,
            'reads': paired_reads_ucsc.id,
        }
        aligned_reads = self.run_process('alignment-hisat2', inputs)
        self.assertFile(aligned_reads, 'bigwig', 'hisat2_paired_ucsc_bigwig.bw')

        inputs = {
            'genome': genome_ncbi.id,
            'reads': paired_reads_ncbi.id,
        }
        aligned_reads = self.run_process('alignment-hisat2', inputs)
        self.assertFile(aligned_reads, 'bigwig', 'hisat2_paired_ncbi_bigwig.bw')

    @tag_process('alignment-hisat2')
    def test_empty_bam(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads(['reads-map-to-nowhere.fastq'])

        aligned_reads = self.run_process('alignment-hisat2', input_={
            'reads': reads.id,
            'genome': genome.id,
        })
        self.assertEqual(aligned_reads.process_warning, ['Bam file has no entries. No bigWig file will be made.'])

    @tag_process('alignment-hisat2')
    def test_no_bigwig_mappings(self):
        """Use dicty genome and reads but declare it as human so no mapping to UCSC can happen."""
        with self.preparation_stage():
            reads = self.prepare_reads()
            genome = self.run_process('upload-genome', {
                'src': 'genome.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'GRCh38.p12',
            })

        aligned_reads = self.run_process('alignment-hisat2', {
            'reads': reads.id,
            'genome': genome.id,
        })
        msg = 'Neither of the chromosomes in the input file has a valid UCSC pair. No mapping will be done.'
        self.assertEqual(aligned_reads.process_warning, [msg])
