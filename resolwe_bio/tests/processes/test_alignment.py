# pylint: disable=missing-docstring
import os

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import KBBioProcessTestCase


class AlignmentProcessorTestCase(KBBioProcessTestCase):

    @tag_process('bowtie-index')
    def test_bowtie_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn='./test_bowtie/input/g.en ome.fasta.gz',
                species='Dictyostelium discoideum',
                build='dd-05-2009',
            )

        bowtie_index = self.run_process('bowtie-index', {'ref_seq': ref_seq.id})
        self.assertDir(bowtie_index, 'index', os.path.join('test_bowtie', 'output', 'bowtie_index.tar.gz'))
        self.assertFile(bowtie_index, 'fasta', os.path.join('test_bowtie', 'output', 'genome.fasta'))
        self.assertFile(
            bowtie_index,
            'fastagz',
            os.path.join('test_bowtie', 'output', 'genome.fasta.gz'),
            compression='gzip'
        )
        self.assertFile(bowtie_index, 'fai', os.path.join('test_bowtie', 'output', 'genome.fasta.fai'))
        self.assertFields(bowtie_index, 'species', 'Dictyostelium discoideum')
        self.assertFields(bowtie_index, 'build', 'dd-05-2009')

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

    @tag_process('bowtie2-index')
    def test_bowtie2_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn='./test_bowtie2/input/g.en ome.fasta.gz',
                species='Dictyostelium discoideum',
                build='dd-05-2009',
            )

        bowtie2_index = self.run_process('bowtie2-index', {'ref_seq': ref_seq.id})
        self.assertDir(bowtie2_index, 'index', os.path.join('test_bowtie2', 'output', 'bowtie2_index.tar.gz'))
        self.assertFile(bowtie2_index, 'fasta', os.path.join('test_bowtie2', 'output', 'genome.fasta'))
        self.assertFile(
            bowtie2_index,
            'fastagz',
            os.path.join('test_bowtie2', 'output', 'genome.fasta.gz'),
            compression='gzip'
        )
        self.assertFile(bowtie2_index, 'fai', os.path.join('test_bowtie2', 'output', 'genome.fasta.fai'))
        self.assertFields(bowtie2_index, 'species', 'Dictyostelium discoideum')
        self.assertFields(bowtie2_index, 'build', 'dd-05-2009')

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

    @with_resolwe_host
    @tag_process('alignment-star-index', 'alignment-star')
    def test_star(self):
        with self.preparation_stage():
            reads = self.prepare_reads(
                [os.path.join('test_star', 'input', 'hs_single bbduk_star_htseq_reads_single.fastq.gz')]
            )
            paired_reads = self.prepare_paired_reads(
                mate1=[os.path.join('test_star', 'input', 'hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz')],
                mate2=[os.path.join('test_star', 'input', 'hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz')],
            )
            annotation = self.prepare_annotation(
                os.path.join('test_star', 'input', 'hs annotation.gtf.gz'),
                source='ENSEMBL',
                species='Homo sapiens',
                build='GRCh38_ens90',
            )

            inputs = {
                'src': os.path.join('test_star', 'input', 'hs genome.fasta.gz'),
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
            }
            star_index_fasta = self.run_process('upload-fasta-nucl', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        gene_counts_star_single = os.path.join('test_star', 'output', 'gene_counts_star_single.tab.gz')
        star_expression_single = os.path.join('test_star', 'output', 'star_expression_single.tab.gz')

        # prepare genome indices
        star_index = self.run_process('alignment-star-index', {
            'annotation': annotation.id,
            'ref_seq': star_index_fasta.id
        })

        star_index_wo_annot = self.run_process('alignment-star-index', {
            'ref_seq': star_index_fasta.id,
            'source': 'ENSEMBL',
        })

        # test STAR alignment
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
            'detect_chimeric': {
                'chimeric': True,
            },
        }
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_single, compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_single, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs['star_sort'] = True
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_single, compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_single, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs['genome'] = star_index_wo_annot.id
        inputs['annotation'] = annotation.id
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_single, compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_single, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        inputs['star_sort'] = False
        aligned_reads = self.run_process('alignment-star', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_single, compression='gzip')
        self.assertFields(aligned_reads, 'species', 'Homo sapiens')
        self.assertFields(aligned_reads, 'build', 'GRCh38_ens90')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_single, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFields(exp, 'species', 'Homo sapiens')
        self.assertFields(exp, 'build', 'GRCh38_ens90')
        self.assertFields(exp, 'feature_type', 'gene')

        bigwig_star_ens_paired = os.path.join('test_star', 'output', 'bigwig_star_ens_paired.bw')
        gene_counts_star_paired = os.path.join('test_star', 'output', 'gene_counts_star_paired.tab.gz')
        star_expression_paired = os.path.join('test_star', 'output', 'star_expression_paired.tab.gz')
        star_out_exp_set = os.path.join('test_star', 'output', 'star_out_exp_set.txt.gz')
        star_exp_set = os.path.join('test_star', 'output', 'star_exp_set.json.gz')

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
            'star_sort': True,
        }
        aligned_reads = self.run_process('alignment-star', inputs)
        self.assertFile(aligned_reads, 'bigwig', bigwig_star_ens_paired)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_paired, compression='gzip')
        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_paired, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFile(exp, 'exp_set', star_out_exp_set, compression='gzip')
        self.assertJSON(exp, exp.output['exp_set_json'], '', star_exp_set)

        inputs['star_sort'] = False
        aligned_reads = self.run_process('alignment-star', inputs)
        self.assertFile(aligned_reads, 'bigwig', bigwig_star_ens_paired)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(aligned_reads, 'gene_counts', gene_counts_star_paired, compression='gzip')
        exp = Data.objects.last()
        self.assertFile(exp, 'exp', star_expression_paired, compression='gzip')
        self.assertFields(exp, 'source', 'ENSEMBL')
        self.assertFile(exp, 'exp_set', star_out_exp_set, compression='gzip')
        self.assertJSON(exp, exp.output['exp_set_json'], '', star_exp_set)

    @tag_process('bwa-index')
    def test_bwa_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn='./test_bwa/input/g.en ome.fasta.gz',
                species='Dictyostelium discoideum',
                build='dd-05-2009',
            )

        bwa_index = self.run_process('bwa-index', {'ref_seq': ref_seq.id})
        self.assertDir(bwa_index, 'index', os.path.join('test_bwa', 'output', 'bwa_index.tar.gz'))
        self.assertFile(bwa_index, 'fasta', os.path.join('test_bwa', 'output', 'genome.fasta'))
        self.assertFile(
            bwa_index,
            'fastagz',
            os.path.join('test_bwa', 'output', 'genome.fasta.gz'),
            compression='gzip'
        )
        self.assertFile(bwa_index, 'fai', os.path.join('test_bwa', 'output', 'genome.fasta.fai'))
        self.assertFields(bwa_index, 'species', 'Dictyostelium discoideum')
        self.assertFields(bwa_index, 'build', 'dd-05-2009')

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

    @tag_process('hisat2-index')
    def test_hisat2_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn='./test_hisat2/input/g.en ome.fasta.gz',
                species='Dictyostelium discoideum',
                build='dd-05-2009',
            )

        hisat2_index = self.run_process('hisat2-index', {'ref_seq': ref_seq.id})
        self.assertDir(hisat2_index, 'index', os.path.join('test_hisat2', 'output', 'hisat2_index.tar.gz'))
        self.assertFile(hisat2_index, 'fasta', os.path.join('test_hisat2', 'output', 'genome.fasta'))
        self.assertFile(
            hisat2_index,
            'fastagz',
            os.path.join('test_hisat2', 'output', 'genome.fasta.gz'),
            compression='gzip'
        )
        self.assertFile(hisat2_index, 'fai', os.path.join('test_hisat2', 'output', 'genome.fasta.fai'))
        self.assertFields(hisat2_index, 'species', 'Dictyostelium discoideum')
        self.assertFields(hisat2_index, 'build', 'dd-05-2009')

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

    @tag_process('subread-index')
    def test_subread_index(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn='./test_subread/input/g.en ome.fasta.gz',
                species='Dictyostelium discoideum',
                build='dd-05-2009',
            )

        subread_index = self.run_process('subread-index', {'ref_seq': ref_seq.id})
        self.assertDirExists(subread_index, 'index')
        self.assertFile(subread_index, 'fasta', os.path.join('test_subread', 'output', 'genome.fasta'))
        self.assertFile(
            subread_index,
            'fastagz',
            os.path.join('test_subread', 'output', 'genome.fasta.gz'),
            compression='gzip'
        )
        self.assertFile(subread_index, 'fai', os.path.join('test_subread', 'output', 'genome.fasta.fai'))
        self.assertFields(subread_index, 'species', 'Dictyostelium discoideum')
        self.assertFields(subread_index, 'build', 'dd-05-2009')

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
        self.assertFileExists(aligned_reads, 'bam')
        self.assertFields(aligned_reads, 'species', 'Dictyostelium discoideum')
        self.assertFields(aligned_reads, 'build', 'dd-05-2009')

        inputs = {'genome': genome_2.id, 'reads': reads_paired.id}
        aligned_reads = self.run_process('alignment-subread', inputs)
        self.assertFile(aligned_reads, 'stats', 'subread_paired_reads_report.txt')
        self.assertFileExists(aligned_reads, 'bam')

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
