# pylint: disable=missing-docstring
from resolwe.flow.models import Data, Process
from resolwe.test import tag_process

from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class SupportProcessorTestCase(KBBioProcessTestCase):

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
        self.assertFile(bam1, 'bigwig', 'hybrid_mm10.bw')
        self.assertFields(bam1, 'species', 'Mus musculus')
        self.assertFields(bam1, 'build', 'mm10')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')
        self.assertFile(bam2, 'bigwig', 'hybrid_dm6.bw')
        self.assertFields(bam2, 'species', 'Drosophila melanogaster')
        self.assertFields(bam2, 'build', 'dm6')

        inputs['header'] = header.id
        inputs['header2'] = header2.id
        bam1 = self.run_process('bam-split', inputs)
        bam2 = Data.objects.last()

        self.assertFile(bam1, 'bam', 'hybrid_mm10.bam')
        self.assertFile(bam1, 'bai', 'hybrid_mm10.bam.bai')
        self.assertFile(bam1, 'bigwig', 'hybrid_mm10.bw')
        self.assertFields(bam1, 'species', 'Mus musculus')
        self.assertFields(bam1, 'build', 'mm10')
        self.assertFile(bam2, 'bam', 'hybrid_dm6.bam')
        self.assertFile(bam2, 'bai', 'hybrid_dm6.bam.bai')
        self.assertFile(bam2, 'bigwig', 'hybrid_dm6.bw')
        self.assertFields(bam2, 'species', 'Drosophila melanogaster')
        self.assertFields(bam2, 'build', 'dm6')

    @tag_process('gff-to-gtf')
    def test_gff_to_gtf(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation_gff()

        gff_to_gtf = self.run_process('gff-to-gtf', {'annotation': annotation.id})
        self.assertFile(gff_to_gtf, 'annot', 'gff_to_gtf_annotation.gtf')
        del gff_to_gtf.output['annot_sorted_track_jbrowse']['total_size']  # Non-deterministic output.
        self.assertFields(gff_to_gtf, 'annot_sorted_track_jbrowse', {'refs': ['tracks/annotation'],
                                                                     'file': 'trackList.json'})

    @with_resolwe_host
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

            multiqc = self.run_process('multiqc', {
                'data': [
                    bam.id,
                    reads.id,
                ],
                'advanced': {
                    'dirs': True,
                    'config': True,
                }
            })

        self.run_process('archive-samples', {
            'data': [txt_file.id, bam.id, reads.id, vcf.id, expression_1.id, expression_2.id, expression_5.id,
                     expression_4.id, expression_3.id, multiqc.id],
            'fields': ['file', 'bam', 'bai', 'fastq', 'fastqc_url', 'fastqc_archive', 'vcf', 'exp_set', 'report']})

    @with_resolwe_host
    @tag_process('archive-samples')
    def test_archive_samples_exp_set(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_exp='exp_1_rc.tab.gz', f_type='RC')
            expression_2 = self.prepare_expression(f_exp='exp_2_tpm.tab.gz', f_type='TPM', f_rc='exp_2_rc.tab.gz')

        inputs = {
            'data': [
                expression_1.pk,
                expression_2.pk,
            ],
            'fields': [
                'exp_set',
            ],
        }
        self.run_process('archive-samples', inputs)
        # Structured zip files are not supported by assertFile. When implemented, add here
        # self.assertFile(_, 'archive', 'test_archive_samples_exp_set.zip', compression='zip').

    @with_resolwe_host
    @tag_process('archive-samples')
    def test_archive_samples_exp(self):
        with self.preparation_stage():
            # Upload expression without exp_set output.
            Process.objects.create(
                slug='exp',
                type='data:expression:exp:',
                contributor=self.contributor,
                requirements={'expression-engine': 'jinja'},
                data_name='Upload expression into exp output field',
                entity_type='sample',
                entity_descriptor_schema='sample',
                input_schema=[
                    {
                        'name': 'exp',
                        'type': 'basic:file:',
                    }
                ],
                output_schema=[
                    {
                        'name': 'exp',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'exp_type',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'build',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'species',
                        'type': 'basic:string:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': '''
                        re-import {{ exp.file_temp|default(exp.file) }} {{ exp.file }} "tab|gz" tab 1.0 compress
                        re-save-file exp "${NAME}.tab.gz"
                        re-save exp_type RC
                        re-save build ens
                        re-save species "Homo sapiens"
                    ''',
                },
            )
            inputs = {
                'exp': 'exp_1_rc.tab.gz',
            }
            expression_1 = self.run_process('exp', inputs)
            expression_2 = self.prepare_expression(f_exp='exp_2_tpm.tab.gz', f_type='TPM', f_rc='exp_2_rc.tab.gz')

        inputs = {
            'data': [
                expression_1.pk,
                expression_2.pk,
            ],
            'fields': [
                'exp_set',
            ],
        }
        self.run_process('archive-samples', inputs)
        # Structured zip files are not supported by assertFile. When implemented, add here
        # self.assertFile(_, 'archive', 'test_archive_samples_exp.zip', compression='zip').

    @tag_process('prepare-geo-chipseq')
    def test_prepare_geo_chipseq(self):
        with self.preparation_stage():
            reads_1 = self.prepare_paired_reads(mate1=['fw reads.fastq.gz', 'fw reads_2.fastq.gz'],
                                                mate2=['rw reads.fastq.gz', 'rw reads_2.fastq.gz'])
            reads_2 = self.prepare_reads()
            reads_3 = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])

            macs14_case_bam = self.prepare_bam(fn='macs14_case.bam', species='Homo sapiens',
                                               build='hg19')
            macs14_control_bam = self.prepare_bam(fn='macs14_control.bam', species='Homo sapiens',
                                                  build='hg19')

            macs2_case_bam = self.prepare_bam(fn='macs2/input/SRR5675973_chr17.bam',
                                              species='Homo sapiens', build='hg19')
            macs2_control_bam = self.prepare_bam(fn='macs2/input/SRR5675974_chr17.bam',
                                                 species='Homo sapiens', build='hg19')

            # Run macs14
            inputs = {'treatment': macs14_case_bam.id,
                      'control': macs14_control_bam.id}
            macs14_1 = self.run_process('macs14', inputs)
            macs14_1.entity_set.all().delete()
            reads_1.entity_set.first().data.add(macs14_1)

            macs14_control_bam.entity_set.all().delete()
            reads_2.entity_set.first().data.add(macs14_control_bam)

            # Run macs14 without control/background sample
            del inputs['control']
            macs14_2 = self.run_process('macs14', inputs)
            macs14_2.entity_set.all().delete()
            reads_3.entity_set.first().data.add(macs14_2)

            # Run macs2
            inputs = {
                'case': macs2_case_bam.id,
                "control": macs2_control_bam.id,
                'settings': {
                    'extsize': 298,
                    'nomodel': True,
                    'bedgraph': True,
                },
            }
            macs2_1 = self.run_process('macs2-callpeak', inputs)
            macs2_1.entity_set.all().delete()
            reads_1.entity_set.first().data.add(macs2_1)

            macs2_control_bam.entity_set.all().delete()
            reads_2.entity_set.first().data.add(macs2_control_bam)

            # Run macs2 without control/background sample
            del inputs['control']
            macs2_2 = self.run_process('macs2-callpeak', inputs)
            macs2_2.entity_set.all().delete()
            reads_3.entity_set.first().data.add(macs2_2)

        inputs = {
            'reads': [reads_1.id, reads_2.id, reads_3.id],
            'macs': [macs14_1.id, macs14_2.id],
            'name': 'prepare_geo',
        }
        prepare_geo_chipseq = self.run_process('prepare-geo-chipseq', inputs)

        self.assertFile(prepare_geo_chipseq, 'table', 'prepare_geo_ChIP-Seq_macs14.txt')

        inputs = {
            'reads': [reads_1.id, reads_2.id, reads_3.id],
            'macs': [macs2_1.id, macs2_2.id],
            'name': 'prepare_geo',
        }
        prepare_geo_chipseq = self.run_process('prepare-geo-chipseq', inputs)

        self.assertFile(prepare_geo_chipseq, 'table', 'prepare_geo_ChIP-Seq_macs2.txt')

    @with_resolwe_host
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

    @tag_process('multiqc')
    def test_multiqc(self):
        with self.preparation_stage():
            reads = self.run_processor('upload-fastq-single', {
                'src': ['hs_single bbduk_star_htseq_reads_single.fastq.gz']
            })

            paired_reads = self.prepare_paired_reads(['hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz'],
                                                     ['hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz'])

            filtered_reads = self.run_process('bbduk-paired', {'reads': paired_reads.id})

            bam_samtools = self.run_process('upload-bam-indexed', {
                'src': 'alignment_position_sorted.bam',
                'src2': 'alignment_position_sorted.bam.bai',
                'species': 'Homo sapiens',
                'build': 'hg19',
            })

            annotation = self.run_process('upload-gtf', {
                'src': 'hs annotation.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            })

            genome_fasta = self.run_process('upload-fasta-nucl', {'src': 'hs genome.fasta.gz'})
            star_index = self.run_process('alignment-star-index', {
                'annotation': annotation.id,
                'genome2': genome_fasta.id,
            })

            star_alignment = self.run_process('alignment-star', {
                'genome': star_index.id,
                'reads': paired_reads.id,
            })

            samtools_idxstats = self.run_process('samtools-idxstats', {
                'alignment': star_alignment.id,
            })

            qorts_report = self.run_process('qorts-qc', {
                'alignment': star_alignment.id,
                'annotation': annotation.id,
                'options': {
                    'maxPhredScore': 42,
                },
            })

        multiqc = self.run_process('multiqc', {
            'data': [
                reads.id,
                paired_reads.id,
                filtered_reads.id,
                bam_samtools.id,
                star_alignment.id,
                samtools_idxstats.id,
                qorts_report.id,
            ],
            'advanced': {
                'dirs': True,
                'config': True,
            }
        })
        self.assertFileExists(multiqc, 'report')

    @tag_process('seqtk-sample-single', 'seqtk-sample-paired')
    def test_seqtk_sample(self):
        with self.preparation_stage():
            reads = self.run_processor('upload-fastq-single', {
                'src': ['hs_single bbduk_star_htseq_reads_single.fastq.gz']
            })

            paired_reads = self.prepare_paired_reads(['hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz'],
                                                     ['hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz'])

        inputs_single = {
            'reads': reads.id,
            'n_reads': 42,
            'advanced': {
                'seed': 42,
            }
        }

        seqtk_single = self.run_process('seqtk-sample-single', inputs_single)
        self.assertFiles(seqtk_single, 'fastq', ['seqtk_subsampled_reads_single_end.fastq.gz'],
                         compression='gzip')

        inputs_paired = {
            'reads': paired_reads.id,
            'advanced': {
                'seed': 42,
                'fraction': 0.25,
            }
        }

        seqtk_paired = self.run_process('seqtk-sample-paired', inputs_paired)
        self.assertFiles(seqtk_paired, 'fastq', ['seqtk_subsampled_reads_paired_end_mate1.fastq.gz'],
                         compression='gzip')
        self.assertFiles(seqtk_paired, 'fastq2', ['seqtk_subsampled_reads_paired_end_mate2.fastq.gz'],
                         compression='gzip')

    @with_resolwe_host
    @tag_process('spikein-qc')
    def test_spikein_pairwise(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp='exp1_cpm_ercc_sirv.tab.gz',
                f_type='CPM',
                name='Sample 1',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_2 = self.prepare_expression(
                f_exp='exp2_cpm_ercc_sirv.tab.gz',
                f_type='CPM',
                name='Sample 2',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_3 = self.prepare_expression(
                f_exp='exp3_cpm_ercc_sirv.tab.gz',
                f_type='CPM',
                name='Sample3.txt',  # Test for sample names that might look like filename extensions
                source='ENSEMBL',
                species='Homo sapiens')
            expression_4 = self.prepare_expression(
                f_exp='exp4_cpm_ercc_sirv.tab.gz',
                f_type='CPM',
                name='Sample without ERCC',
                source='ENSEMBL',
                species='Homo sapiens')

        # SIRV Set 3
        sirv_set3 = self.run_process('spikein-qc', {
            'samples': [expression_1.pk, expression_2.pk, expression_3.pk, expression_4.pk],
            'mix': 'sirv_set3',
        })

        self.assertEqual(sirv_set3.process_warning, [
            'All ERCC spike-ins have zero expression in sample Sample without ERCC',
        ])

        self.assertFilesExist(sirv_set3, 'plots')
        expected_names = [item['file'] for item in sirv_set3.output['plots']]
        self.assertEqual(expected_names, [
            "Sample 1 (ERCC spike-in's).png",
            "Sample 2 (ERCC spike-in's).png",
            "Sample3.txt (ERCC spike-in's).png",
        ])

        self.assertFileExists(sirv_set3, 'report')

        self.assertFileExists(sirv_set3, 'report_zip')

    @tag_process('qorts-qc')
    def test_qorts_qc(self):
        with self.preparation_stage():
            alignment = self.run_process('upload-bam', {
                'src': 'qorts/input/hs paired.bam',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            })

            cds = self.run_process('upload-fasta-nucl', {
                'src': 'qorts/input/salmon_cds.fa.gz'
            })
            inputs = {
                'nucl': cds.id,
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            }
            salmon_index = self.run_process('salmon-index', inputs)

            annotation = self.run_process('upload-gtf', {
                'src': 'qorts/input/hs annotation.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            })

        inputs = {
            'alignment': alignment.id,
            'annotation': annotation.id,
            'options': {
                'stranded': 'auto',
                'cdna_index': salmon_index.id,
                'adjustPhredScore': 31,
            },
        }
        qorts_report = self.run_process('qorts-qc', inputs)
        self.assertFileExists(qorts_report, 'plot')
        self.assertFileExists(qorts_report, 'summary')
        self.assertFileExists(qorts_report, 'qorts_data')

    @tag_process('samtools-idxstats')
    def test_samtools_idxstats(self):
        with self.preparation_stage():
            alignment = self.prepare_bam(fn='alignment_position_sorted.bam')

        idxstats = self.run_process('samtools-idxstats', {'alignment': alignment.id})
        self.assertFile(idxstats, 'report', 'samtools_idxstats_report.txt')
