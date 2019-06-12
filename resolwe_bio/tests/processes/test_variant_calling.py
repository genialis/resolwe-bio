# pylint: disable=missing-docstring
from os.path import join

from resolwe.flow.models import Data
from resolwe.test import tag_process
from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import skipUnlessLargeFiles, BioProcessTestCase


class VariantCallingTestCase(BioProcessTestCase):

    @tag_process('vc-chemut')
    def test_variant_calling_chemut(self):
        with self.preparation_stage():
            inputs = {
                'src': 'chemut_genome.fasta.gz',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            genome = self.run_process('upload-genome', inputs)

            inputs = {'src1': ['AX4_mate1.fq.gz'],
                      'src2': ['AX4_mate2.fq.gz']}

            parental_reads = self.run_process('upload-fastq-paired', inputs)

            inputs = {'src1': ['CM_mate1.fq.gz'],
                      'src2': ['CM_mate2.fq.gz']}

            mut_reads = self.run_process('upload-fastq-paired', inputs)

            inputs = {'genome': genome.id, 'reads': parental_reads.id}
            align_parental = self.run_process('alignment-bwa-mem', inputs)

            inputs = {'genome': genome.id, 'reads': mut_reads.id}
            align_mut = self.run_process('alignment-bwa-mem', inputs)

        inputs = {
            'genome': genome.pk,
            'parental_strains': [align_parental.id],
            'mutant_strains': [align_mut.id],
            'reads_info': {
                'PL': "Illumina",
                'LB': "x",
                'CN': "def",
                'DT': "2014-08-05"},
            'Varc_param': {'stand_emit_conf': 10, 'stand_call_conf': 30}}

        variants = self.run_process('vc-chemut', inputs)
        self.assertFields(variants, 'build', 'dd-05-2009')
        self.assertFields(variants, 'species', 'Dictyostelium discoideum')

    @tag_process('filtering-chemut')
    def test_filtering_chemut(self):
        with self.preparation_stage():
            vcf_input = {
                'src': 'variant_calling_filtering.vcf.gz',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            variants = self.run_process('upload-variants-vcf', vcf_input)

        inputs = {
            'variants': variants.pk,
            'analysis_type': 'snv',
            'parental_strain': 'AX4',
            'mutant_strain': 'mutant',
            'read_depth': 5}

        filtered_variants = self.run_process('filtering-chemut', inputs)
        self.assertFile(filtered_variants, 'vcf', 'variant_calling_filtered_variants.vcf.gz', compression='gzip')
        self.assertFields(filtered_variants, 'build', 'dd-05-2009')
        self.assertFields(filtered_variants, 'species', 'Dictyostelium discoideum')

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    @tag_process('vc-realign-recalibrate')
    def test_vc_preprocess_bam(self):
        with self.preparation_stage():
            bam_input = {
                'src': join('large', '56GSID_10k_mate1_RG.bam'),
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            bam = self.run_process('upload-bam', bam_input)
            inputs = {
                'src': 'hs_b37_chr2_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            genome = self.run_process('upload-genome', inputs)
            vcf_input = {
                'src': '1000G_phase1.indels.b37_chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            indels = self.run_process('upload-variants-vcf', vcf_input)
            dbsnp_input = {
                'src': 'dbsnp_138.b37.chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            dbsnp = self.run_process('upload-variants-vcf', dbsnp_input)

        inputs = {
            'alignment': bam.id,
            'genome': genome.id,
            'known_indels': [indels.id],
            'known_vars': [dbsnp.id]
        }

        variants = self.run_process('vc-realign-recalibrate', inputs)
        self.assertFields(variants, 'build', 'b37')
        self.assertFields(variants, 'species', 'Homo sapiens')

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    @tag_process('picard-pcrmetrics')
    def test_collecttargetedpcrmetrics(self):
        with self.preparation_stage():
            bam_input = {
                'src': join('large', '56GSID_10k_mate1_RG.bam'),
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            bam = self.run_process('upload-bam', bam_input)
            master_file = self.prepare_amplicon_master_file()

            inputs = {
                'src': 'hs_b37_chr2_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            genome = self.run_process('upload-genome', inputs)

        inputs = {
            'alignment': bam.id,
            'master_file': master_file.id,
            'genome': genome.id
        }

        pcrmetrics = self.run_process('picard-pcrmetrics', inputs)
        self.assertFile(pcrmetrics, 'target_coverage', 'picard.perTargetCov.txt')

    @skipUnlessLargeFiles('56GSID_10k.realigned.bqsrCal.bam')
    @tag_process('vc-gatk-hc', 'vc-gatk4-hc')
    def test_gatk_haplotypecaller(self):
        with self.preparation_stage():
            alignment = self.run_process(
                'upload-bam', {
                    'src': join('large', '56GSID_10k.realigned.bqsrCal.bam'),
                    'species': 'Homo sapiens',
                    'build': 'b37'
                }
            )

            inputs = {
                'src': 'hs_b37_chr2_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            genome = self.run_process('upload-genome', inputs)

            master_file = self.prepare_amplicon_master_file()
            bed_file = self.run_process('upload-bed', {
                'src': './haplotypecaller/input/56G_masterfile_test_merged_targets_5col.bed',
                'species': 'Homo sapiens',
                'build': 'b37'
            })

            dbsnp_input = {
                'src': 'dbsnp_138.b37.chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            dbsnp = self.run_process('upload-variants-vcf', dbsnp_input)

        gatk3_vars = self.run_process('vc-gatk-hc', {
            'alignment': alignment.id,
            'intervals': master_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        })
        self.assertFile(
            gatk3_vars,
            'vcf',
            '56GSID_10k.gatkHC.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )
        self.assertFields(gatk3_vars, 'build', 'b37')
        self.assertFields(gatk3_vars, 'species', 'Homo sapiens')

        gatk4_vars = self.run_process('vc-gatk4-hc', {
            'alignment': alignment.id,
            'intervals': master_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        })
        self.assertFile(
            gatk4_vars,
            'vcf',
            '56GSID_10k.gatkHC4.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )
        self.assertFields(gatk4_vars, 'build', 'b37')
        self.assertFields(gatk4_vars, 'species', 'Homo sapiens')

        # This chunk checks that user is specifying only master file or bed
        # file, but not both. Once we remove the dependence on master file,
        # this test will be obsolete.
        gatk3_incol = self.run_process('vc-gatk-hc', {
            'alignment': alignment.id,
            'intervals': master_file.id,
            'intervals_bed': bed_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        }, Data.STATUS_ERROR)
        self.assertEqual(gatk3_incol.process_error[0],
                         'You have specified intervals and intervals_bed, whereas only one is permitted.')

        gatk4_incol = self.run_process('vc-gatk4-hc', {
            'alignment': alignment.id,
            'intervals': master_file.id,
            'intervals_bed': bed_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        }, Data.STATUS_ERROR)
        self.assertEqual(gatk4_incol.process_error[0],
                         'You have specified intervals and intervals_bed, whereas only one is permitted.')

        gatk3_bed = self.run_process('vc-gatk-hc', {
            'alignment': alignment.id,
            'intervals_bed': bed_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        })

        gatk4_bed = self.run_process('vc-gatk4-hc', {
            'alignment': alignment.id,
            'intervals_bed': bed_file.id,
            'genome': genome.id,
            'dbsnp': dbsnp.id
        })

        self.assertFile(
            gatk4_bed,
            'vcf',
            '56GSID_10k.gatkHC4.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )

        self.assertFile(
            gatk3_bed,
            'vcf',
            '56GSID_10k.gatkHC.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )

    @skipUnlessLargeFiles('56GSID_10k.realigned.bqsrCal.bam')
    @tag_process('lofreq')
    def test_lofreq(self):
        with self.preparation_stage():
            alignment = self.run_process(
                'upload-bam', {
                    'src': join('large', '56GSID_10k.realigned.bqsrCal.bam'),
                    'species': 'Homo sapiens',
                    'build': 'b37'
                }
            )

            inputs = {
                'src': 'hs_b37_chr2_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            genome = self.run_process('upload-genome', inputs)

            master_file = self.prepare_amplicon_master_file()

        inputs = {
            'alignment': alignment.id,
            'intervals': master_file.id,
            'genome': genome.id,
        }

        lofreq_vars = self.run_process('lofreq', inputs)
        self.assertFile(
            lofreq_vars,
            'vcf',
            '56GSID_10k.lf.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )
        self.assertFields(lofreq_vars, 'build', 'b37')
        self.assertFields(lofreq_vars, 'species', 'Homo sapiens')

    @tag_process('snpeff')
    def test_snpeff(self):
        with self.preparation_stage():
            variants_lf = self.run_process('upload-variants-vcf', {
                'src': '56GSID_10k.lf.vcf',
                'species': 'Homo sapiens',
                'build': 'b37',
            })
            variants_gatk = self.run_process('upload-variants-vcf', {
                'src': '56GSID_10k0.gatkHC.vcf',
                'species': 'Homo sapiens',
                'build': 'b37',
            })
            dbsnp = self.run_process('upload-variants-vcf', {
                'src': 'dbsnp_138.b37.chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37',
            })

        final_var_lf = self.run_process('snpeff', {
            'variants': variants_lf.id,
            'known_vars_annot': [dbsnp.id],
            'var_source': 'lofreq',
        })
        self.assertFile(final_var_lf, 'annotation', '56GSID.lf.finalvars.txt')
        self.assertRegex(final_var_lf.process_warning[0], r'Inconsistency for entry .*')

        final_var_gatk = self.run_process('snpeff', {
            'variants': variants_gatk.id,
            'known_vars_annot': [dbsnp.id],
            'var_source': 'gatk_hc',
        })
        self.assertFile(final_var_gatk, 'annotation', '56GSID.gatk.finalvars.txt')
