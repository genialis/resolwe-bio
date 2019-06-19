# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class ReadsFilteringProcessorTestCase(BioProcessTestCase):

    @tag_process('trimmomatic-single')
    def test_trimmomatic_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads()
            adapters = self.run_process('upload-fasta-nucl', {'src': 'bbduk_adapters.fasta'})

        inputs = {
            'reads': reads.pk,
            'illuminaclip': {
                'adapters': adapters.pk,
                'seed_mismatches': 2,
                'simple_clip_threshold': 10,
            },
            'maxinfo': {
                'target_length': 10,
                'strictness': 0.6,
            },
            'slidingwindow': {
                'window_size': 4,
                'required_quality': 15,
            },
            'trim_bases': {
                'leading': 20,
                'trailing': 20,
                'crop': 40,
                'headcrop': 3,
            },
            'reads_filtering': {
                'minlen': 22,
                'average_quality': 10,
            }}
        filtered_reads = self.run_processor('trimmomatic-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_trimmomatic_single.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc']}])

    @tag_process('trimmomatic-paired')
    def test_trimmomatic_paired(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['rRNA_forw.fastq.gz'],
                'src2': ['rRNA_rew.fastq.gz']}
            reads = self.run_processor('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk,
                  'trim_bases': {'trailing': 3}}

        filtered_reads = self.run_processor('trimmomatic-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_trimmomatic_paired_fw.fastq.gz'],
                         compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_trimmomatic_paired_rw.fastq.gz'],
                         compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA_forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA_forw_fastqc']}])
        del filtered_reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc']}])

    @tag_process('cutadapt-single')
    def test_cutadapt_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['cutadapt single.fastq.gz', 'cutadapt_single1.fastq.gz'])

            primers_up = self.prepare_adapters('5_prime_adapter.fasta.gz')
            primers_down = self.prepare_adapters('3_prime_adapter.fasta.gz')

        inputs = {
            'reads': reads.id,
            'polya_tail': 5,
            'down_primers_seq': ['AGCACCT'],
            'up_primers_seq': ['AGCTAAA'],
            'minlen': 10
        }

        cutadapt_single = self.run_process('cutadapt-single', inputs)

        self.assertFiles(cutadapt_single, 'fastq', ['cutadapt_single_trimmed.fastq.gz'],
                         compression='gzip')

        inputs = {
            'reads': reads.id,
            'polya_tail': 5,
            'down_primers_file': primers_down.id,
            'up_primers_file': primers_up.id,
            'minlen': 10
        }

        cutadapt_single = self.run_process('cutadapt-single', inputs)

        self.assertFiles(cutadapt_single, 'fastq', ['cutadapt_single_trimmed.fastq.gz'],
                         compression='gzip')

    @tag_process('cutadapt-paired')
    def test_cutadapt_paired(self):
        with self.preparation_stage():
            reads = self.prepare_paired_reads(mate1=['cutadapt mate1.fastq.gz'],
                                              mate2=['cutadapt mate2.fastq.gz'])

            primers_up = self.prepare_adapters('5_prime_adapter.fasta.gz')
            primers_down = self.prepare_adapters('3_prime_adapter.fasta.gz')

        inputs = {
            'reads': reads.id,
            'adapters': {
                'mate1_3prime_seq': ['AGCACCT'],
                'mate2_3prime_seq': ['AGCACCT'],
                'mate1_5prime_seq': ['AGCTAAA'],
                'mate2_5prime_seq': ['AGCTAAA'],
            },
            'filtering': {
                'minlen': 10,
            },
        }

        cutadapt_paired = self.run_process('cutadapt-paired', inputs)

        self.assertFiles(cutadapt_paired, 'fastq', ['cutadapt_paired_forward_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq2', ['cutadapt_paired_reverse_trimmed.fastq.gz'],
                         compression='gzip')

        inputs = {
            'reads': reads.id,
            'adapters': {
                'mate1_3prime_file': primers_down.id,
                'mate2_3prime_file': primers_down.id,
                'mate1_5prime_file': primers_up.id,
                'mate2_5prime_file': primers_up.id,
            },
            'filtering': {
                'minlen': 10,
            }
        }

        cutadapt_paired = self.run_process('cutadapt-paired', inputs)

        self.assertFiles(cutadapt_paired, 'fastq', ['cutadapt_paired_forward_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq2', ['cutadapt_paired_reverse_trimmed.fastq.gz'],
                         compression='gzip')

    @tag_process('cutadapt-custom-single', 'cutadapt-custom-paired')
    def test_cutadapt_custom(self):
        with self.preparation_stage():
            reads_single = self.prepare_reads(['cutadapt single.fastq.gz', 'cutadapt_single1.fastq.gz'])
            reads_paired = self.prepare_paired_reads(
                mate1=['cutadapt mate1.fastq.gz'],
                mate2=['cutadapt mate2.fastq.gz']
            )

        inputs_single = {'reads': reads_single.id}
        inputs_paired = {'reads': reads_paired.id}

        cutadapt_single = self.run_process('cutadapt-custom-single', inputs_single)
        cutadapt_paired = self.run_process('cutadapt-custom-paired', inputs_paired)

        self.assertFiles(cutadapt_single, 'fastq', ['cutadapt_custom_single_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq', ['cutadapt_custom_paired_forward_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq2', ['cutadapt_custom_paired_reverse_trimmed.fastq.gz'],
                         compression='gzip')

    @tag_process('bbduk-single')
    def test_bbduk_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['bbduk test reads.fastq.gz', 'rRNA forw.fastq.gz'])

        inputs = {
            'reads': reads.id,
        }
        filtered_reads = self.run_process('bbduk-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['bbduk_reads.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        report = {
            'file': 'fastqc/bbduk test reads_preprocessed_fastqc/fastqc_report.html',
            'refs': [
                'fastqc/bbduk test reads_preprocessed_fastqc',
            ],
        }
        self.assertFields(filtered_reads, "fastqc_url", [report])

    @tag_process('bbduk-paired')
    def test_bbduk_paired(self):
        with self.preparation_stage():
            reads_paired = self.prepare_paired_reads(['rRNA forw.fastq.gz'], ['rRNA_rew.fastq.gz'])

        inputs = {
            'reads': reads_paired.id,
        }
        filtered_reads = self.run_process('bbduk-paired', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['bbduk_fw_reads.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['bbduk_rv_reads.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        report = {
            'file': 'fastqc/rRNA forw_preprocessed_fastqc/fastqc_report.html',
            'refs': [
                'fastqc/rRNA forw_preprocessed_fastqc',
            ],
        }
        self.assertFields(filtered_reads, "fastqc_url", [report])
        del filtered_reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        report2 = {
            'file': 'fastqc/rRNA_rew_preprocessed_fastqc/fastqc_report.html',
            'refs': [
                'fastqc/rRNA_rew_preprocessed_fastqc',
            ],
        }
        self.assertFields(filtered_reads, "fastqc_url2", [report2])

    @tag_process('bamclipper')
    def test_bamclipper(self):
        species = 'Homo sapiens'
        build = 'fake_genome_RSEM'

        with self.preparation_stage():
            bam = self.prepare_bam(
                fn='./bamclipper/input/STK11.bam',
                species=species,
                build=build
            )

            inputs_bedpe = {'src': './bamclipper/input/STK11.bedpe',
                            'species': species, 'build': build}
            bedpe = self.run_process('upload-bedpe', inputs_bedpe)

        inputs_bamclipper = {'alignment': bam.id, 'bedpe': bedpe.id}
        clipped = self.run_process('bamclipper', inputs_bamclipper)

        self.assertFile(clipped, 'stats', './bamclipper/output/STK11.primerclipped.bam_stats.txt')
        self.assertFile(clipped, 'bigwig', './bamclipper/output/STK11.primerclipped.bw')
        self.assertFields(clipped, 'species', species)
        self.assertFields(clipped, 'build', build)

    @tag_process('markduplicates')
    def test_markduplicates(self):
        species = 'Homo sapiens'
        build = 'custombuild'
        primerclipped = './bamclipper/output/STK11.primerclipped.bam'

        with self.preparation_stage():
            bam = self.prepare_bam(
                fn=primerclipped,
                species=species,
                build=build)

        # Teste if skipped. Input bam should always equal output bam.
        md_inputs = {'bam': bam.id, 'skip': True}
        skipped_md = self.run_process('markduplicates', md_inputs)
        self.assertFile(skipped_md, 'bam', primerclipped)

        # Test that removal of duplicates works.
        md_inputs = {'bam': bam.id, 'remove_duplicates': True}
        removed_md = self.run_process('markduplicates', md_inputs)

        def filter_startedon(line):
            return line.startswith(b'# Started on:') or line.startswith(b'# MarkDuplicates')

        self.assertFileExists(removed_md, 'bam')
        self.assertFileExists(removed_md, 'bai')
        self.assertFile(removed_md, 'stats', './markduplicate/output/STK11.primerclipped.markduplicates.bam_stats.txt')
        self.assertFile(removed_md, 'bigwig', './markduplicate/output/STK11.primerclipped.markduplicates.bw')
        self.assertFile(removed_md, 'metrics_file', './markduplicate/output/STK11.primerclipped_metrics.txt',
                        file_filter=filter_startedon)
        self.assertFields(removed_md, 'species', species)
        self.assertFields(removed_md, 'build', build)

    @tag_process('bqsr')
    def test_bqsr(self):
        species = 'Homo sapiens'
        build = 'custom_build'
        with self.preparation_stage():
            input_genome = {
                # Based on b37 genome, chromosome 19 has been cut from beginning up to position 1207173.
                # This includes an exon of STK11. Cutting from the start of the chromosome was done so that
                # there is no need to shift any subsequent bed and vcf files.
                'src': './bqsr/input/hs_b37_chr19_upto_stk11.fasta.gz',
                'species': species,
                'build': build
            }
            input_bam = {
                'src': './markduplicate/output/STK11.primerclipped.markduplicates.bam',
                'species': species,
                'build': build
            }

            ks_dbsnp = []
            # dbSNP vcf has been cut to only specific region using VCFtools - 0.1.15. dbSNP database version is 138.
            # vcftools --gzvcf dbsnp_138.b37.vcf.gz \
            # --chr 19 \
            # --from-bp 1206985 \
            # --to-bp 1207173 \
            # --recode --recode-INFO-all \
            # --stdout > dbsnp_STK11.vcf
            # Header has been fixed by removing all but contig 19. Length has also been changed
            # to 1207200.
            for i in ['./bqsr/input/dbsnp_STK11.vcf.gz']:  # add more files if needed
                ks_dbsnp.append(
                    self.run_process('upload-variants-vcf', {'src': i, 'species': species, 'build': build})
                )

            intervals = self.run_process('upload-bed', {
                'src': './bqsr/input/STK11.bed',
                'species': species,
                'build': build})

            bam = self.run_process('upload-bam', input_bam)
            reference = self.run_process('upload-genome', input_genome)

            bqsr_inputs = {
                'bam': bam.id,
                'reference': reference.id,
                'known_sites': [i.id for i in ks_dbsnp],
                'intervals': intervals.id
            }
            bqsr = self.run_process('bqsr', bqsr_inputs)

            self.assertFileExists(bqsr, 'bam')
            self.assertFileExists(bqsr, 'bai')
            self.assertFile(bqsr, 'stats', './bqsr/output/STK11.primerclipped.markduplicates.bam_stats.txt')
            self.assertFile(bqsr, 'bigwig', './bqsr/output/STK11.primerclipped.markduplicates.bw')
            self.assertFile(bqsr, 'recal_table',
                            './bqsr/output/STK11.primerclipped.markduplicates_recalibration.table')
            self.assertFields(bqsr, 'species', species)
            self.assertFields(bqsr, 'build', build)

            # Check if read groups has successfully been added.
            bqsr_inputs['read_group'] = '-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1'
            bqsr_rg = self.run_process('bqsr', bqsr_inputs)

            self.assertFileExists(bqsr_rg, 'bam')
            self.assertFileExists(bqsr_rg, 'bai')
