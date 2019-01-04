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
