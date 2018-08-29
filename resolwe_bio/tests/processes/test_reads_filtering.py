# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase, skipDockerFailure


class ReadsFilteringProcessorTestCase(BioProcessTestCase):

    @tag_process('prinseq-lite-single')
    def test_prinseq_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads()
            inputs = {'reads': reads.pk}

        filtered_reads = self.run_process('prinseq-lite-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_single.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc'],
                                                          'size': 314749}])

    @tag_process('prinseq-lite-paired')
    def test_prinseq_paired(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['rRNA forw.fastq.gz'],
                'src2': ['rRNA_rew.fastq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('prinseq-lite-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_prinseq_paired_rw.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA forw_fastqc'],
                                                          'size': 347773}])
        del filtered_reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 340697}])

    @tag_process('fastq-mcf-single')
    def test_fastqmcf_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads()

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('fastq-mcf-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_single.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc'],
                                                          'size': 305101}])

    @tag_process('fastq-mcf-paired')
    def test_fastqmcf_paired(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['rRNA forw.fastq.gz'],
                'src2': ['rRNA_rew.fastq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('fastq-mcf-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_fastqmcf_paired_rw.fastq.gz'], compression='gzip')
        del filtered_reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA forw_fastqc'],
                                                          'size': 347791}])
        del filtered_reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 340715}])

    @skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
    @tag_process('sortmerna-single')
    def test_sortmerna_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['rRNA forw.fastq.gz'])
            rrnadb_1 = self.run_process('upload-fasta-nucl', {'src': 'silva-arc-16s-id95.fasta.gz'})
            rrnadb_2 = self.run_process('upload-fasta-nucl', {'src': 'silva-arc-23s-id98.fasta.gz'})

        inputs = {
            'reads': reads.id,
            'database_selection': [rrnadb_1.id, rrnadb_2.id],
            'options': {'threads': 2, 'sam': True}}
        filtered_reads = self.run_process('sortmerna-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['reads_wo_rRNA_single.fastq.gz'], compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_single.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam', {'file': 'rRNA forw_rRNA.sam'})
        self.assertFields(filtered_reads, 'stats', {'file': 'stats.log'})
        self.assertFields(filtered_reads, "fastqc_url",
                          [{'file': 'fastqc/rRNA forw_filtered_fastqc/fastqc_report.html',
                            'refs': ['fastqc/rRNA forw_filtered_fastqc'],
                            'size': 345492}])

    @skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
    @tag_process('sortmerna-paired')
    def test_sortmerna_paired(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['rRNA forw.fastq.gz'],
                'src2': ['rRNA_rew.fastq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

            rrnadb_1 = self.run_process('upload-fasta-nucl', {'src': 'silva-arc-16s-id95.fasta.gz'})
            rrnadb_2 = self.run_process('upload-fasta-nucl', {'src': 'silva-arc-23s-id98.fasta.gz'})

        inputs = {
            'reads': reads.id,
            'database_selection': [rrnadb_1.id, rrnadb_2.id],
            'options': {'threads': 2, 'sam': True}}

        filtered_reads = self.run_process('sortmerna-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['reads_wo_rRNA_paired_forw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['reads_wo_rRNA_paired_rew.fastq.gz'], compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_paired.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam', {'file': 'rRNA forw_rRNA.sam'})
        self.assertFields(filtered_reads, 'stats', {'file': 'stats.log'})
        self.assertFields(filtered_reads, "fastqc_url",
                          [{'file': 'fastqc/rRNA forw_filtered_fastqc/fastqc_report.html',
                            'refs': ['fastqc/rRNA forw_filtered_fastqc'],
                            'size': 345492}])
        self.assertFields(filtered_reads, "fastqc_url2",
                          [{'file': 'fastqc/rRNA_rew_filtered_fastqc/fastqc_report.html',
                            'refs': ['fastqc/rRNA_rew_filtered_fastqc'],
                            'size': 347212}])

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

    @tag_process('hsqutils-trim')
    def test_hsqutils_trim(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['hsqutils_reads_mate1_paired_filtered.fastq.gz'],
                'src2': ['hsqutils_reads_mate2_paired_filtered.fastq.gz']}
            reads = self.run_processor('upload-fastq-paired', inputs)

            probe = self.run_processor('upload-file', {'src': 'hsqutils_probe_info.txt'})

        inputs = {'reads': reads.id,
                  'probe': probe.id}

        hsqutils_trimm = self.run_processor('hsqutils-trim', inputs)

        self.assertFiles(hsqutils_trimm, 'fastq', ['hsqutils_reads_trimmed_mate1.fastq.gz'],
                         compression='gzip')
        self.assertFiles(hsqutils_trimm, 'fastq2', ['hsqutils_reads_trimmed_mate2.fastq.gz'],
                         compression='gzip')

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
            reads = self.prepare_paired_reads(mate1=['cutadapt forward1.fastq.gz', 'cutadapt_forward2.fastq.gz'],
                                              mate2=['cutadapt_reverse.fastq.gz'])

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
                mate1=['cutadapt forward1.fastq.gz', 'cutadapt_forward2.fastq.gz'],
                mate2=['cutadapt_reverse.fastq.gz']
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
