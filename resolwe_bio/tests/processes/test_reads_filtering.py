# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase, skipDockerFailure


class ReadsFilteringProcessorTestCase(BioProcessTestCase):

    def test_prinseq_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}

        filtered_reads = self.run_process('prinseq-lite-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc'],
                                                          'size': 314749}])

    def test_prinseq_paired(self):
        inputs = {
            'src1': ['rRNA forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('prinseq-lite-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_prinseq_paired_rw.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA forw_fastqc'],
                                                          'size': 347773}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 340697}])

    def test_fastqmcf_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('fastq-mcf-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc'],
                                                          'size': 305101}])

    def test_fastqmcf_paired(self):
        inputs = {
            'src1': ['rRNA forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_process('fastq-mcf-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_fastqmcf_paired_rw.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA forw_fastqc'],
                                                          'size': 347791}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 340715}])

    @skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
    def test_sortmerna_single(self):
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
    def test_sortmerna_paired(self):
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

    def test_trimmomatic_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk,
                  'trailing': 3,
                  'crop': 5}
        filtered_reads = self.run_processor('trimmomatic-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_trimmomatic_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/reads_fastqc'],
                                                          'size': 206718}])

    def test_trimmomatic_paired(self):
        inputs = {
            'src1': ['rRNA_forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_processor('upload-fastq-paired', inputs)
        inputs = {'reads': reads.pk,
                  'trailing': 3}
        filtered_reads = self.run_processor('trimmomatic-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_trimmomatic_paired_fw.fastq.gz'],
                         compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_trimmomatic_paired_rw.fastq.gz'],
                         compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA_forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA_forw_fastqc'],
                                                          'size': 352347}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 340745}])

    def test_hsqutils_trimm(self):
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

    def test_cutadapt_single(self):
        reads = self.prepare_reads(['cutadapt single.fastq.gz', 'cutadapt_single1.fastq.gz'])

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

        primers_up = self.prepare_adapters('5_prime_adapter.fasta.gz')
        primers_down = self.prepare_adapters('3_prime_adapter.fasta.gz')

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

    def test_cutadapt_paired(self):
        reads = self.prepare_paired_reads(mate1=['cutadapt forward1.fastq.gz', 'cutadapt_forward2.fastq.gz'],
                                          mate2=['cutadapt_reverse.fastq.gz'])

        inputs = {
            'reads': reads.id,
            'polya_tail': 5,
            'down_primers_seq_fwd': ['AGCACCT'],
            'down_primers_seq_rev': ['AGCACCT'],
            'up_primers_seq_fwd': ['AGCTAAA'],
            'up_primers_seq_rev': ['AGCTAAA'],
            'minlen': 10
        }

        cutadapt_paired = self.run_process('cutadapt-paired', inputs)

        self.assertFiles(cutadapt_paired, 'fastq', ['cutadapt_paired_forward_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq2', ['cutadapt_paired_reverse_trimmed.fastq.gz'],
                         compression='gzip')

        primers_up = self.prepare_adapters('5_prime_adapter.fasta.gz')
        primers_down = self.prepare_adapters('3_prime_adapter.fasta.gz')

        inputs = {
            'reads': reads.id,
            'polya_tail': 5,
            'down_primers_file_fwd': primers_down.id,
            'down_primers_file_rev': primers_down.id,
            'up_primers_file_fwd': primers_up.id,
            'up_primers_file_rev': primers_up.id,
            'minlen': 10
        }

        cutadapt_paired = self.run_process('cutadapt-paired', inputs)

        self.assertFiles(cutadapt_paired, 'fastq', ['cutadapt_paired_forward_trimmed.fastq.gz'],
                         compression='gzip')

        self.assertFiles(cutadapt_paired, 'fastq2', ['cutadapt_paired_reverse_trimmed.fastq.gz'],
                         compression='gzip')

    def test_cutadapt_custom(self):
        reads_single = self.prepare_reads(['cutadapt single.fastq.gz', 'cutadapt_single1.fastq.gz'])
        reads_paired = self.prepare_paired_reads(mate1=['cutadapt forward1.fastq.gz', 'cutadapt_forward2.fastq.gz'],
                                                 mate2=['cutadapt_reverse.fastq.gz'])

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

    def test_cutadapt_amplicon(self):
        inputs = {
            'src1': ['56GSID_1k_mate1.fastq.gz'],
            'src2': ['56GSID_1k_mate2.fastq.gz']}
        reads = self.run_processor('upload-fastq-paired', inputs)

        inputs = {'src': '5ptrim_new56Gprimers.fa.gz'}
        primers_1 = self.run_processor('upload-fasta-nucl', inputs)

        inputs = {'src': '3ptrim_new56Gprimers.fa.gz'}
        primers_2 = self.run_processor('upload-fasta-nucl', inputs)

        inputs = {
            'reads': reads.id,
            'up_primers': primers_1.id,
            'down_primers': primers_2.id
        }

        filtered_reads = self.run_processor('cutadapt-amplicon', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['cutadapt_trimmed_mate1.fastq.gz'],
                         compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['cutadapt_trimmed_mate2.fastq.gz'],
                         compression='gzip')

    def test_bbduk_single(self):
        inputs = {'src': ['bbduk test reads.fastq.gz']}
        reads = self.run_processor('upload-fastq-single', inputs)

        inputs = {'src': 'bbduk_adapters.fasta'}
        adapters = self.run_processor('upload-fasta-nucl', inputs)

        inputs = {
            'reads': reads.id,
            'adapters': {
                'reference': [adapters.id]},
            'trimming_par': {
                'kmask': 'lc'},
            'barcode_par': {
                'barcode_seq': ['GCGTTT', 'TATGNNN']}
        }

        filtered_reads = self.run_process('bbduk-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_bbduk_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/bbduk test reads_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/bbduk test reads_fastqc'],
                                                          'size': 303594}])

        inputs = {
            'reads': reads.id,
            'adapters': {
                'bbduk_adapters': 'truseq'},
            'trimming_par': {
                'kmask': 'lc'},
            'barcode_par': {
                'barcode_seq': ['GCGTTT', 'TATGNNN']}
        }
        filtered_reads_1 = self.run_process('bbduk-single', inputs)
        self.assertFiles(filtered_reads_1, 'fastq', ['filtered_reads_bbduk_single.fastq.gz'], compression='gzip')

    def test_bbduk_paired(self):
        inputs = {'src': 'bbduk_adapters.fasta'}
        adapters = self.run_processor('upload-fasta-nucl', inputs)

        reads_paired = self.prepare_paired_reads(mate1=['rRNA forw.fastq.gz'],
                                                 mate2=['rRNA_rew.fastq.gz'])
        inputs = {
            'reads': reads_paired.id,
            'adapters': {
                'reference': [adapters.id],
                'literal_reference': ['AAAAAAAAAAAAAAAAAAAAAA', 'CCCCCCC'],
            },
            'trimming_par': {
                'ktrim': 'r'
            }}

        filtered_reads = self.run_process('bbduk-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_bbduk_fw.fastq.gz'],
                         compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_bbduk_rw.fastq.gz'],
                         compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                          'refs': ['fastqc/rRNA forw_fastqc'],
                                                          'size': 255848}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                           'refs': ['fastqc/rRNA_rew_fastqc'],
                                                           'size': 244724}])
