"""Prepare MultiQC report."""
import json
import os
import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DirField,
    FileHtmlField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    StringField
)


def create_symlink(src, dst):
    """Create a symbolic link."""
    return Cmd['ln']('-s', '--backup=numbered', src, dst)


def create_summary_table(samples, species, build):
    """Prepare sample summary MultiQC table."""
    sample_summary_json = {
        'id': 'sample_info',
        'section_name': 'Sample Info',
        'plot_type': 'table',
        'file_format': 'json',
        'data': {}
    }

    for sample_name, sample_species, sample_build in zip(samples, species, build):
        if sample_build not in ['rRNA', 'globin']:
            sample_summary_json['data'][sample_name] = {
                'Species': sample_species,
                'Genome Build': sample_build,
            }

    with open('sample_data_mqc.json', 'w') as out_file:
        json.dump(sample_summary_json, out_file)


def parse_chip_qc_report(report):
    """Parse ChiP-seq QC report file."""
    df = pd.read_csv(report, sep='\t')
    df.fillna('', inplace=True)
    return df.to_dict(orient='records')[0]


def create_prepeak_table(sample_names, reports):
    """Prepare ChIP-seq pre-peak MultiQC table."""
    prepeak_qc_json = {
        'id': 'chip_seq_prepeak_qc',
        'section_name': 'ChIP-seq pre-peak QC',
        'plot_type': 'table',
        'file_format': 'json',
        'data': {}
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_chip_qc_report(report)
        prepeak_qc_json['data'][sample_name] = report_data

    with open('chipseq_prepeak_qc_mqc.json', 'w') as out_file:
        json.dump(prepeak_qc_json, out_file)


def create_postpeak_table(sample_names, reports):
    """Prepare ChIP-seq pre-peak MultiQC table."""
    postpeak_qc_json = {
        'id': 'chip_seq_postpeak_qc',
        'section_name': 'ChIP-seq post-peak QC',
        'plot_type': 'table',
        'file_format': 'json',
        'data': {}
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_chip_qc_report(report)
        postpeak_qc_json['data'][sample_name] = report_data

    with open('chipseq_postpeak_qc_mqc.json', 'w') as out_file:
        json.dump(postpeak_qc_json, out_file)


def create_lib_strand_table(samples, reports):
    """Prepare library strandedness MultiQC table."""
    strand_codes = {
        'IU': 'Strand non-specific (paired-end; -fr-unstranded)',
        'U': 'Strand non-specific (single-end; -fr-unstranded)',
        'ISF': 'Strand-specific forward (paired-end; -fr-secondstrand)',
        'OSF': 'Strand-specific forward (paired-end; outward facing reads)',
        'SF': 'Strand-specific forward (single-end; -fr-secondstrand)',
        'ISR': 'Strand-specific reverse (paired-end; -fr-firststrand)',
        'OSR': 'Strand-specific reverse (paired-end; outward facing reads)',
        'SR': 'Strand-specific reverse (single-end; -fr-firststrand)',
    }

    lib_strand_json = {
        'id': 'lib_strandedness',
        'section_name': 'Library Strandedness',
        'plot_type': 'table',
        'file_format': 'json',
        'data': {}
    }

    for sample_name, report in zip(samples, reports):
        with open(report) as infile:
            data = json.load(infile)
            if 'expected_format' in data:
                strandedness = data['expected_format']
            else:
                self.error("Cannot parse library type information file.")

        lib_strand_json['data'][sample_name] = {
            'Strandedness code': strandedness,
            'Description': strand_codes[strandedness],
        }

    with open('lib_strandedness_mqc.json', 'w') as out_file:
        json.dump(lib_strand_json, out_file)


def process_strand_report_file(data, lib_type_samples, lib_type_reports):
    """Process Strandedness report file if it exists as Data output file."""
    try:
        if os.path.isfile(data.strandedness_report.path):
            lib_type_samples.append(data.entity_name)
            lib_type_reports.append(data.strandedness_report.path)
    except AttributeError:
        pass


class MultiQC(Process):
    """Aggregate results from bioinformatics analyses across many samples into a single report.

    [MultiQC](http://www.multiqc.info) searches a given directory for analysis logs and compiles a
    HTML report. It's a general purpose tool, perfect for summarising the output from numerous
    bioinformatics tools.
    """

    slug = 'multiqc'
    process_type = 'data:multiqc'
    name = 'MultiQC'
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/common:1.3.1'
            },
        },
        'resources': {
            'cores': 1,
            'memory': 8192,
        },
    }
    entity = {
        'type': 'sample',
    }
    category = 'Other'
    data_name = 'MultiQC report'
    version = '1.5.0'

    class Input:
        """Input fields to process MultiQC."""

        data = ListField(
            DataField(
                data_type='',
                description='Select multiple data objects for which the MultiQC report is to be '
                            'generated.'
            ),
            label='Input data'
        )

        class Advanced:
            """Options."""

            dirs = BooleanField(
                label='--dirs',
                default=True,
                description='Prepend directory to sample names.',
            )

            dirs_depth = IntegerField(
                label='--dirs-depth',
                default=-1,
                description='Prepend a specified number of directories to sample names. Enter a '
                            'negative number (default) to take from start of path.',
            )

            fullnames = BooleanField(
                label='--fullnames',
                default=False,
                description='Disable the sample name cleaning (leave as full file name).',
            )

            config = BooleanField(
                label='Use configuration file',
                default=True,
                description='Use Genialis configuration file for MultiQC report.',
            )

            cl_config = StringField(
                label='--cl-config',
                required=False,
                description='Enter text with command-line configuration options to override the '
                            'defaults (e.g. custom_logo_url: https://www.genialis.com).',
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields."""

        report = FileHtmlField(label='MultiQC report')
        report_data = DirField(label='Report data')

    def run(self, inputs, outputs):
        """Run the analysis."""
        samples = []
        species = []
        build = []
        lib_type_samples = []
        lib_type_reports = []
        chip_seq_samples = []
        chip_seq_prepeak_reports = []
        chip_seq_postpeak_samples = []
        chip_seq_postpeak_reports = []
        unsupported_data = []

        for d in inputs.data:
            sample_dir = d.entity_name

            if sample_dir:
                os.makedirs(sample_dir, exist_ok=True)

            try:
                if d.entity_name and d.species and d.build:
                    samples.append(d.entity_name)
                    species.append(d.species)
                    build.append(d.build)
            except AttributeError:
                pass

            if d.type.startswith('data:reads:fastq:single'):
                for fq_report in d.fastqc_archive:
                    name = os.path.basename(fq_report.path)
                    create_symlink(fq_report.path, os.path.join(sample_dir, name))

            elif d.type.startswith('data:reads:fastq:paired'):
                for fq_report in d.fastqc_archive + d.fastqc_archive2:
                    name = os.path.basename(fq_report.path)
                    create_symlink(fq_report.path, os.path.join(sample_dir, name))

            elif d.type == 'data:alignment:bam:star:':
                stats_file = os.path.basename(d.stats.path)
                assert stats_file.endswith('_stats.txt')
                bam_name = stats_file[:-10]

                if d.build == 'rRNA':
                    rrna_report = f'{bam_name}.rRNA.Log.final.out'
                    create_symlink(d.stats.path, os.path.join(sample_dir, rrna_report))
                elif d.build == 'globin':
                    globin_report = f'{bam_name}.globin.Log.final.out'
                    create_symlink(d.stats.path, os.path.join(sample_dir, globin_report))
                else:
                    report = f'{bam_name}.Log.final.out'
                    create_symlink(d.stats.path, os.path.join(sample_dir, report))

            elif d.type.startswith('data:alignment:bam'):
                name = os.path.basename(d.stats.path)
                create_symlink(d.stats.path, os.path.join(sample_dir, name))

            elif d.type == 'data:expression:featurecounts:':
                name = os.path.basename(d.counts_summary.path)
                create_symlink(d.counts_summary.path, os.path.join(sample_dir, name))
                # Strandedness report exists only if auto detection was enabled
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.type == 'data:chipseq:callpeak:macs2:':
                name = os.path.basename(d.called_peaks.path)
                create_symlink(d.called_peaks.path, os.path.join(sample_dir, name))
                chip_seq_samples.append(d.entity_name)
                chip_seq_prepeak_reports.append(d.case_prepeak_qc.path)
                chip_seq_postpeak_samples.append(d.entity_name)
                chip_seq_postpeak_reports.append(d.chip_qc.path)
                # MACS2 analysis can be run without the background sample,
                # thus the associated report might not exits
                try:
                    if os.path.isfile(d.control_prepeak_qc.path):
                        chip_seq_samples.append(f'Background of {d.entity_name}')
                        chip_seq_prepeak_reports.append(d.control_prepeak_qc.path)
                except AttributeError:
                    pass

            elif d.type == 'data:samtools:idxstats:':
                name = os.path.basename(d.report.path)
                create_symlink(d.report.path, os.path.join(sample_dir, name))

            elif d.type == 'data:qorts:qc:':
                qorts_path = os.path.join(sample_dir, 'QoRTs')
                os.mkdir(qorts_path)
                name = os.path.basename(d.summary.path)
                create_symlink(d.summary.path, os.path.join(qorts_path, name))

            elif d.type == 'data:expression:salmon:':
                create_symlink(d.salmon_output.path, sample_dir)
                # Strandedness report might not exist in legacy Salmon objects
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.type.startswith('data:alleyoop'):
                create_symlink(d.report.path, sample_dir)
                # Alleyoop summary may contain PCA plot data
                try:
                    if os.path.isfile(d.plot_data.path):
                        create_symlink(d.plot_data.path, sample_dir)
                except AttributeError:
                    pass

            elif d.type.startswith('data:picard'):
                name = os.path.basename(d.report.path)
                create_symlink(d.report.path, os.path.join(sample_dir, name))
            else:
                unsupported_data.append(d.name)

        if unsupported_data:
            ext = ', ...' if len(unsupported_data) > 5 else ''
            self.warning(f"The Input data {', '.join(unsupported_data[:5])}{ext} is not supported "
                         f"by the MultiQC analysis.")

        create_summary_table(samples, species, build)

        if lib_type_samples and lib_type_reports:
            create_lib_strand_table(lib_type_samples, lib_type_reports)

        if chip_seq_samples and chip_seq_prepeak_reports:
            create_prepeak_table(chip_seq_samples, chip_seq_prepeak_reports)

        if chip_seq_postpeak_samples and chip_seq_postpeak_reports:
            create_postpeak_table(chip_seq_postpeak_samples, chip_seq_postpeak_reports)

        args = [
            '-dd', inputs.advanced.dirs_depth,
        ]

        if inputs.advanced.dirs:
            args.append('-d')

        if inputs.advanced.fullnames:
            args.append('-s')

        if inputs.advanced.config:
            args.extend(['-c', '/opt/resolwebio/assets/multiqc_config.yml'])

        if inputs.advanced.cl_config:
            args.extend(['--cl-config ', inputs.advanced.cl_config])

        with Cmd.env(LC_ALL='C.UTF-8'):
            return_code, _, _ = Cmd['multiqc']['.'][args] & TEE(retcode=None)
            if return_code:
                self.error('MultiQC analysis failed.')

        if not os.path.isdir('multiqc_data') and not os.path.isfile('multiqc_report.html'):
            self.error('MultiQC finished without creating outputs.')

        outputs.report = 'multiqc_report.html'
        outputs.report_data = 'multiqc_data'
