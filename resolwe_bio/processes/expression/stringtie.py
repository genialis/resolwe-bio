"""StringTie."""
import os

from plumbum import TEE

import pandas as pd

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    GroupField,
    JsonField,
    Process,
    SchedulingClass,
    StringField,
)


def reformat_stringtie_output(infile, outfile):
    """Reformat StringTie gene abundance file."""
    exp = pd.read_csv(
        infile,
        sep='\t',
        usecols=['Gene ID', 'TPM'],
        index_col='Gene ID',
        dtype={
            'Gene ID': str,
            'TPM': float,
        },
        squeeze=True,
    )
    return exp.to_csv(outfile, index_label='Gene', header=['Expression'], sep='\t', compression='gzip')


class StringTie(Process):
    """Quantify transcript/gene abundance using StringTie."""

    slug = 'stringtie'
    name = "StringTie"
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/rnaseq:4.8.0',
            },
        },
        'resources': {
            'cores': 8,
            'memory': 16384,
            'network': True,
        },
    }
    data_name = "Stringtie ({{alignment|sample_name}})"
    version = '1.0.0'
    process_type = 'data:expression:stringtie'
    category = 'Quantify'
    entity = {
        'type': 'sample',
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        alignment = DataField('alignment:bam', label="Alignment")
        annotation = DataField('annotation:gtf', label="GTF annotation")

        class Options:
            """Options."""

            stranded = StringField(
                label="Assay type",
                default='non_specific',
                choices=[
                    ('non_specific', 'Strand non-specific'),
                    ('forward', 'Strand-specific forward'),
                    ('reverse', 'Strand-specific reverse'),
                ],
            )

        options = GroupField(Options, label="Options")

    class Output:
        """Output fields."""

        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        rc = FileField(label="Read counts", required=False)
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        abundance = FileField(label='StringTie abundance file')
        transcripts = FileField(label='StringTie transcripts GTF file')
        ctab = FileField(label='StringTie transcripts ctab file')
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if inputs.alignment.species != inputs.annotation.species:
            self.error('Sample species ({}) must match GTF annotation '
                       'file species ({})'.format(inputs.alignment.species, inputs.annotation.species))

        if inputs.alignment.build != inputs.annotation.build:
            self.error('Sample alignment file build ({}) must match GTF annotation '
                       'file build ({})'.format(inputs.alignment.build, inputs.annotation.build))

        name = os.path.basename(inputs.alignment.bam.path).strip('.bam')
        abundance_file = '{}_gene_abundance.txt'.format(name)
        transcripts_file = '{}_transcripts.gtf'.format(name)
        tpm_outfile = '{}_tpm.txt.gz'.format(name)

        args = [
            '-p', self.requirements.resources.cores,
            '-e',
            '-B',
            '-G', inputs.annotation.annot.path,
            '-A', abundance_file,
            '-o', transcripts_file,
        ]

        # If applicable, set the strandedness flag
        if inputs.options.stranded == 'forward':
            args.append('--fr')
        elif inputs.options.stranded == 'reverse':
            args.append('--rf')

        return_code, _, _ = Cmd['stringtie'][inputs.alignment.bam.path][args] & TEE(retcode=None)
        if return_code:
            self.error('StringTie analysis failed.')

        # Reformat StringTie output to Genialis expressions format
        reformat_stringtie_output(abundance_file, tpm_outfile)

        # Rename the ctab output file used of DE analysis
        os.rename('t_data.ctab', name + '_t_data.ctab')

        # Save the abundance estimates to JSON storage
        Cmd['expression2storage.py']('--output', 'json.txt', tpm_outfile)

        # Prepare expression set file with feature_id -> gene_id mappings
        exp_set_args = [
            '--expressions', tpm_outfile,
            '--source_db', inputs.annotation.source,
            '--species', inputs.annotation.species,
            '--output_name', name + '_expressions',
            '--expressions_type', 'TPM',
        ]
        return_code, _, _ = Cmd['create_expression_set.py'][exp_set_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing the expression set file.")

        # Save all the outputs
        outputs.exp = tpm_outfile
        outputs.exp_json = 'json.txt'
        outputs.exp_set = name + '_expressions.txt.gz'
        outputs.exp_set_json = name + '_expressions.json'
        outputs.abundance = abundance_file
        outputs.transcripts = transcripts_file
        outputs.ctab = name + '_t_data.ctab'
        outputs.exp_type = 'TPM'
        outputs.feature_type = 'gene'
        outputs.source = inputs.annotation.source
        outputs.species = inputs.alignment.species
        outputs.build = inputs.alignment.build
