"""Run Alleyoop collapse tool on Slamdunk results."""
import os

from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, Process, StringField


class AlleyoopCollapse(Process):
    """Run Alleyoop collapse tool on Slamdunk results."""

    slug = 'alleyoop-collapse'
    process_type = 'data:alleyoop:collapse'
    name = 'Alleyoop collapse'
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/slamdunk:1.0.0'
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
    category = 'Slamdunk'
    data_name = '{{ slamdunk|sample_name|default("?") }}'
    version = '1.0.0'

    class Input:
        """Input fields for SlamdunkAllPaired."""

        slamdunk = DataField('alignment:bam:slamdunk', label='Slamdunk results')

    class Output:
        """Output fields to process SlamdunkAllPaired."""

        tcount = FileField(label='Count report containing SLAMSeq statistics')
        species = StringField(label='Species')
        build = StringField(label='Build')

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.slamdunk.tcount.path)
        assert basename.endswith('.tsv')
        name = basename[:-4]

        args = [
            '-o', '.',
            '-t', self.requirements.resources.cores,
        ]

        return_code, _, _ = Cmd['alleyoop']['collapse'][args][inputs.slamdunk.tcount.path] & TEE(retcode=None)
        if return_code:
            self.error('Alleyoop collapse analysis failed.')

        collapsed_output = name + '_collapsed.txt'
        os.rename(name + '_collapsed.csv', collapsed_output)

        outputs.tcount = collapsed_output
        outputs.species = inputs.slamdunk.species
        outputs.build = inputs.slamdunk.build
