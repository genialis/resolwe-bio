"""Salmon Quant."""
import json
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DirField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    JsonField,
    Process,
    SchedulingClass,
    StringField,
)


class SalmonQuant(Process):
    """Perform mapping-based estimation of transcript abundance from RNA-seq reads.

    Final abundance estimates are summarized to the gene-level using
    [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html).
    """

    slug = 'salmon-quant'
    name = "Salmon Quant"
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/rnaseq:4.5.0',
            },
        },
        'resources': {
            'cores': 8,
            'memory': 16384,
            'network': True,
        },
    }
    data_name = "{{ reads.fastq.0.file|default('?') }}"
    version = '1.0.0'
    process_type = 'data:expression:salmon'
    category = 'Quantify'
    entity = {
        'type': 'sample',
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        reads = DataField('reads:fastq', label="Input sample(s)")
        salmon_index = DataField('index:salmon', label="Salmon index")
        annotation = DataField('annotation:gtf', label="GTF annotation")

        advanced = BooleanField(
            label="Show advanced options",
            description="Inspect and modify parameters.",
            default=False,
        )

        class Options:
            """Options."""

            stranded = StringField(
                label="Assay type",
                default='A',
                choices=[
                    ('A', 'Detect automatically'),
                    ('U', 'Strand non-specific (U)'),
                    ('SF', 'Strand-specific forward (SF)'),
                    ('SR', 'Strand-specific reverse (SR)'),
                    ('IU', 'Strand non-specific (paired-end IU)'),
                    ('ISF', 'Strand-specific forward (paired-end ISF)'),
                    ('ISR', 'Strand-specific reverse (paired-end (ISR)'),
                ],
            )

            seq_bias = BooleanField(
                label="--seqBias",
                default=False,
                description="Perform sequence-specific bias correction.",
            )

            gc_bias = BooleanField(
                label="--gcBias",
                default=False,
                description="[beta for single-end reads] Perform fragment GC bias correction.",
            )

            discard_orphans_quasi = BooleanField(
                label="--discardOrphansQuasi",
                default=False,
                description="Discard orphan mappings in quasi-mapping mode. "
                            "If this flag is passed then only paired "
                            "mappings will be considered toward "
                            "quantification estimates. The default "
                            "behavior is to consider orphan mappings "
                            "if no valid paired mappings exist.",
            )

            no_length_correction = BooleanField(
                label="--noLengthCorrection",
                default=False,
                description="[Experimental] Entirely disables "
                            "length correction when estimating the "
                            "abundance of transcripts. The abundance "
                            "estimates are reported in CPM (counts per "
                            "million) unit. This option can be used "
                            "with protocols where one expects that "
                            "fragments derive from their underlying "
                            "targets without regard to that target's  "
                            "length (e.g. QuantSeq).",
            )

            validate_mappings = BooleanField(
                label="--validateMappings",
                default=True,
                description="Validate mappings using alignment-based "
                            "verification. If this flag is passed, "
                            "quasi-mappings will be validated to "
                            "ensure that they could give rise to a "
                            "reasonable alignment before they are "
                            "further used for quantification.",
            )

            consensus_slack = FloatField(
                label="--consensusSlack",
                required=False,
                hidden="!options.validate_mappings",
                description="The amount of slack allowed in the quasi-mapping "
                            "consensus mechanism.  Normally, a transcript must "
                            "cover all hits to be considered for mapping.  "
                            "If this is set to a fraction, X, greater than 0 "
                            "(and in [0,1)), then a transcript can fail "
                            "to cover up to (100 * X)% of the hits before it "
                            "is discounted as a mapping candidate. The default "
                            "value of this option is 0.2 if --validateMappings "
                            "is given and 0 otherwise",
            )

            min_score_fraction = FloatField(
                label="--minScoreFraction",
                default=0.65,
                hidden="!options.validate_mappings",
                description="The fraction of the optimal possible alignment "
                            "score that a mapping must achieve in order to be "
                            "considered valid - should be in (0,1]",
            )

            incompat_prior = FloatField(
                label="---incompatPrior",
                default=0,
                description="This option sets the prior probability "
                            "that an alignment that disagrees with "
                            "the specified library type (--libType) "
                            "results from the true fragment origin. "
                            "Setting this to 0 specifies that "
                            "alignments that disagree with the "
                            "library type should be impossible, "
                            "while setting it to 1 says that "
                            "alignments that disagree with the "
                            "library type are no less likely than "
                            "those that do.",
            )

            range_factorization_bins = IntegerField(
                label="--rangeFactorizationBins",
                default=0,
                description="Factorizes the likelihood used in "
                            "quantification by adopting a new notion "
                            "of equivalence classes based on the "
                            "conditional probabilities with which "
                            "fragments are generated from different "
                            "transcripts.  This is a more "
                            "fine-grained factorization than the "
                            "normal rich equivalence classes.  The "
                            "default value (0) corresponds to the "
                            "standard rich equivalence classes, and "
                            "larger values imply a more fine-grained "
                            "factorization.  If range factorization "
                            "is enabled, a common value to select "
                            "for this parameter is 4."
            )

            min_assigned_frag = IntegerField(
                label="--minAssignedFrags",
                default=10,
                description="The minimum number of fragments that "
                            "must be assigned to the transcriptome "
                            "for quantification to proceed."
            )

        options = GroupField(Options, label="Options", hidden="!advanced")

    class Output:
        """Output fields."""

        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        rc = FileField(label="Read counts", required=False)
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        quant = FileField(label="Salmon quant file")
        salmon_output = DirField(label='Salmon output')
        txdb = FileField(label="Transcript to gene mapping")
        strandedness = StringField(label='Strandedness code')
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if inputs.salmon_index.species != inputs.annotation.species:
            self.error('Salmon index file species ({}) must match GTF annotation '
                       'file species ({})'.format(inputs.salmon_index.species, inputs.annotation.species))

        if inputs.salmon_index.build != inputs.annotation.build:
            self.error('Salmon index file build ({}) must match GTF annotation '
                       'file build ({})'.format(inputs.salmon_index.build, inputs.annotation.build))

        if inputs.salmon_index.source != inputs.annotation.source:
            self.error('Salmon index file source ({}) must match GTF annotation '
                       'file source ({})'.format(inputs.salmon_index.source, inputs.annotation.source))

        if inputs.options.no_length_correction:
            abundance_unit = 'CPM'
            output_suffix = '_cpm.txt.gz'
        else:
            abundance_unit = 'TPM'
            output_suffix = '_tpm.txt.gz'

        args = [
            '-i', inputs.salmon_index.index.path,
            '-l', inputs.options.stranded,
            '--incompatPrior', inputs.options.incompat_prior,
            '--minAssignedFrags', inputs.options.min_assigned_frag,
            '--rangeFactorizationBins', inputs.options.range_factorization_bins,
            '-p', self.requirements.resources.cores,
            '-o', 'salmon_output',
        ]

        # Prepare .FASTQ file inputs based on the reads input type
        if inputs.reads.type.startswith('data:reads:fastq:single:'):
            args.extend(['-r'] + [lane.path for lane in inputs.reads.fastq])
        else:
            args.extend(['-1'] + [lane.path for lane in inputs.reads.fastq])
            args.extend(['-2'] + [lane.path for lane in inputs.reads.fastq2])

        # Prepare optional inputs
        if inputs.options.seq_bias and not inputs.options.no_length_correction:
            args.append('--seqBias')
        else:
            self.warning("Since bias correction relies on modifying effective lengths, "
                         "you cannot enable bias correction simultaneously with the "
                         "--noLengthCorrection option. Skipping --seqBias option.")

        if inputs.options.gc_bias and not inputs.options.no_length_correction:
            args.append('--gcBias')
        else:
            self.warning("Since bias correction relies on modifying effective lengths, "
                         "you cannot enable bias correction simultaneously with the "
                         "--noLengthCorrection option. Skipping --gcBias option.")

        if inputs.options.discard_orphans_quasi:
            args.append('--discardOrphansQuasi')

        if inputs.options.no_length_correction:
            args.append('--noLengthCorrection')

        if inputs.options.validate_mappings:
            args.append('--validateMappings')

        if inputs.options.validate_mappings and inputs.options.min_score_fraction > 0:
            args.extend(['--minScoreFraction', inputs.options.min_score_fraction])

        if inputs.options.consensus_slack is not None:
            args.extend(['--consensusSlack', inputs.options.consensus_slack])

        # Run Salmon Quant
        return_code, _, _ = Cmd['salmon']['quant'][args] & TEE(retcode=None)
        if return_code:
            self.error("Error while running Salmon Quant.")

        # Use tximport to produce gene-level TPM values
        reads_name = os.path.basename(inputs.reads.fastq[0].path).strip('.fastq.gz')
        annot_name = os.path.basename(inputs.annotation.annot.path).strip('.gtf')
        tx2gene = 'tx2gene_{}.txt'.format(annot_name)
        if os.path.exists('salmon_output/quant.sf'):
            tximport_args = ['salmon_output/quant.sf', inputs.annotation.annot.path, reads_name, tx2gene]
            return_code, _, _ = Cmd['tximport_summarize.R'][tximport_args] & TEE(retcode=None)
            if return_code:
                self.error("Error while running tximport.")
        else:
            self.error('Salmon Quant results file quant.sf does not exists.')

        # Zip the gene-level abundance estimates
        (Cmd['gzip']['-c', reads_name] > reads_name + output_suffix)()

        # Save the abundance estimates to JSON storage
        Cmd['expression2storage.py']('--output', 'json.txt', reads_name + output_suffix)

        # Prepare expression set file with feature_id -> gene_id mappings
        exp_set_args = [
            '--expressions', reads_name + output_suffix,
            '--source_db', inputs.salmon_index.source,
            '--species', inputs.salmon_index.species,
            '--output_name', reads_name + '_expressions',
            '--expressions_type', abundance_unit,
        ]
        return_code, _, _ = Cmd['create_expression_set.py'][exp_set_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing the expression set file.")

        Cmd['ln']['-s', 'salmon_output/quant.sf', reads_name + '.sf']()

        strandedness = json.load(open('salmon_output/lib_format_counts.json')).get('expected_format', '')

        # Save all the outputs
        outputs.salmon_output = 'salmon_output'
        outputs.quant = reads_name + '.sf'
        outputs.txdb = tx2gene
        outputs.exp = reads_name + output_suffix
        outputs.exp_json = 'json.txt'
        outputs.exp_set = reads_name + '_expressions.txt.gz'
        outputs.exp_set_json = reads_name + '_expressions.json'
        outputs.strandedness = strandedness
        outputs.exp_type = abundance_unit
        outputs.feature_type = 'gene'
        outputs.source = inputs.salmon_index.source
        outputs.species = inputs.salmon_index.species
        outputs.build = inputs.salmon_index.build
