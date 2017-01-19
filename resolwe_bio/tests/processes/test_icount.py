# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase, skipDockerFailure


@skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
class ICountImport(BioProcessTestCase):

    def test_genome(self):
        genome = self.run_process('icount-genome', {
            'species': 'homo_sapiens',
            'release': 84,
            'chromosomes': ['MT'],
        })

        self.assertFile(genome, "fasta", "icount.genome.out.fasta.gz", compression='gzip')
        self.assertFields(genome, "number", 1)
        self.assertFile(genome, "fai", "icount.genome.out.fasta.gz.fai")

    def test_annotation(self):
        # TODO: when iCount updates, change to saccharomyces_cerevisiae
        ann = self.run_process('icount-annotation', {
            'species': 'homo_sapiens',
            'release': 84,
        })

        self.assertFileExists(ann, "gtf")
        self.assertFields(ann, "source", "ENSEMBL")

    def test_import_bed_icount(self):
        bed = self.run_process('upload-bed-icount', {'src': 'icount.group.in1.bed'})

        self.assertFileExists(bed, "BED")
        self.assertFile(bed, "BED", 'icount.group.in1.bed')


@skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
class ICountPreprocess(BioProcessTestCase):

    def test_segment(self):
        inputs = {'src': 'icount.segment.in.gtf', 'source': 'ENSEMBL'}
        annotation = self.run_process('upload-gtf', inputs)
        inputs = {'species': 'homo_sapiens', 'release': 84, 'chromosomes': ['Y']}
        genome = self.run_process('icount-genome', inputs)

        segment = self.run_process('icount-segment', {
            'annotation': annotation.pk,
            'genome': genome.pk,
        })

        # Since the order of GTF attributes is arbitrary, we cant check for file
        # equality. Just check that file exists.
        self.assertFileExists(segment, "gtf")
        self.assertFile(segment, "types_length", 'icount.segment.out.txt')
        self.assertFile(segment, "fai", 'icount.segment.out.fai')
        self.assertFields(segment, "source", 'ENSEMBL')

    def test_demultiplex(self):
        inputs = {'src': ['icount.demultiplex.in.fq']}
        raw_reads = self.run_process('upload-fastq-single', inputs)

        demux_reads = self.run_process('icount-demultiplex', {
            'reads': raw_reads.pk,
            'adapter': 'CCCCCCCCC',
            'barcodes': ['NNNGGTTNN', 'NNNTTGTNN'],
        })

        self.assertFiles(demux_reads, "fastq", compression='gzip', fn_list=[
            "icount.demultiplex.out.NNNGGTTNN.fastq.gz",
            "icount.demultiplex.out.NNNTTGTNN.fastq.gz",
            "icount.demultiplex.out.nomatch.fastq.gz",
        ])
        # Check also that new data object(s) were created from demultiplex and
        # that they have correct files:
        nomatch_reads = Data.objects.last()
        self.assertFiles(nomatch_reads, "fastq", ["icount.demultiplex.out.nomatch.fastq.gz"], compression='gzip')

    def test_xlsites(self):
        inputs = {'src': 'icount.xlsites.in.bam'}
        bam = self.run_process('upload-bam', inputs)

        xlsites = self.run_process('icount-xlsites', {'bam': bam.pk})

        self.assertFile(xlsites, "BED", "icount.xlsites.out.unique.bed.gz", compression='gzip')
        self.assertFile(xlsites, "BED_multi", "icount.xlsites.out.multi.bed.gz", compression='gzip')
        self.assertFile(xlsites, "skipped", "icount.xlsites.out.skipped.bam")


@skipDockerFailure("Skip until Docker image with iCount is supported on Travis.")
class ICountAnalyses(BioProcessTestCase):

    def test_annotate(self):
        inputs = {'src': 'icount.annotate.in.fasta'}
        genome = self.run_process('upload-fasta-nucl', inputs)
        inputs = {'src': 'icount.annotate.in.gtf', 'source': 'ENSEMBL'}
        ann_input = self.run_process('upload-gtf', inputs)
        inputs = {'annotation': ann_input.pk, 'genome': genome.pk}
        seg_input = self.run_process('icount-segment', inputs)

        inputs = {'src': 'icount.annotate.in.bed'}
        bed_input = self.run_process('upload-bed-icount', inputs)

        annotate = self.run_process('icount-annotate', {
            'segmentation': seg_input.pk,
            'sites': bed_input.pk,
            'subtype': 'biotype',
            'excluded_types': ['ncRNA'],
        })

        self.assertFile(annotate, "sites_annotated", "icount.annotate.out.bed")

    def test_clusters(self):
        inputs = {'src': 'icount.clusters.in.bed'}
        bed_input = self.run_process('upload-bed-icount', inputs)

        inputs = {'sites': bed_input.pk, 'dist': 20}
        clusters = self.run_process('icount-clusters', inputs)

        self.assertFile(clusters, "clusters", "icount.clusters.out.bed")

    def test_group(self):
        inputs = {'src': 'icount.group.in1.bed'}
        bed_input1 = self.run_process('upload-bed-icount', inputs)
        inputs = {'src': 'icount.group.in2.bed'}
        bed_input2 = self.run_process('upload-bed-icount', inputs)

        inputs = {'sites': [bed_input1.pk, bed_input2.pk]}
        grouped = self.run_process('icount-group', inputs)

        self.assertFile(grouped, "sites_grouped", "icount.group.out.bed")

    def test_peaks(self):
        inputs = {'src': 'icount.peaks.in.fa'}
        genome = self.run_process('upload-fasta-nucl', inputs)
        inputs = {'src': 'icount.peaks.in.gtf', 'source': 'ENSEMBL'}
        ann_input = self.run_process('upload-gtf', inputs)
        inputs = {'annotation': ann_input.pk, 'genome': genome.pk}
        seg_input = self.run_process('icount-segment', inputs)

        inputs = {'src': 'icount.peaks.in.bed'}
        bed_input = self.run_process('upload-bed-icount', inputs)

        peaks = self.run_process('icount-peaks', {
            'segmentation': seg_input.pk,
            'sites': bed_input.pk,
            'features': ['gene'],
            'merge_features': True,
            'report_progress': True,
        })

        self.assertFile(peaks, "peaks", "icount.peaks.out1.bed.gz", compression='gzip')
        self.assertFile(peaks, "scores", "icount.peaks.out2.tsv.gz", compression='gzip')

    def test_summary(self):
        inputs = {'src': 'icount.summary.in.fa'}
        genome = self.run_process('upload-fasta-nucl', inputs)
        inputs = {'src': 'icount.summary.in.gtf', 'source': 'ENSEMBL'}
        ann_input = self.run_process('upload-gtf', inputs)
        inputs = {'annotation': ann_input.pk, 'genome': genome.pk}
        seg_input = self.run_process('icount-segment', inputs)

        inputs = {'src': 'icount.summary.in.bed'}
        bed_input = self.run_process('upload-bed-icount', inputs)

        inputs = {'segmentation': seg_input.pk, 'sites': bed_input.pk}
        summary = self.run_process('icount-summary', inputs)

        self.assertFile(summary, "summary", "icount.summary.out.tsv")

    def test_workflow(self):
        inputs = {'src': 'icount.workflow.in.fasta.gz'}
        genome = self.run_process('upload-fasta-nucl', inputs)
        inputs = {'src': 'icount.workflow.in.gtf', 'source': 'ENSEMBL'}
        ann = self.run_process('upload-gtf', inputs)
        inputs = {'annotation': ann.pk, 'genome': genome.pk}
        segmentation = self.run_process('icount-segment', inputs)

        inputs = {'genome2': genome.pk, 'annotation': ann.pk}
        star_index = self.run_process('alignment-star-index', inputs)

        inputs = {'src': ['icount.workflow.in.fastq']}
        reads = self.run_process('upload-fastq-single', inputs)

        self.run_process('workflow-icount', {
            'reads': reads.pk,
            'index': star_index.pk,
            'segmentation': segmentation.pk,
        })

        peaks = Data.objects.last()
        self.assertFile(peaks, "scores", "icount.workflow.out.tsv.gz", compression='gzip')
