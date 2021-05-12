from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class MicroRNATestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("workflow-mirna")
    def test_mirna_workflow(self):
        # Prepare data for aligning the reads with bowtie2 and annotation file for featureCounts.
        with self.preparation_stage():
            inputs = {
                "src": "genome_rsem.fa.gz",
                "species": "Homo sapiens",
                "build": "fake_genome_RSEM",
            }
            ref_seq = self.run_process("upload-fasta-nucl", inputs)
            bowtie2_index = self.run_process("bowtie2-index", {"ref_seq": ref_seq.id})
            single_reads = self.prepare_reads(["reads rsem.fq.gz"])
            annotation = self.prepare_annotation(
                "annotation_rsem.gtf.gz",
                species="Homo sapiens",
                build="fake_genome_RSEM",
            )

        inputs = {
            "preprocessing": {
                "reads": single_reads.pk,
                "adapters": {"down_primers_seq": ["TAATGAACAATGCAAGTTTGA"]},
                "filtering": {"minlen": 15, "maxlen": 35, "error_rate": 0.2},
            },
            "alignment": {
                "genome": bowtie2_index.pk,
                "alignment_options": {
                    "mode": "--local",
                    "speed": "--very-sensitive",
                    "L": 8,
                    "rep_mode": "k",
                    "k_reports": 5,
                },
            },
            "quant_options": {
                "annotation": annotation.pk,
                "id_attribute": "gene_id",
                "feature_class": "exon",
                "normalization_type": "CPM",
                "count_multi_mapping_reads": True,
                "allow_multi_overlap": True,
            },
            "assay_type": "non_specific",
        }

        # Run process and assert.
        self.run_process("workflow-mirna", inputs)
        workflow = Data.objects.filter(process__slug="feature_counts").last()

        # check featureCount summary
        self.assertFile(
            workflow, "rc", "mirna_featurecounts_rc.tab.gz", compression="gzip"
        )
        self.assertFile(
            workflow, "exp", "mirna_featurecounts_cpm.tab.gz", compression="gzip"
        )
