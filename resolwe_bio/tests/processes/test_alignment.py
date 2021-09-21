import os
from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import KBBioProcessTestCase


class AlignmentProcessorTestCase(KBBioProcessTestCase):
    @tag_process("bowtie-index", "alignment-bowtie")
    def test_bowtie(self):
        input_folder = Path("test_bowtie") / "input"
        output_folder = Path("test_bowtie") / "output"
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn=str(input_folder / "g.en ome.fasta.gz"),
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )

            reads_single = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        bowtie_index = self.run_process("bowtie-index", {"ref_seq": ref_seq.id})
        self.assertDir(bowtie_index, "index", output_folder / "bowtie_index.tar.gz")
        self.assertFile(bowtie_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            bowtie_index,
            "fastagz",
            output_folder / "genome.fasta.gz",
            compression="gzip",
        )
        self.assertFile(
            bowtie_index,
            "fai",
            output_folder / "genome.fasta.fai",
        )
        self.assertFields(bowtie_index, "species", "Dictyostelium discoideum")
        self.assertFields(bowtie_index, "build", "dd-05-2009")

        alignment_single = self.run_process(
            "alignment-bowtie",
            {
                "genome": bowtie_index.id,
                "reads": reads_single.id,
                "trimming": {"trim_iter": 2, "trim_nucl": 4},
                "reporting": {"r": "-a -m 1 --best --strata"},
            },
        )
        self.assertFile(
            alignment_single,
            "stats",
            output_folder / "bowtie_single_reads_report.tab.gz",
            compression="gzip",
        )
        self.assertFields(alignment_single, "species", "Dictyostelium discoideum")
        self.assertFields(alignment_single, "build", "dd-05-2009")

        alignment_paired = self.run_process(
            "alignment-bowtie",
            {
                "genome": bowtie_index.id,
                "reads": reads_paired.id,
                "trimming": {"trim_iter": 2, "trim_nucl": 4},
                "reporting": {"r": "-a -m 1 --best --strata"},
            },
        )
        self.assertFile(
            alignment_paired,
            "stats",
            output_folder / "bowtie_paired_reads_report.tab.gz",
            compression="gzip",
        )
        sample = Sample.objects.get(data=alignment_paired)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

        alignment_use_se = self.run_process(
            "alignment-bowtie",
            {
                "genome": bowtie_index.id,
                "reads": reads_paired.id,
                "trimming": {"trim_iter": 2, "trim_nucl": 4},
                "reporting": {"r": "-a -m 1 --best --strata"},
                "use_se": True,
            },
        )
        self.assertFile(
            alignment_use_se,
            "stats",
            output_folder / "bowtie_use_SE_report.tab.gz",
            compression="gzip",
        )

    @tag_process("bowtie2-index", "alignment-bowtie2")
    def test_bowtie2(self):
        input_folder = Path("test_bowtie2") / "input"
        output_folder = Path("test_bowtie2") / "output"
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn=str(input_folder / "g.en ome.fasta.gz"),
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        bowtie2_index = self.run_process("bowtie2-index", {"ref_seq": ref_seq.id})
        self.assertDir(bowtie2_index, "index", output_folder / "bowtie2_index.tar.gz")
        self.assertFile(bowtie2_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            bowtie2_index,
            "fastagz",
            output_folder / "genome.fasta.gz",
            compression="gzip",
        )
        self.assertFile(bowtie2_index, "fai", output_folder / "genome.fasta.fai")
        self.assertFields(bowtie2_index, "species", "Dictyostelium discoideum")
        self.assertFields(bowtie2_index, "build", "dd-05-2009")

        # Values for alignment options are default according to the documentation. However, L may not be set
        # correctly as there is some incongruency. See https://github.com/BenLangmead/bowtie2/issues/215
        inputs = {
            "genome": bowtie2_index.pk,
            "reads": reads.pk,
            "trimming": {"trim_iter": 2, "trim_nucl": 4},
            "reporting": {"rep_mode": "def"},
            "alignment_options": {
                "N": 0,
                "gbar": 4,
                "L": 22,
                "mp": "6",
                "rdg": "5,3",
                "rfg": "5,3",
                "score_min": "L,-0.6,-0.6",
            },
            "misc_opts": {
                "bw_binsize": 50,
                "bw_timeout": 30,
            },
        }
        single_end = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(single_end, "stats", output_folder / "bowtie2_reads_report.txt")
        self.assertFields(single_end, "species", "Dictyostelium discoideum")
        self.assertFields(single_end, "build", "dd-05-2009")

        inputs["reads"] = reads_paired.id
        paired_end = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(
            paired_end, "stats", output_folder / "bowtie2_paired_end_report.txt"
        )
        sample = Sample.objects.get(data=paired_end)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

        inputs["PE_options"] = {"use_se": True}
        paired_end_se_mode = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(
            paired_end_se_mode, "stats", output_folder / "bowtie2_use_SE_report.txt"
        )

        inputs["PE_options"] = {"no_overlap": True}  # this overwrites use_se from above
        paired_end_no_overlap = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(
            paired_end_no_overlap,
            "stats",
            output_folder / "bowtie2_no_overlap_report.txt",
        )

        inputs["PE_options"]["no_overlap"] = False
        inputs["output_opts"] = {"no_unal": True}
        paired_end_no_unal = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(
            paired_end_no_unal, "stats", output_folder / "bowtie2_no_unal_report.txt"
        )

        del inputs["output_opts"]
        inputs["PE_options"]["dovetail"] = True
        paired_end_dove = self.run_process("alignment-bowtie2", inputs)
        self.assertFile(
            paired_end_dove, "stats", output_folder / "bowtie2_paired_end_report.txt"
        )

    @with_resolwe_host
    @tag_process("alignment-star-index", "alignment-star")
    def test_star(self):
        with self.preparation_stage():
            reads = self.prepare_reads(
                [
                    os.path.join(
                        "test_star",
                        "input",
                        "hs_single bbduk_star_htseq_reads_single.fastq.gz",
                    )
                ]
            )
            paired_reads = self.prepare_paired_reads(
                mate1=[
                    os.path.join(
                        "test_star",
                        "input",
                        "hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz",
                    )
                ],
                mate2=[
                    os.path.join(
                        "test_star",
                        "input",
                        "hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz",
                    )
                ],
            )
            annotation = self.prepare_annotation(
                os.path.join("test_star", "input", "hs annotation.gtf.gz"),
                source="ENSEMBL",
                species="Homo sapiens",
                build="GRCh38_ens90",
            )

            inputs = {
                "src": os.path.join("test_star", "input", "hs genome.fasta.gz"),
                "species": "Homo sapiens",
                "build": "GRCh38_ens90",
            }
            star_index_fasta = self.run_process("upload-fasta-nucl", inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        gene_counts_star_single = os.path.join(
            "test_star", "output", "gene_counts_star_single.tab.gz"
        )
        star_expression_single = os.path.join(
            "test_star", "output", "star_expression_single.tab.gz"
        )

        # prepare genome indices
        star_index = self.run_process(
            "alignment-star-index",
            {"annotation": annotation.id, "ref_seq": star_index_fasta.id},
        )

        star_index_wo_annot = self.run_process(
            "alignment-star-index",
            {
                "ref_seq": star_index_fasta.id,
                "source": "ENSEMBL",
            },
        )

        # test STAR alignment
        inputs = {
            "genome": star_index.id,
            "reads": reads.id,
            "t_coordinates": {
                "quantmode": True,
                "gene_counts": True,
            },
            "two_pass_mapping": {
                "two_pass_mode": True,
            },
            "detect_chimeric": {
                "chimeric": True,
            },
        }
        aligned_reads = self.run_process("alignment-star", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_single, compression="gzip"
        )
        self.assertFields(aligned_reads, "species", "Homo sapiens")
        self.assertFields(aligned_reads, "build", "GRCh38_ens90")
        sample = Sample.objects.get(data=aligned_reads)
        self.assertEqual(sample.descriptor["general"]["species"], "Homo sapiens")

        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_single, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")
        self.assertFields(exp, "build", "GRCh38_ens90")
        self.assertFields(exp, "feature_type", "gene")

        inputs["star_sort"] = True
        aligned_reads = self.run_process("alignment-star", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_single, compression="gzip"
        )
        self.assertFields(aligned_reads, "species", "Homo sapiens")
        self.assertFields(aligned_reads, "build", "GRCh38_ens90")

        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_single, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")
        self.assertFields(exp, "build", "GRCh38_ens90")
        self.assertFields(exp, "feature_type", "gene")

        inputs["genome"] = star_index_wo_annot.id
        inputs["annotation"] = annotation.id
        aligned_reads = self.run_process("alignment-star", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_single, compression="gzip"
        )
        self.assertFields(aligned_reads, "species", "Homo sapiens")
        self.assertFields(aligned_reads, "build", "GRCh38_ens90")

        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_single, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")
        self.assertFields(exp, "build", "GRCh38_ens90")
        self.assertFields(exp, "feature_type", "gene")

        inputs["star_sort"] = False
        aligned_reads = self.run_process("alignment-star", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_single, compression="gzip"
        )
        self.assertFields(aligned_reads, "species", "Homo sapiens")
        self.assertFields(aligned_reads, "build", "GRCh38_ens90")

        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_single, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")
        self.assertFields(exp, "build", "GRCh38_ens90")
        self.assertFields(exp, "feature_type", "gene")

        bigwig_star_ens_paired = os.path.join(
            "test_star", "output", "bigwig_star_ens_paired.bw"
        )
        gene_counts_star_paired = os.path.join(
            "test_star", "output", "gene_counts_star_paired.tab.gz"
        )
        star_expression_paired = os.path.join(
            "test_star", "output", "star_expression_paired.tab.gz"
        )
        star_out_exp_set = os.path.join(
            "test_star", "output", "star_out_exp_set.txt.gz"
        )
        star_exp_set = os.path.join("test_star", "output", "star_exp_set.json.gz")

        inputs = {
            "genome": star_index.id,
            "reads": paired_reads.id,
            "t_coordinates": {
                "quantmode": True,
                "gene_counts": True,
            },
            "two_pass_mapping": {
                "two_pass_mode": True,
            },
            "star_sort": True,
        }
        aligned_reads = self.run_process("alignment-star", inputs)
        self.assertFile(aligned_reads, "bigwig", bigwig_star_ens_paired)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_paired, compression="gzip"
        )
        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_paired, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFile(exp, "exp_set", star_out_exp_set, compression="gzip")
        self.assertJSON(exp, exp.output["exp_set_json"], "", star_exp_set)

        inputs["star_sort"] = False
        aligned_reads = self.run_process("alignment-star", inputs)
        self.assertFile(aligned_reads, "bigwig", bigwig_star_ens_paired)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFile(
            aligned_reads, "gene_counts", gene_counts_star_paired, compression="gzip"
        )
        exp = Data.objects.last()
        self.assertFile(exp, "exp", star_expression_paired, compression="gzip")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFile(exp, "exp_set", star_out_exp_set, compression="gzip")
        self.assertJSON(exp, exp.output["exp_set_json"], "", star_exp_set)

    @tag_process(
        "bwa-index", "alignment-bwa-aln", "alignment-bwa-sw", "alignment-bwa-mem"
    )
    def test_bwa(self):
        input_folder = Path("test_bwa") / "input"
        output_folder = Path("test_bwa") / "output"
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn=str(input_folder / "g.en ome.fasta.gz"),
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})
        self.assertDir(bwa_index, "index", output_folder / "bwa_index.tar.gz")
        self.assertFile(bwa_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            bwa_index, "fastagz", output_folder / "genome.fasta.gz", compression="gzip"
        )
        self.assertFile(bwa_index, "fai", output_folder / "genome.fasta.fai")
        self.assertFields(bwa_index, "species", "Dictyostelium discoideum")
        self.assertFields(bwa_index, "build", "dd-05-2009")

        single_end_aln = self.run_process(
            "alignment-bwa-aln", {"genome": bwa_index.id, "reads": reads.id}
        )
        self.assertFile(
            single_end_aln, "stats", output_folder / "bwa_bt_reads_report.txt"
        )

        paired_end_aln = self.run_process(
            "alignment-bwa-aln", {"genome": bwa_index.id, "reads": reads_paired.id}
        )
        self.assertFile(
            paired_end_aln, "stats", output_folder / "bwa_bt_paired_reads_report.txt"
        )
        self.assertFields(paired_end_aln, "species", "Dictyostelium discoideum")
        self.assertFields(paired_end_aln, "build", "dd-05-2009")
        sample = Sample.objects.get(data=paired_end_aln)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

        single_end_sw = self.run_process(
            "alignment-bwa-sw", {"genome": bwa_index.id, "reads": reads.id}
        )
        self.assertFile(single_end_sw, "bam", output_folder / "bwa_sw_reads_mapped.bam")
        self.assertFile(
            single_end_sw, "stats", output_folder / "bwa_sw_reads_report.txt"
        )

        paired_end_sw = self.run_process(
            "alignment-bwa-sw", {"genome": bwa_index.id, "reads": reads_paired.id}
        )
        self.assertFile(
            paired_end_sw, "bam", output_folder / "bwa_sw_paired_reads_mapped.bam"
        )
        self.assertFile(
            paired_end_sw, "stats", output_folder / "bwa_sw_paired_reads_report.txt"
        )
        self.assertFields(paired_end_sw, "species", "Dictyostelium discoideum")
        self.assertFields(paired_end_sw, "build", "dd-05-2009")
        sample = Sample.objects.get(data=paired_end_sw)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

        single_end_mem = self.run_process(
            "alignment-bwa-mem", {"genome": bwa_index.id, "reads": reads.id}
        )
        self.assertFile(
            single_end_mem, "stats", output_folder / "bwa_mem_reads_report.txt"
        )

        paired_end_mem = self.run_process(
            "alignment-bwa-mem", {"genome": bwa_index.id, "reads": reads_paired.id}
        )
        self.assertFile(
            paired_end_mem, "stats", output_folder / "bwa_mem_paired_reads_report.txt"
        )
        self.assertFile(
            paired_end_mem,
            "unmapped",
            output_folder / "bwa_mem_unmapped_reads.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFields(paired_end_mem, "species", "Dictyostelium discoideum")
        self.assertFields(paired_end_mem, "build", "dd-05-2009")
        sample = Sample.objects.get(data=paired_end_mem)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

    @tag_process("bwamem2-index", "alignment-bwa-mem2")
    def test_bwa2(self):
        input_folder = Path("test_bwa") / "input"
        output_folder = Path("test_bwa") / "output"
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn=str(input_folder / "g.en ome.fasta.gz"),
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        bwa2_index = self.run_process("bwamem2-index", {"ref_seq": ref_seq.id})
        self.assertDir(bwa2_index, "index", output_folder / "bwa-mem2_index.tar.gz")
        self.assertFile(bwa2_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            bwa2_index, "fastagz", output_folder / "genome.fasta.gz", compression="gzip"
        )
        self.assertFile(bwa2_index, "fai", output_folder / "genome.fasta.fai")
        self.assertFields(bwa2_index, "species", "Dictyostelium discoideum")
        self.assertFields(bwa2_index, "build", "dd-05-2009")

        single_end_mem = self.run_process(
            "alignment-bwa-mem2", {"genome": bwa2_index.id, "reads": reads.id}
        )
        self.assertFile(
            single_end_mem, "stats", output_folder / "bwa_mem_reads_report.txt"
        )

        paired_end_mem = self.run_process(
            "alignment-bwa-mem2", {"genome": bwa2_index.id, "reads": reads_paired.id}
        )
        self.assertFile(
            paired_end_mem, "stats", output_folder / "bwa_mem_paired_reads_report.txt"
        )
        self.assertFile(
            paired_end_mem,
            "unmapped",
            output_folder / "bwa_mem_unmapped_reads.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFields(paired_end_mem, "species", "Dictyostelium discoideum")
        self.assertFields(paired_end_mem, "build", "dd-05-2009")
        sample = Sample.objects.get(data=paired_end_mem)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

    @tag_process("upload-bwamem2-index", "alignment-bwa-mem2")
    def test_bwa2_upload(self):
        input_folder = Path("test_bwa") / "input"
        output_folder = Path("test_bwa") / "output"
        with self.preparation_stage():
            reads = self.prepare_reads()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        bwa2_index = self.run_process(
            "upload-bwamem2-index",
            {
                "index_name": input_folder / "bwa-mem2_index.tar.gz",
                "ref_seq": input_folder / "g.en ome.fasta.gz",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            },
        )

        self.assertFile(bwa2_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            bwa2_index, "fastagz", output_folder / "genome.fasta.gz", compression="gzip"
        )
        self.assertFile(bwa2_index, "fai", output_folder / "genome.fasta.fai")
        self.assertFields(bwa2_index, "species", "Dictyostelium discoideum")
        self.assertFields(bwa2_index, "build", "dd-05-2009")

        single_end_mem = self.run_process(
            "alignment-bwa-mem2", {"genome": bwa2_index.id, "reads": reads.id}
        )
        self.assertFile(
            single_end_mem, "stats", output_folder / "bwa_mem_reads_report.txt"
        )

        paired_end_mem = self.run_process(
            "alignment-bwa-mem2", {"genome": bwa2_index.id, "reads": reads_paired.id}
        )
        self.assertFile(
            paired_end_mem, "stats", output_folder / "bwa_mem_paired_reads_report.txt"
        )
        self.assertFile(
            paired_end_mem,
            "unmapped",
            output_folder / "bwa_mem_unmapped_reads.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFields(paired_end_mem, "species", "Dictyostelium discoideum")
        self.assertFields(paired_end_mem, "build", "dd-05-2009")
        sample = Sample.objects.get(data=paired_end_mem)
        self.assertEqual(
            sample.descriptor["general"]["species"], "Dictyostelium discoideum"
        )

    @tag_process("hisat2-index", "alignment-hisat2")
    def test_hisat2(self):
        input_folder = Path("test_hisat2") / "input"
        output_folder = Path("test_hisat2") / "output"
        with self.preparation_stage():
            # Use dicty genome and reads but declare it as human so no mapping to UCSC chr names can happen.
            ref_seq = self.prepare_ref_seq(
                fn=str(input_folder / "g.enome.fasta.gz"),
                species="Homo sapiens",
                build="GRCh38.p12",
            )
            reads = self.prepare_reads()
            sample = Sample.objects.get(data=reads)
            sample.name = "Single reads"
            sample.save()
            reads_paired = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )
            sample_paired = Sample.objects.get(data=reads_paired)
            sample_paired.name = "Paired-end reads"
            sample_paired.save()
            no_mapping_reads = self.prepare_reads(
                fn=[str(input_folder / "reads-map-to-nowhere.fastq.gz")]
            )

        hisat2_index = self.run_process("hisat2-index", {"ref_seq": ref_seq.id})
        self.assertDir(hisat2_index, "index", output_folder / "hisat2_index.tar.gz")
        self.assertFile(hisat2_index, "fasta", output_folder / "genome.fasta")
        self.assertFile(
            hisat2_index,
            "fastagz",
            output_folder / "genome.fasta.gz",
            compression="gzip",
        )
        self.assertFile(hisat2_index, "fai", output_folder / "genome.fasta.fai")
        self.assertFields(hisat2_index, "species", "Homo sapiens")
        self.assertFields(hisat2_index, "build", "GRCh38.p12")

        single_end = self.run_process(
            "alignment-hisat2", {"genome": hisat2_index.id, "reads": reads.id}
        )
        self.assertFile(single_end, "stats", output_folder / "hisat2_report.txt")
        self.assertFile(
            single_end,
            "unmapped_f",
            output_folder / "hisat2_unmapped.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFileExists(single_end, "splice_junctions")
        msg = "Neither of the chromosomes in the input file has a valid UCSC pair. No mapping will be done."
        self.assertEqual(single_end.process_warning, [msg])

        paired_end = self.run_process(
            "alignment-hisat2", {"genome": hisat2_index.id, "reads": reads_paired.id}
        )
        self.assertFile(paired_end, "stats", output_folder / "hisat2_paired_report.txt")
        self.assertFile(paired_end, "bigwig", output_folder / "hisat2_paired_bigwig.bw")
        self.assertFile(
            paired_end,
            "unmapped_f",
            output_folder / "hisat2_unmapped_1.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFile(
            paired_end,
            "unmapped_r",
            output_folder / "hisat2_unmapped_2.fastq.gz",
            compression="gzip",
            sort=True,
        )
        self.assertFileExists(paired_end, "splice_junctions")
        self.assertFields(paired_end, "species", "Homo sapiens")
        self.assertFields(paired_end, "build", "GRCh38.p12")

        no_mapping = self.run_process(
            "alignment-hisat2",
            {"reads": no_mapping_reads.id, "genome": hisat2_index.id},
        )
        self.assertEqual(
            no_mapping.process_warning,
            ["Bam file has no entries. No bigWig file will be made."],
        )

    @tag_process("alignment-hisat2")
    def test_hisat2_bigwig(self):
        input_folder = Path("test_hisat2") / "input"
        output_folder = Path("test_hisat2") / "output"
        with self.preparation_stage():
            paired_reads_ucsc = self.prepare_paired_reads(
                mate1=["SRR2124780_1 1k.fastq.gz"], mate2=["SRR2124780_2 1k.fastq.gz"]
            )
            paired_reads_ncbi = self.prepare_paired_reads(
                mate1=["GRCh38.p12_NCBIchr21_R1.fq.gz"],
                mate2=["GRCh38.p12_NCBIchr21_R2.fq.gz"],
            )

            ref_seq_ucsc = self.prepare_ref_seq(
                fn=str(input_folder / "hg38_chr21_9M.fa.gz"),
                species="Homo sapiens",
                build="hg38",
            )
            hisat2_index_ucsc = self.run_process(
                "hisat2-index", {"ref_seq": ref_seq_ucsc.id}
            )

            ref_seq_ncbi = self.prepare_ref_seq(
                fn=str(input_folder / "GRCh38.p12_NCBIchr21_9M.fasta.gz"),
                species="Homo sapiens",
                build="GRCh38.p12",
            )
            hisat2_index_ncbi = self.run_process(
                "hisat2-index", {"ref_seq": ref_seq_ncbi.id}
            )

        aligned_reads = self.run_process(
            "alignment-hisat2",
            {"genome": hisat2_index_ucsc.id, "reads": paired_reads_ucsc.id},
        )
        self.assertFile(
            aligned_reads, "bigwig", output_folder / "hisat2_paired_ucsc_bigwig.bw"
        )

        aligned_reads = self.run_process(
            "alignment-hisat2",
            {"genome": hisat2_index_ncbi.id, "reads": paired_reads_ncbi.id},
        )
        self.assertFile(
            aligned_reads, "bigwig", output_folder / "hisat2_paired_ncbi_bigwig.bw"
        )
        sample = Sample.objects.get(data=aligned_reads)
        self.assertEqual(sample.descriptor["general"]["species"], "Homo sapiens")
