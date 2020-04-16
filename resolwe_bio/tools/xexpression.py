#!/usr/bin/env python2
# XXX: Refactor to a comand line tool and remove pylint disable
"""Create gene expression profiles."""
from __future__ import absolute_import, division, print_function

import argparse
import gzip
import logging
import math
import os

import biox

parser = argparse.ArgumentParser(description="Create gene expression profiles.")
parser.add_argument("gff_file", help="GFF file")
parser.add_argument("bam_file", help="aligned BAM file")
parser.add_argument(
    "--rc", action="store_true", help="reads uniquely mapped to gene exons"
)
parser.add_argument("--rpkm", action="store_true", help="reads scaled by exon length")
parser.add_argument(
    "--rpkum", help="reads scaled by uniquely mappable part of exons <mappability_file>"
)
parser.add_argument(
    "--mrna", action="store_true", help="scale with reads that map to polyA transcripts"
)
parser.add_argument(
    "--ncrna",
    action="store_true",
    help="Exclude reads that map to chrR from scaling factor",
)
parser.add_argument("--stranded", action="store_true", help="Stranded library type")
parser.add_argument("-v", "--verbose", action="store_true", help="verbose output")

args = parser.parse_args()

if os.path.splitext(args.bam_file)[1] != ".bam":
    raise ValueError("Expected .bam file, got {}.".format(args.bam_file))

if os.path.splitext(args.gff_file)[1] not in [".gff", ".gff3"]:
    raise ValueError("Expected .gff file, got {}.".format(args.gff_file))

if args.verbose:
    biox.utils.verbosity(logging.INFO)

gff_file = args.gff_file
gtf_file = "foo.gtf"
bam_file = args.bam_file
suffix = "_polya" if args.mrna else ""

f_gtf = open(gtf_file, "w")
parents = {}
with open(gff_file) as f:
    for line in f:
        if line.strip().startswith("#"):
            continue
        vals = line.split("\t")
        attrlist = vals[-1].strip().split(";")
        attrs = {}
        for att in attrlist:
            key, val = att.split("=")
            attrs[key] = val

        if vals[2] == "CDS" or vals[2] == "exon":
            gene_id = attrs["Parent"]
            f_gtf.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    vals[0],
                    "",
                    "exon",
                    vals[3],
                    vals[4],
                    ".",
                    vals[6],
                    ".",
                    'gene_id "{}"; gene_name "{}"; gene_description "{}"; gene_type "{}";'.format(
                        gene_id,
                        parents[gene_id]["name"],
                        parents[gene_id]["description"],
                        parents[gene_id]["gene_type"],
                    ),
                )
            )
        else:
            parents[attrs["ID"]] = {"gene_type": vals[2]}
            parents[attrs["ID"]]["name"] = (
                attrs["Name"] if "Name" in attrs else attrs["ID"]
            )
            parents[attrs["ID"]]["description"] = (
                attrs["Note"] if "Note" in attrs else ""
            )

f_gtf.close()

gtf = biox.data.Gtf(gtf_file)

genes = None
if args.mrna:
    genes = set()
    for gene_id, gene in gtf.genes.items():
        if gene.attrs["gene_type"] == "mRNA" or gene.attrs["gene_type"] == "'mRNA'":
            genes.add(gene_id)


def gene_expression_overlap_stranded(gtf_file, bam_file, quality=30):
    """Compute gene expression overlap."""
    gtf = biox.data.Gtf(gtf_file)
    genes_exp = {}

    for gene_id in gtf.genes:
        genes_exp[gene_id] = 0

    current = 0
    for gene_id, gene in gtf.genes.items():
        current += 1
        if current % 300 == 0:
            print("%.2f" % (float(current) / len(gtf.genes)), bam_file)
        for feature in gene.features:
            if feature.type != "exon":
                continue
            assert feature.start <= feature.stop
            if gene.strand == "-":
                # warning: for std. illumina library prep, this statements would count reads in "plus" strand
                # alignments of the second in pair if they map to the forward strand
                command = "{samtools} view -f 128 -F 16 -q {quality} -c {bam_file} {chr}:{start}-{stop}".format(
                    samtools=os.path.join(biox.samtools_folder, "samtools"),
                    bam_file=bam_file,
                    quality=30,
                    chr=gene.chr,
                    start=feature.start,
                    stop=feature.stop,
                )
                output_second_in_pair, _ = biox.utils.cmd(command)
                output_second_in_pair = (
                    output_second_in_pair if output_second_in_pair != "" else 0
                )

                # alignments of the first in pair if they map to the reverse  strand
                command = "{samtools} view -F 80 -q {quality} -c {bam_file} {chr}:{start}-{stop}".format(
                    samtools=os.path.join(biox.samtools_folder, "samtools"),
                    bam_file=bam_file,
                    quality=30,
                    chr=gene.chr,
                    start=feature.start,
                    stop=feature.stop,
                )
                output_first_in_pair, _ = biox.utils.cmd(command)
                output_first_in_pair = (
                    output_first_in_pair if output_first_in_pair != "" else 0
                )

                output = int(output_second_in_pair) + int(output_first_in_pair)
                genes_exp[gene_id] = genes_exp.get(gene_id, 0) + int(output)
            else:
                # warning: for std. illumina library prep, this statements would count reads in "minus" strand
                # alignments of the second in pair if they map to the reverse strand
                command = "{samtools} view -f 144 -q {quality} -c {bam_file} {chr}:{start}-{stop}".format(
                    samtools=os.path.join(biox.samtools_folder, "samtools"),
                    bam_file=bam_file,
                    quality=30,
                    chr=gene.chr,
                    start=feature.start,
                    stop=feature.stop,
                )
                output_second_in_pair, _ = biox.utils.cmd(command)
                output_second_in_pair = (
                    output_second_in_pair if output_second_in_pair != "" else 0
                )

                # alignments of the first in pair if they map to the forward strand
                command = "{samtools} view -f 64 -F 16 -q {quality} -c {bam_file} {chr}:{start}-{stop}".format(
                    samtools=os.path.join(biox.samtools_folder, "samtools"),
                    bam_file=bam_file,
                    quality=30,
                    chr=gene.chr,
                    start=feature.start,
                    stop=feature.stop,
                )
                output_first_in_pair, _ = biox.utils.cmd(command)
                output_first_in_pair = (
                    output_first_in_pair if output_first_in_pair != "" else 0
                )

                output = int(output_second_in_pair) + int(output_first_in_pair)
                genes_exp[gene_id] = genes_exp.get(gene_id, 0) + int(output)

    return genes_exp


print("Expression profile overlap...")
if args.stranded:
    results = gene_expression_overlap_stranded(gtf_file, bam_file)
else:
    results = biox.expression.gene_expression_overlap(gtf_file, bam_file, genes=genes)


def save_expression_profile(file_name, exp=results.get):
    """Save expression profile."""
    with gzip.open(file_name, "wb") as f:
        f.write("Gene\tExpression\n")
        gene_ids = results.keys()
        gene_ids.sort()
        for gene_id in gene_ids:
            _exp = str(exp(gene_id))
            f.write("{}\t{}\n".format(gene_id, _exp))


if args.rc:
    print("Writing read counts...")
    save_expression_profile("expression_rc{}.tab.gz".format(suffix))

if args.rpkm or args.rpkum:
    if not args.mrna and not args.ncrna:
        command = "{samtools} view -F 4 -q {quality} -c {bam_file}".format(
            samtools=os.path.join(biox.samtools_folder, "samtools"),
            bam_file=bam_file,
            quality=30,
        )
        output, error = biox.utils.cmd(command)
        N = int(output)
    if args.ncrna:
        command = "{samtools} view -F 4 -q {quality} -c {bam_file}".format(
            samtools=os.path.join(biox.samtools_folder, "samtools"),
            bam_file=bam_file,
            quality=30,
        )
        output, error = biox.utils.cmd(command)
        all_reads = int(output)
        print("Number of all reads: ", all_reads)
        command = "{samtools} view -F 4 -q {quality} -c {bam_file} chrR".format(
            samtools=os.path.join(biox.samtools_folder, "samtools"),
            bam_file=bam_file,
            quality=30,
        )
        output, error = biox.utils.cmd(command)
        chrR_reads = int(output)
        print("Number of reads that map to chrR: ", chrR_reads)
        N = all_reads - chrR_reads
        print("Effective library size: ", N)
    else:
        N = sum(results.values())

if args.rpkm:
    print("Writing RPKM...")
    gene_exon_lens = {}
    for gene_id, gene in gtf.genes.items():
        exon_len = 0
        for feature in gene.features:
            if feature.type != "exon":
                continue
            exon_len += feature.stop - feature.start + 1
        gene_exon_lens[gene_id] = exon_len

    def exp(gene_id):
        """Compute RPKM."""
        return (math.pow(10, 9) * results[gene_id]) / (N * gene_exon_lens[gene_id])

    save_expression_profile("expression_rpkm{}.tab.gz".format(suffix), exp)

if args.rpkum:
    print("Processing RPKUM...")

    # find read length
    mapability_file = args.rpkum
    data_mapability = {}
    f = biox.data.TabReader(mapability_file)
    while f.readline():
        data_mapability[f.gene_id] = f.coverage

    def exp_rpkum(gene_id):
        """Compute RPKUM."""
        if data_mapability[gene_id] == 0:
            return 0
        else:
            return (math.pow(10, 9) * results[gene_id]) / (N * data_mapability[gene_id])

    save_expression_profile("expression_rpkum{}.tab.gz".format(suffix), exp_rpkum)
