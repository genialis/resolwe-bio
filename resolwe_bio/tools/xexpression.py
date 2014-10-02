import argparse
import logging
import math
import os
import gzip

import biox


parser = argparse.ArgumentParser(description='Create gene expression profiles.')
parser.add_argument('gff_file', help='GFF file')
parser.add_argument('bam_file', help='aligned BAM file')
parser.add_argument('--rc', action='store_true', help='reads uniquely mapped to gene exons')
parser.add_argument('--rpkm', action='store_true', help='reads scaled by exon length')
parser.add_argument('--rpkum', help='reads scaled by uniquely mappable part of exons <mappability_file>')
parser.add_argument('--mrna', action='store_true', help='scale with reads that map to polyA transcripts')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')

args = parser.parse_args()

if os.path.splitext(args.bam_file)[1] != '.bam':
    raise ValueError("Expected .bam file, got {}.".format(args.bam_file))

if os.path.splitext(args.gff_file)[1] not in ['.gff', '.gff3']:
    raise ValueError("Expected .gff file, got {}.".format(args.gff_file))

if args.verbose:
    biox.utils.verbosity(logging.INFO)

gff_file = args.gff_file
gtf_file = 'foo.gtf'
bam_file = args.bam_file
suffix = '_polya' if args.mrna else ''

f_gtf = open(gtf_file, 'w')
parents = {}
with open(gff_file) as f:
    for line in f:
        if line.strip().startswith('#'):
            continue
        vals = line.split('\t')
        attrlist = vals[-1].strip().split(';')
        attrs = {}
        for att in attrlist:
            key, val = att.split('=')
            attrs[key] = val

        if vals[2] == 'CDS' or vals[2] == 'exon':
            gene_id = attrs['Parent']
            f_gtf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                vals[0], '', 'exon', vals[3], vals[4], '.', vals[6], '.',
                'gene_id "{}"; gene_name "{}"; gene_description "{}"; gene_type "{}";'.format(
                    gene_id, parents[gene_id]['name'],
                    parents[gene_id]['description'], parents[gene_id]['gene_type'])))
        else:
            parents[attrs['ID']] = {'gene_type': vals[2]}
            parents[attrs['ID']]['name'] = attrs['Name'] if 'Name' in attrs else attrs['ID']
            parents[attrs['ID']]['description'] = attrs['Note'] if 'Note' in attrs else ''

f_gtf.close()

gtf = biox.data.Gtf(gtf_file)

genes = None
if args.mrna:
    genes = set()
    for gene_id, gene in gtf.genes.items():
        if gene.attrs['gene_type'] == 'mRNA' or gene.attrs['gene_type'] == '\'mRNA\'':
            genes.add(gene_id)

print "Expression profile overlap..."
results = biox.expression.gene_expression_overlap(gtf_file, bam_file, genes=genes)


def save_expression_profile(file_name, exp=results.get):
    with gzip.open(file_name, 'wb') as f:
        f.write('Gene\tExpression\n')
        gene_ids = results.keys()
        gene_ids.sort()
        for gene_id in gene_ids:
            _exp = str(exp(gene_id))
            f.write('{}\t{}\n'.format(gene_id, _exp))


if args.rc:
    print "Writing read counts..."
    save_expression_profile('expression_rc{}.tab.gz'.format(suffix))

if args.rpkm or args.rpkum:
    if not args.mrna:
        command = "{samtools} view -F 4 -q {quality} -c {bam_file}".format(
            samtools=os.path.join(biox.samtools_folder, "samtools"), bam_file=bam_file, quality=30)
        output, error = biox.utils.cmd(command)
        N = int(output)
    else:
        N = sum(results.values())

if args.rpkm:
    print "Writing RPKM..."
    gene_exon_lens = {}
    for gene_id, gene in gtf.genes.items():
        exon_len = 0
        for feature in gene.features:
            if feature.type != 'exon':
                continue
            exon_len += (feature.stop - feature.start + 1)
        gene_exon_lens[gene_id] = exon_len

    def exp(gene_id):
        return (math.pow(10, 9) * results[gene_id]) / (N * gene_exon_lens[gene_id])

    save_expression_profile('expression_rpkm{}.tab.gz'.format(suffix), exp)

if args.rpkum:
    print "Processing RPKUM..."

    # find read length
    mapability_file = args.rpkum
    data_mapability = {}
    f = biox.data.TabReader(mapability_file)
    while f.readline():
        data_mapability[f.gene_id] = f.coverage

    def exp(gene_id):
        if data_mapability[gene_id] == 0:
            return 0
        else:
            return (math.pow(10, 9) * results[gene_id]) / (N * data_mapability[gene_id])

    save_expression_profile('expression_rpkum{}.tab.gz'.format(suffix), exp)
