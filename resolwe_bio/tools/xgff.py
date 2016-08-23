#!/usr/bin/env python2
# pylint: disable=invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Create gene expression profiles.

All transcripts of a gene are of the same type ('mRNA' or 'tRNA' etc.)
Each gene has a 'dictyBase Curator' transcript (among many) or (if no
curator transcript) it has exactly one non-curated transcript.
Exons in each transcript do NOT overlap.
The result of the parse is ONE transcript per gene.
"""
import argparse
import sys


gene_types = ['mRNA', 'SRP_RNA', 'snRNA', 'tRNA', 'C_D_box_snoRNA', 'class_I_RNA',
              'class_II_RNA', 'rRNA', 'H_ACA_box_snoRNA', 'RNase_P_RNA', 'ncRNA']

parser = argparse.ArgumentParser(description='Create gene expression profiles.')
parser.add_argument('gff_file', help='gff file')
parser.add_argument('chr', help='chromosome')
parser.add_argument('--gene_types', default=gene_types, nargs='+', help='')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')

args = parser.parse_args()

genes = {}
tr_genes = {}
gene_type = None
regstart = sys.maxsize
regstop = 0

with open(args.gff_file, 'rt') as f:
    r = f.readline()

    while r:
        r = r.replace('\r', '').replace('\n', '').split('\t')
        if len(r) < 3:
            r = f.readline()
            continue
        seqid = r[0]
        source = r[1]
        rec_type = r[2]
        if rec_type not in ['gene', 'exon', 'pseudogene'] + gene_types:
            r = f.readline()
            continue
        start = int(r[3])
        stop = int(r[4])
        regstart = min(regstart, start)
        regstop = max(regstop, stop)

        strand = r[6]
        attr = {}
        for item in r[-1].split(';'):
            item = item.split('=')
            attr[item[0]] = item[1]
        if rec_type == 'gene':
            gene_id = attr['ID']
            genes[gene_id] = {'chr': args.chr, 'strand': strand, 'attr': attr, 'gene_type': '', 'transcripts': {}}
        if rec_type in gene_types:
            gene_id = attr['Parent']
            t_id = attr['ID']
            gene_type = rec_type
            genes[gene_id]['gene_type'] = gene_type
            genes[gene_id]['transcripts'][t_id] = {}
            genes[gene_id]['transcripts'][t_id]['source'] = source
            genes[gene_id]['transcripts'][t_id]['type'] = rec_type
            genes[gene_id]['transcripts'][t_id]['exons'] = []
            tr_genes[t_id] = gene_id
        if rec_type == 'pseudogene':
            gene_type = None
        if gene_type is not None and rec_type == 'exon':
            t_id = attr['Parent']
            gene_id = tr_genes[t_id]
            genes[gene_id]['transcripts'][t_id]['exons'].append((start, stop))
        r = f.readline()

f_gff3 = open('%s.gff3' % args.chr, 'wt')
f_gff3.write('##sequence-region %s %d %d\n' % (args.chr, regstart, regstop))

nocodons = set(['tRNA', 'snRNA', 'class_I_RNA', 'C_D_box_snoRNA', 'class_II_RNA'])

for gene_id, data in genes.items():
    gene_ID = data['attr']['ID']
    gene_type = data['gene_type']
    desc = data['attr'].get('description', '')
    name = data['attr']['Name']
    strand = data['strand']
    _chr = data['chr']
    transcripts = data['transcripts']
    tr_types = set([t_data['type'] for t_data in transcripts.values()])
    assert len(tr_types) <= 1  # all transcripts of a gene need to be of the same type, i.e. mRNA, tRNA etc.
    # sort transcripts by source
    for t_id, t_data in transcripts.items():
        transcripts[t_data['source']] = (t_id, t_data)
        del transcripts[t_id]
    if len(transcripts.keys()) == 0:
        continue
    if 'dictyBase Curator' not in transcripts.keys():
        assert len(transcripts.keys()) == 1
    # or curated or the single non-curated transcript:
    t_id, t_data = transcripts.get('dictyBase Curator', transcripts[transcripts.keys()[0]])
    t_type = t_data['type']
    t_source = t_data['source']
    exons = t_data['exons']

    if len(exons) > 0:
        phases = []
        exons.sort(reverse=True if strand == '-' else False)
        if gene_type not in nocodons:
            phases.append('0')
            phase = 0
            for (exon_start, exon_stop) in exons:
                phase = exon_stop - exon_start + 1 + (2 if phase == 1 else (1 if phase == 2 else 0))
                phase %= 3
                phase = 2 if phase == 1 else (1 if phase == 2 else 0)
                phases.append(str(phase))

            phases = phases[:-1]
            if strand == '-':
                phases = list(reversed(phases))

        exons.sort()
        attrs = 'ID=%s;' % gene_ID
        if name:
            attrs += 'Name=%s;' % name
        if desc:
            attrs += 'Note=%s;' % desc

        row = [_chr, 'biox', gene_type, str(exons[0][0]), str(exons[-1][1]), '0', strand, '.', attrs[:-1]]
        f_gff3.write('\t'.join(row) + '\n')
        if gene_type not in nocodons:
            for i, (exon_start, exon_stop) in enumerate(exons):
                row = [_chr, 'biox', 'CDS', str(exon_start), str(exon_stop),
                       '0', strand, phases[i], 'ID=cds%s;Parent=%s' % (gene_ID, gene_ID)]
                f_gff3.write('\t'.join(row) + '\n')

f_gff3.close()
