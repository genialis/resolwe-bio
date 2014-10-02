#!/usr/bin/env python
import sys
import json
import numpy as np
import os

min_filter = float(sys.argv[1]) if len(sys.argv) > 1 else 0.
variants_f = sys.argv[2] if len(sys.argv) > 2 else None

vars = {}
if variants_f:
    with open(variants_f, 'rt') as variants:
        vars = [l.split()[:3] for l in variants if not l.startswith('#')]

    vars = dict(((ch, pos), name) for ch, pos, name in vars)

exon_coverage = open('exon_coverage.tsv', 'w')
exon_coverage.write(
    "Gene\tTranscript\tChromosome\tExon\tLocation\tCoverage Total\tCoverage Mean\tCoverage Median"
    "\tBases All\tBases Covered (>{0}x)\tBases Covered (>{0}x) Ratio\tVariants All\tVariants "
    "Covered (>{0}x)\tVariants Covered (>{0}x) Ratio\n".format(min_filter))

transcript_coverage = open('transcript_coverage.tsv', 'w')
transcript_coverage.write(
    "Chromosome\tGene\tTranscript\tCoverage Total\tCoverage Mean\tCoverage Median\tBases All\t"
    "Bases Covered (>{0}x)\tBases Covered (>{0}x) Ratio\tExons All\tExons Covered (>{0}x)\t"
    "Exons Covered (>{0}x) Ratio\tVariants All\tVariants Covered (>{0}x)\t"
    "Variants Covered (>{0}x) Ratio\n".format(min_filter))

report = {'report': []}

for exonage in os.listdir('.'):
    if not os.path.isfile(exonage):
        continue

    v = exonage.split('_')
    if len(v) < 2 or v[0] != 'exonage':
        continue

    coverage, coverage_map, exons, transcripts, variantset = [], {}, {}, [], []
    ts = set()
    gene, chromosome = "", ""

    nbases_above_filter = 0

    chrs = set(l.split()[0] for l in open(exonage) if not l.startswith("EXON"))
    keyz = set(k for k in vars.keys() if k[0] in chrs)

    with open(exonage) as f:
        for line in f:
            if not line.startswith("EXON"):
                cov = float(line.split()[-1])
                l = line.split()
                chromosome = l[0]
                if (l[0], l[1]) in keyz:
                    variantset.append({
                        "chromosome": chromosome,
                        "location": l[1],
                        "coverage": cov,
                        "x-ref": vars[l[0], l[1]],
                        "above_filter": True if cov > min_filter else False
                    })

                if cov > min_filter:
                    nbases_above_filter += 1
                coverage.append(cov)

            else:
                l = line.split()
                coverage_total = sum(coverage)
                bases = int(l[-1]) + 1
                gene = l[1].strip("\";")
                nvariants = len(variantset)
                nvariants_above_filter = len([v for v in variantset if v["above_filter"]])

                exons[l[3]] = ({
                    "chromosome": chromosome,
                    "gene": gene,
                    "transcript": l[2].strip("\";"),
                    "location": l[4],
                    "coverage_total": int(coverage_total),
                    "coverage_max": int(max(coverage or [0])),
                    "coverage_mean": (coverage_total / float(bases)) if bases else 0.,
                    "coverage_median": np.median(coverage) if sum(coverage) > 0 else 0.,
                    "bases_all": bases,
                    "bases_above_filter": nbases_above_filter,
                    "bases_above_filter_ratio": (nbases_above_filter / float(bases)) if bases else 0.,
                    "variants_all": nvariants,
                    "variants_above_filter": nvariants_above_filter,
                    "variants_above_filter_ratio": (nvariants_above_filter / float(nvariants)) if nvariants else 0.,
                    "variants": variantset
                })

                coverage_map[l[3]] = coverage
                ts.add(l[2].strip("\";"))
                coverage, variantset = [], []
                nbases_above_filter = 0

    for e in ts:
        coverages = []
        coverage_max, coverage_total, variants_all, variants_above_filter = 0, 0, 0, 0
        bases_all, bases_above_filter, cov_exons, exons_all = 0, 0, 0, 0
        genelocfrom, genelocto = sys.maxsize, 0
        gene, chromosome = '', ''
        gene_exons = [(exon, meta) for exon, meta in exons.iteritems() if e == meta["transcript"]]

        for exon, meta in gene_exons:
            bases_all += meta['bases_all']
            bases_above_filter += meta['bases_above_filter']
            coverage_max = max(coverage_max, meta['coverage_total'])
            coverage_total += meta['coverage_total']
            cov_exons += 1 if meta['bases_above_filter'] > 0 else 0
            exons_all += 1
            variants_above_filter += meta['variants_above_filter']
            variants_all += meta['variants_all']
            gene = meta['gene']
            chromosome = meta['chromosome']
            coverages.extend(coverage_map[exon])
            locfrom, locto = meta['location'].split('-')
            genelocfrom = min(genelocfrom, int(locfrom))
            genelocto = max(genelocto, int(locto))

        transcripts.append({
            'chromosome': chromosome,
            'location': '{}-{}'.format(genelocfrom, genelocto),
            'gene': gene,
            'transcript': e,
            'exons': [exons[j] for j in exons if exons[j]['transcript'] == e],
            'coverage_total': coverage_total,
            'coverage_max': coverage_max,
            'coverage_mean': coverage_total / float(bases_all),
            'coverage_median': np.median(coverages) if sum(coverages) > 0 else 0.,
            'bases_all': bases_all,
            'bases_above_filter': bases_above_filter,
            'bases_above_filter_ratio': (bases_above_filter / float(bases_all)) if bases_all else 0.,
            'exons_all': exons_all,
            'exons_above_filter': cov_exons,
            'exons_above_filter_ratio': (cov_exons / float(exons_all)) if exons_all else 0.,
            'variants_all': variants_all,
            'variants_above_filter': variants_above_filter,
            'variants_above_filter_ratio': ((variants_above_filter / float(variants_all))
                                            if variants_above_filter else 0.),
        })

    report['report'].extend(transcripts)

    for exon, meta in exons.iteritems():
        exon_coverage.write("\t".join([
            meta['gene'],
            meta['transcript'],
            meta['chromosome'],
            exon,
            meta['location'],
            str(meta['coverage_total']),
            str(meta['coverage_mean']),
            str(meta['coverage_median']),
            str(meta['bases_all']),
            str(meta['bases_above_filter']),
            str(meta['bases_above_filter_ratio']),
            str(meta['variants_all']),
            str(meta['variants_above_filter']),
            str(meta['variants_above_filter_ratio'])]) + "\n")

    for meta in transcripts:
        transcript_coverage.write("\t".join([
            meta['chromosome'],
            meta['gene'],
            meta['transcript'],
            str(meta['coverage_total']),
            str(meta['coverage_mean']),
            str(meta['coverage_median']),
            str(meta['bases_all']),
            str(meta['bases_above_filter']),
            str(meta['bases_above_filter_ratio']),
            str(meta['exons_all']),
            str(meta['exons_above_filter']),
            str(meta['exons_above_filter_ratio']),
            str(meta['variants_all']),
            str(meta['variants_above_filter']),
            str(meta['variants_above_filter_ratio'])]) + "\n")


exon_coverage.close()
transcript_coverage.close()
print json.dumps({'report': report}, separators=(',', ':'))
