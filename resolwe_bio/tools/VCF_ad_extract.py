#!/usr/bin/env python2
import sys

vcf_file = file(sys.argv[1], 'r')
annotation_file = file(sys.argv[2], 'r')

variant_allel_depth = {}

for entry in vcf_file:
    if entry[0] == '#':
        continue
    entry = entry.split()
    IDs = entry[-2].split(':')
    vals = entry[-1].split(':')
    variant_allel_depth[(entry[0].lower(), entry[1])] = vals[IDs.index('AD')]

annotation_file.next()

print "Allelic depths for the ref and alt alleles"
for entry in annotation_file:
    entry = entry.split()
    print variant_allel_depth.get((entry[0].lower(), entry[1]), "No info")
