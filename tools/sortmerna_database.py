#!/usr/bin/env python2
import argparse

parser = argparse.ArgumentParser(description="Prepare rRNA database paths")
parser.add_argument('databases', help="Comma separated list of rRNA databases")
parser.add_argument('root', help="Root path")
parser.add_argument('slugs_path', help="Slugs path")
args = parser.parse_args()

database_path = args.slugs_path + "/sortmerna-2.0/rRNA_databases/"
databases = args.databases.split(',')
tmp = []

for d in databases:
    tmp.append(database_path + d + ',' + args.root + "/index/" + d.split('.')[0])

output = ":".join(tmp)

print output
