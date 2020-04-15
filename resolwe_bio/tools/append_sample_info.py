#!/usr/bin/env python3
"""Extract reference (FASTA) and sample names from the VCF file."""

import argparse
import os

from pysam import VariantFile
from resolwe_runtime_utils import error, warning

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("vcf_file", help="VCF file (can be compressed using gzip/bgzip).")
parser.add_argument("summary", help="Summary file to append to.")
args = parser.parse_args()

try:
    vcf = VariantFile(args.vcf_file)
except (OSError, ValueError) as error_msg:
    proc_error = "Input VCF file does not exist or could not be correctly opened."
    print(error(proc_error))
    raise ValueError(error_msg)

vcf_header = vcf.header
header_records = {record.key: record.value for record in vcf_header.records}

with open(args.summary, "a") as out_file:
    try:
        fasta_name = os.path.basename(header_records["reference"])
    except KeyError:
        fasta_name = ""
        print(
            warning(
                "Reference sequence (FASTA) name could not be recognized from the VCF header."
            )
        )

    out_file.write("\nReference (genome) sequence:\n{}\n".format(fasta_name))
    out_file.write("\nSamples:\n{}".format("\n".join(list(vcf_header.samples))))
