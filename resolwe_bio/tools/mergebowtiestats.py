#!/usr/bin/env python3
"""Merge Bowtie statistics."""

import os
import sys

if len(sys.argv) < 2:
    sys.stderr.write("No stats file given.\n")
    exit(1)

for i in range(2, len(sys.argv), 2):
    if not os.path.isfile(sys.argv[i]):
        sys.stderr.write("Stats file {} not found.\n".format(sys.argv[i]))
        exit(1)

stats = []
tprocessed, tonevalid, tfailed, tsuppressed = 0, 0, 0, 0

for i in range(1, len(sys.argv), 2):
    trimmed = sys.argv[i]
    with open(sys.argv[i + 1]) as f:
        processed, onevalid, failed, suppressed = -1, -1, -1, -1

        for line in f:
            vals = line.strip().split(" ")

            if "reads processed" in line:
                processed = int(vals[-1])

            if "reads with at least one reported alignment" in line:
                onevalid = int(vals[-2])

            if "reads that failed to align" in line:
                failed = int(vals[-2])

            if "reads with alignments suppressed due to -m" in line:
                suppressed = int(vals[-2])

        if onevalid > 0:
            tonevalid += onevalid

        mapped = round((onevalid / processed) * 100, 1)

        stats.append((trimmed, processed, onevalid, failed, suppressed, mapped))

with open("stats.tab", "w") as f:
    f.write(
        "Trim3 size\tReads processed\t"
        "Reads with at least one reported alignment\t"
        "Reads that failed to align\t"
        "Reads with alignments suppressed due to -m\tMapped (%)\n"
    )
    for vals in stats:
        f.write("\t".join(map(str, vals)) + "\n")

    tmapped = round((tonevalid / stats[0][1]) * 100, 1)

    f.write(
        "\t".join(
            map(
                str,
                ("Total", stats[0][1], tonevalid, stats[-1][3], stats[-1][4], tmapped),
            )
        )
        + "\n"
    )
