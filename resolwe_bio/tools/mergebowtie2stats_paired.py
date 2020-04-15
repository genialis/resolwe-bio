#!/usr/bin/env python3
"""Merge Bowtie2 statistics from paired-end reads."""

import sys

if len(sys.argv) < 2:
    sys.stderr.write("No stats file given.\n")
    exit(1)

stats = []
t_unique, t_multiple, t_discord = 0, 0, 0

with open(sys.argv[1]) as f:
    iteration, processed, not_aligned, unique, multiple, discord, mapped = (
        "Initial alignment",
        0,
        0,
        0,
        0,
        0,
        0,
    )
    for line in f:
        vals = line.strip().split(" ")

        if "reads; of these" in line:
            processed = int(vals[0])

        if "pairs aligned 0 times" in line:
            not_aligned = int(vals[0])

        if "aligned concordantly exactly 1 time" in line:
            unique = int(vals[0])
            t_unique += unique

        if "aligned concordantly >1 times" in line:
            multiple = int(vals[0])
            t_multiple += multiple

        if "aligned discordantly 1 time" in line:
            discord = int(vals[0])
            t_discord += discord

        if "Trimming iteration" in line:
            iteration = line.strip()

        if "overall alignment rate" in line:
            mapped = float(vals[0].strip("%"))
            stats.append(
                (iteration, processed, not_aligned, unique, multiple, discord, mapped)
            )
            processed, not_aligned, unique, multiple, discord, mapped = (
                0,
                0,
                0,
                0,
                0,
                0,
            )

with open("stats.tab", "w") as f:
    f.write(
        "Alignment\tPairs processed\tPairs aligned 0 times\t"
        "Pairs aligned exactly 1 time\tPairs aligned >1 times\t"
        "Pairs aligned discordantly\tOverall alignment rate (%)\n"
    )
    for vals in stats:
        f.write("\t".join(map(str, vals)) + "\n")

    t_mapped = round(((t_unique + t_multiple + t_discord) / stats[0][1]) * 100, 1)

    f.write(
        "\t".join(
            map(
                str,
                (
                    "Total",
                    stats[0][1],
                    stats[-1][2],
                    t_unique,
                    t_multiple,
                    t_discord,
                    t_mapped,
                ),
            )
        )
        + "\n"
    )
