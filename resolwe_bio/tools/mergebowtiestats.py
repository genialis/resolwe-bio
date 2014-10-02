import os
import sys

if len(sys.argv) < 2:
    sys.stderr.write('No stats file given.')
    print
    exit(1)

for i in range(2, len(sys.argv), 2):
    if not os.path.isfile(sys.argv[i]):
        sys.stderr.write('Stats file {} not found.'.format(sys.argv[i]))
        print
        exit(1)

stats = []
tprocessed, tonevalid, tfailed, tsuppressed = 0, 0, 0, 0

for i in range(1, len(sys.argv), 2):
    trimmed = sys.argv[i]
    with open(sys.argv[i + 1]) as f:
        processed, onevalid, failed, suppressed = -1, -1, -1, -1

        for line in f:
            vals = line.strip().split(' ')

            if 'reads processed' in line:
                processed = int(vals[-1])

            if 'reads with at least one reported alignment' in line:
                onevalid = int(vals[-2])

            if 'reads that failed to align' in line:
                failed = int(vals[-2])

            if 'reads with alignments suppressed due to -m' in line:
                suppressed = int(vals[-2])

        if processed > 0:
            tprocessed += processed
        if onevalid > 0:
            tonevalid += onevalid
        if failed > 0:
            tfailed += failed
        if suppressed > 0:
            tsuppressed += suppressed

        stats.append((trimmed, processed, onevalid, failed, suppressed))

with open('stats.tab', 'w') as f:
    f.write("Trim3 size\tReads processed\tReads with at least one reported alignment\t"
            "Reads that failed to align\tReads with alignments suppressed due to\n")
    for vals in stats:
        f.write("\t".join(map(str, vals)) + "\n")

    f.write("\t".join(map(str, ("Total", tprocessed, tonevalid, tfailed, tsuppressed))) + "\n")
