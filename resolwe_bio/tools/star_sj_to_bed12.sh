#! /bin/bash

source /etc/profile.d/resolwe-base.sh

star_sj_file="$1"

# recalculate from STAR SJ.out.tab file to BED12 format
# BED12 consists of 12 columns:
# 1: chromosome - $1: copy first column from STAR SJ.out.tab
# 2: SJ start (0-based) - `$2-$9-1`: (first base of the intron (1-based) - maximum spliced alignment overhang) -1 (to recalculate to 0 based system)
# 3: SJ end (0-based) - `$3+$9`: (last base of the intron (1-based) + maximum spliced alignment overhang) (to recalculate to 0 based system)
# 4: SJ name - `JUNC000"NR"`
# 5: score - `$7`: number of uniquely mapping reads crossing the junction
# 6: strand - depends on condition:
#               `if($4=="2")`, strand is '-'
#               `if($4=="1")`, strand is '+'
# 7: thick start is the same as SJ start : $2-$9-1
# 8: thick end is the same as SJ end : $3+$9
# 9: item RGB - `255,0,0`
# 10: block counts - `2`
# 11: block sizes - ``"$9","$9"`: maximum spliced alignment overhang, maximum spliced alignment overhang
# 12: block starts - `0,"$3-$2+$9+1`: 0, (SJ end - SJ start + maximum spliced alignment overhang -1 )

awk \
    {'if($4=="2") print \
""$1"\t\
"$2-$9-1"\t\
"$3+$9"\t\
JUNC000"NR"\t\
"$7"\t\
-\t\
"$2-$9-1"\t\
"$3+$9"\t\
255,0,0\t\
2\t\
"$9","$9"\t\
0,"$3-$2+$9+1""; \
    else \
    if($4=="1") print \
""$1"\t\
"$2-$9-1"\t\
"$3+$9"\t\
JUNC000"NR"\t\
"$7"\t\
+\t\
"$2-$9-1"\t\
"$3+$9"\t\
0,0,255\t\
2\t\
"$9","$9"\t\
0,"$3-$2+$9+1""'} \
    "${star_sj_file}" > junctions_unsorted.bed
re-checkrc "Calculation from STAR SJ.out.tab file to bed file failed."
