#! /bin/bash
"""
Trim UMI-sequences from single or paired-end reads.

The script appends the UMI-tag sequence of length 'match_len' to the end of read names.
Sequence of length 'trim_size' is trimmed from the start of each mate1 read.
In the case of paired-end reads, the script checks the end of mate2 reads for the presence
of UMI-tags (in reverse complement), and trims the reads accordingly.

Awk commands adapted after Corall analysis protocol script developed by Lexogen GmbH (July 2019).
"""

source /etc/profile.d/resolwe-base.sh


if [[ $# -lt 3 || $# -gt 4 ]]; then
  re-error "Usage: extract_umi.sh match_len trim_size mate1.fastq [mate2.fastq]."
fi

match_len="$1"
trim_size="$2"
mate1="$3"
mate2="$4"

outName=$(basename "${mate1}" .fastq.gz)

if [ "${mate2}" == "" ]
then
    echo "Extracting UMI-tags from single-end input file ${mate1}."
    awk -v match_len="${match_len}" -v trim_size="${trim_size}" \
    'function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
    NR%4==1{ rd_name=$1; rd_info=$2 }
    NR%4==2{ umi=substr($1,1,match_len); rd_seq=substr($1,trim_size) }
    NR%4==0{ id_line=rd_name"_"umi" "rd_info;print rtrim(id_line); print rd_seq; print "+"; print substr($1,13) }' <(zcat "${mate1}") | pigz > "${outName}_umi.fastq.gz"
else
    echo "Extracting UMI-tags from paired-end input files ${mate1} and ${mate2}."

	outName2=$(basename "${mate2}" .fastq.gz)

    # Get a number of lines in mate1 file
	mate1_nlines=$(wc -l < <(zcat "${mate1}"))

	# check the defined window size at end of read2 for UMI trimming
	regExWin=24

    # requires gawk >= 3.1.3 built with "--enable-switch" option
    # Splits the awk output stream into mate1/mate2 files, based on the number of lines in the input mate1 file
	split -l $mate1_nlines --numeric-suffixes=1 --additional-suffix=".fastq.gz" --filter='pigz -c - > $FILE' <( \
    awk -v win="${regExWin}" -v match_len="${match_len}" -v trim_size="${trim_size}" \
    'function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
		function RC_match(match1,seq){
			m1=RC_fct(match1);
				return index(seq,m1);
		}
		function RC_fct(umi){
			split(umi,str,"")
			RC="";
			for(i=0;i<=length(umi);++i) {
				switch (str[i]) {
					case "A":
						RC="T"RC;
						break
					case "T":
						RC="A"RC;
						break
					case "G":
						RC="C"RC;
						break
					case "C":
						RC="G"RC;
						break
					default:
						break
				}
			}
			return RC;
		}

		BEGIN{count=0;winM=win-1;winP=win+1-match_len}
		{if(NR==FNR){
		if(FNR%4==1){ rd_name=$1; rd_info=$2 }
		if(FNR%4==2){ umi=substr($1,1,match_len); matchRd=substr($1,trim_size,match_len);rd_seq=substr($1,trim_size); rd_name2umi[rd_name]=umi;matchPt[rd_name]=matchRd;}
		if(FNR%4==0){ id_line=rd_name"_"umi" "rd_info; print rtrim(id_line); print rd_seq; print "+"; print substr($1,trim_size) }
		next }
		{
		if(FNR%4==1){ if(!rd_name2umi[$1]){ print "ERROR: no read1 for read2="$1 > "/dev/stderr"; exit 1 }else{ umi=rd_name2umi[$1];matchRd2=matchPt[$1]; id_line=$1"_"umi" "$2; print rtrim(id_line)} }
		if(FNR%4==2){val=RC_match(matchRd2,substr($1,length($1)-winM)); if(val!=0){print substr($1,1,length($1)-winP+val);count++}else {print}}
		if(FNR%4==3){print}
		if(FNR%4==0){if(val!=0){print substr($1,1,length($1)-winP+val)}else {print}}
		}}END{print "removed UMI seqences in rd2: " count > "/dev/stderr"}' <(zcat "${mate1}") <(zcat "${mate2}")) "mate_"

    # Compress the paired-end output, if it exists
    if [[ -s "mate_01.fastq.gz" && -s "mate_02.fastq.gz" ]]
    then
        mv "mate_01.fastq.gz" "${outName}_umi.fastq.gz"
        mv "mate_02.fastq.gz" "${outName2}_umi.fastq.gz"
    else
      re-error "Failed to extract UMI-sequences from paired-end input. Output .fastq.gz files do not exist or are empty."
    fi
fi
