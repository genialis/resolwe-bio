#!/bin/bash

###############################################################################
#                                                                             #
# INPUT ARGUMENTS                                                             #
#                                                                             #
# TEMP: source file location                                                  #
# FILE: real file name                                                        #
# IN_FORMAT:                                                                  #
#   fa|fasta     1.) matches files that end with fa or fasta                  #
#                2.) matches a combination of:                                #
#                    (fa|fasta).(gz|bz2|zip|rar|7z|tgz|tar.gz|tar.bz2)        #
#   fa|fasta|zip 1.) matches files that end with fa or fasta                  #
#                2.) matches a combination of:                                #
#                    (fa|fasta).(gz|bz2|zip|rar|7z|tgz|tar.gz|tar.bz2)        #
#                3.) matches if the file end just with zip                    #
#                    supported (gz|bz2|zip|rar|7z)                            #
# OUT_FORMAT:    the desired output format before compression, ie. fasta      #
# LEAVE_UNCOMP:  if any string given (uncomp string is prefered),             #
#                leave the uncompressed file as well                          #
#                                                                             #
###############################################################################

TEMP=$1
FILE=$2
IN_FORMAT=$3
OUT_FORMAT=$4
LEAVE_UNCOMP=$5

echo "Importing and compressing..."
shopt -s nocasematch

function testrc {
  RC=$?
  if [ $1 ]; then
    if [ $RC -eq $1 ]; then
      echo "{\"proc.rc\":$RC}"
      exit $RC
    fi
  else
    if [ $RC -gt 0 ]; then
      echo "{\"proc.rc\":$RC}"
      exit $RC
    fi
  fi
}

function importGz {
  mv "${TEMP}" "${NAME}.${OUT_FORMAT}.gz"
  testrc

  # Return code 2 "trailing garbage ignored" is OK
  gzip -t "${NAME}.${OUT_FORMAT}.gz"
  testrc 1

  if [[ $LEAVE_UNCOMP ]]; then
    # Return code 2 "trailing garbage ignored" is OK
    gzip -dc "${NAME}.${OUT_FORMAT}.gz" > "${NAME}.${OUT_FORMAT}"
    testrc 1
  fi
}

function import7z {
  if [[ ".${FILE}" =~ \.(tgz|tar\.gz|tar\.bz2)$ ]]; then
    7z x -y -so "${TEMP}" | tar -xO > "${NAME}.${OUT_FORMAT}"
    testrc
  else
    7z x -y -so "${TEMP}" > "${NAME}.${OUT_FORMAT}"
    testrc
  fi

  rm "${TEMP}"
  testrc

  gzip -c "${NAME}.${OUT_FORMAT}" > "${NAME}.${OUT_FORMAT}.gz"
  testrc

  if [[ ! $LEAVE_UNCOMP ]]; then
    rm "${NAME}.${OUT_FORMAT}"
    testrc
  fi
}

function importUncompressed {
  gzip -c "${TEMP}" > "${NAME}.${OUT_FORMAT}.gz"
  testrc

  if [[ $LEAVE_UNCOMP ]]; then
    mv "${TEMP}" "${NAME}.${OUT_FORMAT}"
    testrc
  else
    rm "${TEMP}"
    testrc
  fi
}

# add a dot to all input formats except for no extension
# txt|csv -> .txt|.csv, txt|csv| -> .txt|.csv|
IN_FORMAT=`python -c "print '|'.join(['.' + a if a else a for a in '$IN_FORMAT'.split('|')])"`

# Decide which import to use based on the $FILE extension and the $IN_FORMAT
if [[ ".${FILE}" =~ (${IN_FORMAT})\.gz$ ]]; then
  export NAME=`echo "$FILE" | sed -E "s/(${IN_FORMAT})\.gz$//g"`
  importGz

elif [[ ".${FILE}" =~ (${IN_FORMAT})\.(bz2|zip|rar|7z|tgz|tar\.gz|tar\.bz2)$ ]]; then
  export NAME=`echo "$FILE" | sed -E "s/(${IN_FORMAT})\.(bz2|zip|rar|7z|tgz|tar\.gz|tar\.bz2)$//g"`
  import7z

elif [[ ".${FILE}" =~ (${IN_FORMAT})$ ]]; then

  if [[ ".${FILE}" =~ \.gz$ ]]; then
    export NAME=`echo "$FILE" | sed -E "s/\.gz$//g"`
    importGz

  elif [[ ".${FILE}" =~ \.(bz2|zip|rar|7z|tgz|tar\.gz|tar\.bz2)$ ]]; then
    export NAME=`echo "$FILE" | sed -E "s/\.(bz2|zip|rar|7z|tgz|tar\.gz|tar\.bz2)$//g"`
    import7z

  else
    export NAME=`echo "$FILE" | sed -E "s/(${IN_FORMAT})$//g"`
    importUncompressed

  fi
else
  echo "{\"proc.rc\":1}"
  exit 1
fi
