#!/bin/bash

echo "DeprecationWarning: This function will be removed in near future. Use gen-import instead!"

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
# MAX_PROGRES:   maximum progress at the end of transfer (0.0 to 1.0)         #
# COMPRESSION:   if "compress" return compressed data                         #
#                if "extract" return extracted data                           #
#                else return both                                             #
#                                                                             #
###############################################################################
TEMP=$1
FILE=$2
IN_FORMAT=$3
OUT_FORMAT=$4
MAX_PROGRES=${5:-1.0}
COMPRESSION=$6
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

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

  if [ "$COMPRESSION" != extract ]; then
    gzip -t "${NAME}.${OUT_FORMAT}.gz"
    testrc 1  # RC 2 "trailing garbage ignored" is OK
  fi

  if [ "$COMPRESSION" != compress ]; then
    gzip -dc "${NAME}.${OUT_FORMAT}.gz" > "${NAME}.${OUT_FORMAT}"
    testrc 1  # RC 2 "trailing garbage ignored" is OK

    if [ "$COMPRESSION" = extract ]; then
      rm "${NAME}.${OUT_FORMAT}.gz"
      testrc
    fi
  fi
}

function import7z {
  # Uncompress original file
  if [[ ".${FILE}" =~ \.(tgz|tar\.gz|tar\.bz2)$ ]]; then
    7z x -y -so "${TEMP}" | tar -xO > "${NAME}.${OUT_FORMAT}"
    testrc
  else
    7z x -y -so "${TEMP}" > "${NAME}.${OUT_FORMAT}"
    testrc
  fi

  # Remove original file
  rm "${TEMP}"
  testrc

  if [ "$COMPRESSION" != extract ]; then
    # Compress uncompressed file
    gzip -c "${NAME}.${OUT_FORMAT}" > "${NAME}.${OUT_FORMAT}.gz"
    testrc

    if [ "$COMPRESSION" = compress ]; then
      # Remove uncompressed file
      rm "${NAME}.${OUT_FORMAT}"
      testrc
    fi
  fi
}

function importUncompressed {
  if [ "$COMPRESSION" = compress ]; then
    gzip -c "${TEMP}" > "${NAME}.${OUT_FORMAT}.gz"
    testrc
    rm "${TEMP}"
    testrc
  elif [ "$COMPRESSION" = extract ]; then
    mv "${TEMP}" "${NAME}.${OUT_FORMAT}"
    testrc
  else
    gzip -c "${TEMP}" > "${NAME}.${OUT_FORMAT}.gz"
    testrc
    mv "${TEMP}" "${NAME}.${OUT_FORMAT}"
    testrc
  fi
}

regex='^(https?|ftp)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'
if [[ "$TEMP" =~ $regex ]]
then
    URL=${TEMP}
    FILE=${FILE:-`basename "${URL%%\?*}"`}
    TEMP=download_`basename "${URL%%\?*}"`
    curl --connect-timeout 10 -a --retry 10 -# -L -o "${TEMP}" "${URL}" 2>&1 | stdbuf -oL tr '\r' '\n' | grep -o '[0-9]*\.[0-9]' | python -u $DIR/curlprogress.py --scale $MAX_PROGRES
    testrc
fi

# Check if a temporary file exists
if [ ! -f "${TEMP}" ]; then
    echo "{\"proc.error\":\"File transfer failed: temporary file not found\"}"
fi

# Set FILE to extracted filename from TEMP if FILE not set
FILE=${FILE:-`basename "${TEMP%%\?*}"`}

# Take basename if FILE not nice
FILE=`basename "${FILE%%\?*}"`

# Add a dot to all input formats except for no extension
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

echo "{\"proc.progress\":$MAX_PROGRES}"
