#!/bin/bash

# Adapted from https://github.com/learnedsystems/SOSD/blob/master/scripts/download.sh
# and https://gist.github.com/fkraeutli/66fa741d9a8c2a6a238a01d17ed0edc5

function main() {
  mkdir -p SOSD

  echo "Compiling Binaries..."
  make clean && make all

  echo "Downloading SOSD Datasets..."
  cd SOSD
  download_file_zst books_800M_uint64 8708eb3e1757640ba18dcd3a0dbb53bc https://www.dropbox.com/s/y2u3nbanbnbmg7n/books_800M_uint64.zst?dl=1
  download_file_zst fb_200M_uint64 3b0f820caa0d62150e87ce94ec989978 https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/EATHF7

}

# Calculate md5 checksum of FILE and stores it in MD5_RESULT
function get_checksum() {
  FILE=$1

  if [ -x "$(command -v md5sum)" ]; then
    # Linux
    MD5_RESULT=$(md5sum ${FILE} | awk '{ print $1 }')
  else
    # OS X
    MD5_RESULT=$(md5 -q ${FILE})
  fi
}

function download_file_zst() {
  FILE=$1
  CHECKSUM=$2
  URL=$3

  # Check if file already exists
  if [ -f ${FILE} ]; then
    # Exists -> check the checksum
    get_checksum ${FILE}
    if [ "${MD5_RESULT}" != "${CHECKSUM}" ]; then
      wget -O - ${URL} | zstd -d >${FILE}
    fi
  else
    # Does not exists -> download
    wget -O - ${URL} | zstd -d >${FILE}
  fi

  # Validate (at this point the file should really exist)
  get_checksum ${FILE}
  if [ "${MD5_RESULT}" != "${CHECKSUM}" ]; then
    echo "error checksum does not match: run download again"
    exit -1
  else
    echo ${FILE} "checksum ok"
  fi
}

# Run
main
