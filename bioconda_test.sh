#!/bin/bash
set -e
set -o pipefail
set -x

# Runs test for PhyloAcc, including on a small simulated dataset that contains a fasta file, mod file,
# bed file, id subset file, and config file.

# TMP=$(mktemp -d)
# trap 'rm -rf $TMP' EXIT
# export TMPDIR=$TMP
# cd $TMP

echo " ** DOWNLOADING TEST DATA."

if ! wget -q "https://github.com/gwct/grampa/raw/main/bioconda-test-data.zip"; then
    echo "Failed to download $file" >&2
exit 1
fi
unzip bioconda-test-data.zip

echo " ** TEST DATA DOWNLOAD OK."

echo " ** BEGIN GRAMPA TEST."
if ! python grampa.py --tests; then
  echo " ** ERROR: GRAMPA tests failed." >&2
  exit 1
fi
echo " ** GRAMPA TEST OK."