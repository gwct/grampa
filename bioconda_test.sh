#!/bin/bash
set -e
set -o pipefail
set -x

# Runs test for PhyloAcc, including on a small simulated dataset that contains a fasta file, mod file,
# bed file, id subset file, and config file.

TMP=$(mktemp -d)
trap 'rm -rf $TMP' EXIT
export TMPDIR=$TMP
cd $TMP

echo " ** DOWNLOADING TEST DATA."
files=(
  "bioconda-test-cfg.yaml"
  "id-subset.txt"
  "ratite.mod"
  "simu_500_200_diffr_2-1.bed"
  "simu_500_200_diffr_2-1.noanc.fa"
)

for file in "${files[@]}"; do
  if ! wget -q "https://github.com/phyloacc/PhyloAcc-test-data/raw/main/bioconda-test-data/$file"; then
    echo "Failed to download $file" >&2
    exit 1
  fi
done
echo " ** TEST DATA DOWNLOAD OK."

# echo " ** BEGIN BINARY TEST."