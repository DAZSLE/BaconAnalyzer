#!/usr/bin/env bash

# This does the same thing as makeList.sh but slightly faster

eoscmd="eos root://cmseos.fnal.gov/"
prefix="/store/group/lpcbacon/"

if [[ $@ < 1 ]]; then
  echo "Usage: $0 [BaconProd number]"
  echo "  Assumes production directory: $prefix"
  exit 1
fi

prod=$1

if [[ -d production${prod} ]]; then
  echo "prodution${prod} exists already!"
  exit 1
fi

mkdir -p production${prod}

echo "Scanning $prefix/${prod}..."
$eoscmd ls $prefix/$prod | while read dir; do
  fout=production${prod}/${dir}.txt
  $eoscmd find --xurl -name \*.root $prefix/$prod/$dir | grep -v failed > $fout
  if [[ ! -s $fout ]]; then
    echo "Warning: $dir has no root files in it!"
    rm $fout
  fi
done
echo "Done"
