#! /bin/bash

FASTA=$1
OUTPUT=$2

echo Python >> $OUTPUT
/usr/bin/time -v python python/computePairwiseDistances.py $1 /dev/null 2>> $OUTPUT

echo OCaml >> $OUTPUT
/usr/bin/time -v ./ocaml/computePairwiseDistances --fasta $1 --output /dev/null 2>> $OUTPUT
