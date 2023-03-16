#!/bin/bash

python ./scripts/primer_blast_local.py -g ./test/genome_test.fasta \
  -p ./test/primers_test.fasta \
  -o ./test/output/test0 \
  -m 50.

python ./scripts/primer_blast_local.py -g ./test/genome_test.fasta \
  -p ./test/IDT_PrimerQuest_Export.xls \
  -o ./test/output/test1 \
  -m 50.
