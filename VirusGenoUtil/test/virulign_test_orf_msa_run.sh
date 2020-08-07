#!/bin/bash
source /etc/profile

EXP_PATH=$HOME/Documents/L3S/projects/sars_cov2/data
OUT_PATH=$EXP_PATH/test_alignments
REMOVED_PATH=$OUT_PATH/removed_genomes
mkdir $OUT_PATH
mkdir $REMOVED_PATH
virulign $EXP_PATH/test_orf_xml/test_orf.xml $EXP_PATH/test_orf_sequences/test_sequences.fasta --exportKind GlobalAlignment --exportWithInsertions no --exportAlphabet Nucleotides --nt-debug $REMOVED_PATH > $OUT_PATH/test_mfa_nt.fasta 2> $OUT_PATH/test_mfa_nt.err
