#!/bin/bash
source /etc/profile

#=====================================================================#
# Task: Map target sequence to reference for each ORF using aminoacid #
# sequence, all parameters in default                                 #
#=====================================================================#

#=====================================================================#
# in data folder containing a folder with the target sequences fasta  #
# and assuming virulign is installed in selected path,                #
# run the following commands                                          #
#=====================================================================#

# set up the bin path
virulign_bin="/home/melidis/bin/virulign/bin/virulign"
# set up data folder
data_path="/home/melidis/sars_cov2/data"

# follow the ORF order: ORF1ab,S,ORF3a,E,M,ORF6,ORF7a,ORF7b,ORF8,N,ORF10
# from NCBI annotation track: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=graph

echo "Process ORF1ab"
mkdir ORF1ab_removed
$virulign_bin $data_path"/xmls/ORF1ab.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF1ab_removed > ORF1ab_mfa_aa.fasta 2> ORF1ab_mfa_aa.err

echo "Process S"
mkdir S_removed
$virulign_bin $data_path"/xmls/S.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug S_removed > S_mfa_aa.fasta 2> S_mfa_aa.err

echo "Process ORF3a"
mkdir ORF3a_removed
$virulign_bin $data_path"/xmls/ORF3a.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF3a_removed > ORF3a_mfa_aa.fasta 2> ORF3a_mfa_aa.err

echo "Process E"
mkdir E_removed
$virulign_bin $data_path"/xmls/E.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug E_removed > E_mfa_aa.fasta 2> E_mfa_aa.err

echo "Process M"
mkdir M_removed
$virulign_bin $data_path"/xmls/M.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug M_removed > M_mfa_aa.fasta 2> M_mfa_aa.err

echo "Process ORF6"
mkdir ORF6_removed
$virulign_bin $data_path"/xmls/ORF6.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF6_removed > ORF6_mfa_aa.fasta 2> ORF6_mfa_aa.err

echo "Process ORF7a"
mkdir ORF7a_removed
$virulign_bin $data_path"/xmls/ORF7a.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF7a_removed > ORF7a_mfa_aa.fasta 2> ORF7a_mfa_aa.err

echo "Process ORF7b"
mkdir ORF7b_removed
$virulign_bin $data_path"/xmls/ORF7b.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF7b_removed > ORF7b_mfa_aa.fasta 2> ORF7b_mfa_aa.err

echo "Process ORF8"
mkdir ORF8_removed
$virulign_bin $data_path"/xmls/ORF8.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF8_removed > ORF8_mfa_aa.fasta 2> ORF8_mfa_aa.err

echo "Process N"
mkdir N_removed
$virulign_bin $data_path"/xmls/N.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug N_removed > N_mfa_aa.fasta 2> N_mfa_aa.err

echo "Process ORF10"
mkdir ORF10_removed
$virulign_bin $data_path"/xmls/ORF10.xml" $data_path"/targets/sequences.fasta" --exportKind GlobalAlignment --exportAlphabet AminoAcids --exportReferenceSequence yes --nt-debug ORF10_removed > ORF10_mfa_aa.fasta 2> ORF10_mfa_aa.err
