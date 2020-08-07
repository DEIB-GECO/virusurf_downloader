#!/bin/bash
source /etc/profile

#=============================================================================#
# Task: BLAST B-cells immunodominant regions from sars cov1 to sars cov2      #
#=============================================================================#

#=============================================================================#
# input folder: general output folder=                                        #
# /home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Bcells           #
# output folder: blast folder =                                               #
# /home/damian/Documents/L3S/projects/sars_cov2/data/exp_epitopes_input/blast #
# Note: the output folder should be created prior to blast run	              #
#=============================================================================#

# set up input and output folders

input_dir="/home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Bcells"
output_dir="/home/damian/Documents/L3S/projects/sars_cov2/sars_cov2_data/exp_epitopes_input/blast"

echo "BLASTP immunodominant regions of P59594 (Spike) to sars cov2 structural proteins"
blastp -query $input_dir"/immunodom_regions_P59594.fasta" -db $output_dir"/sars_cov2_SMN" -max_target_seqs 2 -max_hsps 1 -evalue 1e-6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident evalue' -out $output_dir"/P59594_blastp_sars_cov2.tsv" -num_threads 2
sed -i '1 i\qseqid\tsseqid\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tpident\tevalue' $output_dir"/P59594_blastp_sars_cov2.tsv"

echo "BLASTP immunodominant regions of P59595 (Nucleoprotein) to sars cov2 structural proteins"
blastp -query $input_dir"/immunodom_regions_P59595.fasta" -db $output_dir"/sars_cov2_SMN" -max_target_seqs 2 -max_hsps 1 -evalue 1e-6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident evalue' -out $output_dir"/P59595_blastp_sars_cov2.tsv" -num_threads 2
sed -i '1 i\qseqid\tsseqid\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tpident\tevalue' $output_dir"/P59595_blastp_sars_cov2.tsv"

echo "BLASTP immunodominant regions of P59596 (Membrane) to sars cov2 structural proteins"
blastp -query $input_dir"/immunodom_regions_P59596.fasta" -db $output_dir"/sars_cov2_SMN" -max_target_seqs 2 -max_hsps 1 -evalue 1e-6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident evalue' -out $output_dir"/P59596_blastp_sars_cov2.tsv" -num_threads 2
sed -i '1 i\qseqid\tsseqid\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tpident\tevalue' $output_dir"/P59596_blastp_sars_cov2.tsv"
