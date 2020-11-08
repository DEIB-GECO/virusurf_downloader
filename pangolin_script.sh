source ~/miniconda3/etc/profile.d/conda.sh
conda activate pangolin
echo "USAGE: bash pangoline_script.sh input.fasta outputDir outputFile"
pangolin $1 -o $2 --outfile $3 
