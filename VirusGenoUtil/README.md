# VirusGenoUtil
This is a project to contribute utilities for a viral genomic data explorer.

## Dependencies
* Python 3.7.6
* BioPython 1.74
* Pandas 1.0.1
* NumPy 1.19.1
* Requests 2.24.0
## Tested functionalities
1. Convert nucleic-acid sequence ORF multiple sequence alignemnts, provived by virulign aligner, saved as multi-fasta alignments (MFA) to variants per target sequence

2. Convert amino-acid sequence ORF multiple sequence alignemnts, provived by virulign aligner, saved as multi-fasta alignments (MFA) to variants per target sequence

3. Identify homology based immunodominant regions for structural proteins of Sars-Cov-2 following ["Grifoni et al."](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(20)30166-9)

4. Extract IEDB epitopes from IEDB data portal

## Convert virulign MFA to variants csv - How to use
a. install virulign (on a server if you need to run for many sequences)

b. depending on the sequence **type**:
* for amino-acid sequence:
  check the example [virulign_msa_aa_run.sh](bash/virulign_msa_aa_run.sh) used for running virulign for all ncbi sars cov2 as of May 15, 
* change this bash appropriately and run it for your own paths

c. put all output MFA fasta (one per ORF) to an alignments folder

d. go to [main.py](code/main.py) and change appropriately for the input and output paths,
* **input** should be the path of such alignment folder
* please change the refernce ncbi id
* create a file that contains your email (will be needed for programmatically accessing the Entrez service and fetch the genbank file of the reference sequence)
* **output** will be a folder containing a variant csv file for all ORFs of one target sequence, for example see [output](test/output/test_variants/)

## Identify immunodominant regions using homology - How to use
a. download the spike, membrane and nucleoprotein proteins of SARS-CoV(NC\_004718.3):
* create input root:
  mkdir data/exp\_epitopes\_input
* create SARS-CoV proteins folder:
  mkdir data/exp\_epitopes\_input/sars\_cov1\_proteins
* save each protein on a separate protein\_id fasta file

b. download the same proteins for SARS-CoV2(NC\_045512.2):
* create SARS-CoV2 proteins folder:
  mkdir data/exp\_epitopes\_input/sars\_cov2\_proteins
* save all in sars\_cov2\_proteins.fasta

c. prepare BLAST db for SARS-CoV2 structural proteins:
* create blast folder:
  mkdir data/exp\_epitopes\_input/blast
* concatenate all SARS-CoV2 proteins in one fasta file:
  cat data/exp\_epitopes\_input/sars_cov2\_proteins/*.fasta > data/exp\_epitopes\_input/blast/sars\_cov2\_proteins.fasta
* build BLAST db:
  makeblastdb -in sars\_cov2\_proteins.fasta -out sars\_cov2\_SMN -dbtype prot -title "Sars Cov2 Spike Membrane Nucleoprotein" -parse_seqids

d. prepare response frequency data:
* get response frequency assays for B-cells and T-cells from ["IEDB"](https://www.iedb.org/home_v3.php)
* save to data/exp\_epitopes\_input/Bcells and data/exp\_epitopes\_inpyt/Tcells respectively

e. go to [exp_epitopes_run.py](code/exp_epitopes_run.py)
* prepare and run the Immunodominance class
* in the specified output folder, the immunodominant regions for SARS-CoV will be saved and optionally the sliding average RF plots per protein

f. BLASTP these regions to SARS-CoV2 proteins
* see folder [epitopes_homology_based](epitopes\_homology\_based) for example blast bash scripts

g. go to [exp_epitopes_run.py](code/exp_epitopes_run.py)
* prepare and run the HomologyBasedEpitopes class
* in the specified output folder, the homology based identified epitopes tsv will be saved

## Extract epitopes from prediction tools - How to use

### Bepipred 2.0
1. concatenate structural proteins into one fasta file (time-efficient way to use bepipred predictor):
cd sars\_cov2\_data/pred\_epitopes\_input/sars\_cov2\_proteins
cat *.fasta > sars_cov2_struct_prot.fasta

2. run Bepipred 2.0 for the SARS-CoV-2 structural proteins:
cd sars\_cov2/sars\_cov2\_data/pred\_epitopes\_input
BepiPred-2.0 -t 0.55 sars\_cov2\_proteins/sars\_cov2\_struct\_prot.fasta > bepipred\_out/bepipred\_sars\_cov2.tsv

3. set parameters on [pred_epitopes_run.py](code/pred_epitopes.py) and run section for BepipredEpitopes:
 * the output tsv file will be placed on your specified output folder:

   ls -lh sars\_cov2/pred\_epitopes/bepipred/bepipred\_sars\_cov2\_epi.tsv

## Extract epitopes from IEDB
IEDB data
* go to iedb database export ["page"](http://www.iedb.org/database_export_v3.php)
* download the B and T cells and the mhc_ligan_full(single file) csv.zip files.
* decompress each of them and then compress it again using gzip. 
* finally, place the three gzip files in [iedb_data](data/iedb_input/cell_epitopes) folder.

Viruses data
* create a folder containing the viruses of interest. Save each virus as a separate subfolder and name as taxon_ncbi_id, where ncbi id is the taxonomic ncbi id of this virus. For example for Sars-CoV-2, the subfolder should be named taxon_2697049.
* in each subfolder add each viral protein as a separate .fa file. Then name the protein file using the uniprot id of this protein. For example, for Spike protein name it as P59594.fa. For more information on how to create a viral folder, consult an example download [log](data/iedb_input/viruses_proteins/taxon_2697049/how2download.log) file.
* finally, you have created a viruses_proteins folder as the example [one](data/iedb_input/viruses_proteins).
 
ONTIE downloads
 * re-use the already downloaded ONTIE .ttl files from ONTIE ["page"](https://ontology.iedb.org/ontology), through already executions of the code.
 
Email for Entrez services
 * please add your email at the [email.txt](data/iedb_input/cell_epitopes/email.txt). It will be needed to run Entrez services through bio-python.
 
Extracting IEDB epitopes

a. extract epitopes, save them in csv files
 * set up the paths for your machine in [iedb_epitopes_run.py](code/iedb_epitopes_run.py)
 * run [iedb_epitopes_run.py](code/iedb_epitopes_run.py)
 * in specified output folder two csv files will be created. The imported_iedb_epitopes.tsv is the realization of the Epitope table and the imported_epitope_fragments.tsv is the realization of EpitopeFragment table in ViruSurf db.

b. extract epitopes, access them through a list of tuples
 * set up the paths for your machine in [integrate_iedb_epitopes.py](code/integrate_iedb_epitopes.py)
 * run [integrate_iedb_epitopes.py](code/integrate_iedb_epitopes.py) and get the list of tuples, where each tuple is one row of Epitope or EpitopeFragment table.
