from os.path import join, exists
from os import scandir
from code.utils import is_fasta_file_extension
from pandas import read_csv
from Bio import SeqIO


class HomologyBasedEpitopes:
	"""
	Class to extract epitopes for subject virus sequence
	based on homology searching of experimental identified immunodominant regions

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: BLAST file of immunodominant regions of relative virus (sars cov1) to target virus (sars cov2)
	output: file containing the aligned sequences, response frequency and identity
	"""

	def __init__(self, ncbi_ids, data_path, target_proteins_folder, blast_folder, out_path):
		"""
		HomologyBasedEpitopes constructor

		Parameters
		----------
		ncbi_ids : dict of str: str
			dictionary to keep the NCBI ids for target organism and its close-related
		data_path : str
			input data root
		target_proteins_folder : str
			target virus proteins folder
		blast_folder : str
			blast folder
		out_path : str
			output data root

		Returns
		-------
		None
		"""
		self.ncbi_ids = ncbi_ids
		self.data_path = data_path
		self.target_proteins_folder = join(data_path, target_proteins_folder)
		self.blast_folder = join(data_path, blast_folder)
		self.out_path = out_path

	@staticmethod
	def process_relative_protein_info(protein_info, process_tcells):
		"""
		Process close relative virus protein info tag to get
		protein id, start end of immunodominant region, max RF (response frequency)

		Parameters
		----------
		protein_info : str

		Returns
		-------
		protein_id : str
			relative virus protein id
		protein_start : str
			start of immunodominant region
		protein_end : str
			end of immunodominant region
		protein_RF : str
			immunodominant region RF score (max of sliding window values: B-cells, RF score calculated by assays: T-cells)
		process_tcells : bool
			process T-cells (True), otherwise process B-cells
		"""
		protein_info_parts = protein_info.split("|")
		protein_id = protein_info_parts[0].split("_")[0]
		protein_start_end = protein_info_parts[1].split("=")[1].split("-")
		protein_start, protein_end = protein_start_end[0], protein_start_end[1]
		protein_RF = protein_info_parts[2].split("=")[1]
		if process_tcells:
			hla_allele = protein_info_parts[3].split("=")[1]
		else:
			hla_allele = "not applicable"
		return protein_id, protein_start, protein_end, protein_RF, hla_allele

	def blast_out2csv(self, blast_out_file, process_tcells):
		"""
		Convert a blast alignment (output format = 7, commented tabular) to csv file

		Parameters
		----------
		blast_out_file : str
			blast output filename
		process_tcells : bool
			process T-cells (True), otherwise process B-cells (False)

		Returns
		-------
		None
		"""
		print("Convert blast output file to homology based epitopes csv")
		epitopes_name = "homology_based_epitopes.tsv"
		if not exists(join(self.out_path, epitopes_name)):
			print("Create file: {}".format(epitopes_name))
			with open(join(self.out_path, epitopes_name), "w") as epitopes_out:
				epitopes_out.write(
					"target_organism\ttarget_prot_id\ttarget_epi_start\ttarget_epi_end\ttarget_epi_sequence\trelative_organism\trelative_prot_id\trelative_epi_start\trelative_epi_end\trelative_epi_sequence\trelative_RF_score\tHLA restriction\tblast_identity\tprediction_method\n")

		with open(join(self.out_path, epitopes_name), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			blast_alignments_df = read_csv(join(self.blast_folder, blast_out_file), sep="\t", header=0)
			target_prot_is_loaded, relative_prot_is_loaded = False, False

			for index, blast_alignment in blast_alignments_df.iterrows():
				relative_prot_id, relative_epi_start, relative_epi_end, relative_RF_score, hla_allele = HomologyBasedEpitopes.process_relative_protein_info(
					blast_alignment["qseqid"], process_tcells)
				target_prot_id = blast_alignment["sseqid"].split("|")[3]
				target_start, target_end = str(blast_alignment["sstart"]), str(blast_alignment["send"])
				blast_identity = str(blast_alignment["pident"])

				# read target organism protein record
				if not target_prot_is_loaded:
					target_prot_rec = SeqIO.read(join(self.target_proteins_folder, target_prot_id + ".fasta"), "fasta")
					target_prot_is_loaded = True
				assert int(
					target_start) >= 1, "AssertionError: BLAST should report 1-index alignments, current alignment start = {}".format(
					target_start)
				target_epi_seq = str(target_prot_rec.seq[int(target_start) - 1:int(target_end)])
				assert target_epi_seq != "", "AssertionError: For relative protein id: {}, region: {}-{}, target epitope sequence should not be empty string".format(
					relative_prot_id, relative_epi_start, relative_epi_end)

				# read close relative organism protein record
				if not relative_prot_is_loaded:
					relative_epi_seqs = {}
					for rel_epi in SeqIO.parse(join(self.out_path, "immunodom_regions_" + relative_prot_id + ".fasta"),
					                           "fasta"):
						relative_epi_seqs[rel_epi.id] = str(rel_epi.seq)
					relative_prot_is_loaded = True
				assert blast_alignment[
					       "qseqid"] in relative_epi_seqs, "AssertionError: Close relative virus protein should contain epitope with id: {}".format(
					blast_alignment["qseqid"])
				relative_epi_seq = relative_epi_seqs[blast_alignment["qseqid"]]

				print("Write homology based epitope")
				epitope_row = "\t".join(
					[self.ncbi_ids["target_organism"], target_prot_id, target_start, target_end, target_epi_seq,
					 self.ncbi_ids["relative_organism"], relative_prot_id, relative_epi_start, relative_epi_end,
					 relative_epi_seq, relative_RF_score, hla_allele, blast_identity, "homology"])
				epitopes_out.write(epitope_row + "\n")

	def find_epitopes(self, process_tcells):
		"""
		Find homology based epitopes

		Parameters
		----------
		process_tcells : bool
			process T-cells (True), otherwise process B-cells (False)

		Returns
		-------
		None
		"""
		with scandir(self.blast_folder) as blast_dir:
			for blast_content in blast_dir:
				if blast_content.name[-3:] == "tsv" and blast_content.is_file():
					print("Process blast output file: {}".format(blast_content.name))
					self.blast_out2csv(blast_content.name, process_tcells)
