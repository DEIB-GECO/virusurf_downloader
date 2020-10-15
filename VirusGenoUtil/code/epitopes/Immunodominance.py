from VirusGenoUtil.code.epitopes.Protein import Protein
from VirusGenoUtil.code.utils import is_fasta_file_extension
from os.path import join
from os import scandir
from pandas import read_csv, Series, unique
import numpy as np
import VirusGenoUtil.matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from math import isnan, sqrt


class Immunodominance:
	"""
	Class to calculate the immunodominance regions per amino acid position
	for a virus with experimentally identified epitopes based on the IEDB immunome browser

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: immunobrowser csv file
	output: sliding average of lower bound of 95% confidence interval (CI) of response frequency (RF)
	"""

	def __init__(self, data_path, proteins_folder, exp_epitopes_folder, assay_species, out_path):
		"""
		Immunodominance constructor

		Parameters
		----------
		data_path : str
			input data root
		proteins_folder : str
			proteins folder
		exp_epitopes_folder : str
			experimental epitopes folder
		assay_species : str
			name of species on which epitopes response frequency was measured
		out_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.proteins_path = join(data_path, proteins_folder)
		self.exp_epitopes_path = join(data_path, exp_epitopes_folder)
		self.assay_species = assay_species
		self.out_path = out_path
		self.proteins = {}

	def load_proteins(self):
		"""
		Load all proteins of virus with experimental epitopes

		Parameters
		----------

		Returns
		-------
		None
		"""
		with scandir(self.proteins_path) as proteins_dir:
			for proteins_content in proteins_dir:
				print(proteins_content.name)
				if is_fasta_file_extension(proteins_content.name) and proteins_content.is_file():
					protein = Protein(join(self.proteins_path, proteins_content.name))
					protein.load()
					protein_id = proteins_content.name.split(".")[0]
					if protein_id not in self.proteins:
						print("Saving protein with id: {}".format(protein_id))
						protein_rec = protein.get_record()
						self.proteins[protein_id] = protein_rec
		print(self.proteins)

	def load_immunome_csv(self, protein_id):
		"""
		Load immunome csv

		Parameters
		----------

		Returns
		-------
		Pandas.DataFrame
			loaded immunome dataframe
		"""
		print("Load immunome csv")
		immunome_name = "RF_" + self.assay_species + "_" + protein_id + ".csv"
		immunome_df = read_csv(join(self.exp_epitopes_path, immunome_name), sep=",", header=0)
		return immunome_df

	def sum_per_position(self, immunome, protein_id):
		"""
		Sum RF per amino acid position using the start-end positions and lower bound values

		Parameters
		----------
		immunome : Pandas.DataFrame
			immunome pandas dataframe
		protein_id : str
			id of protein currently processed

		Returns
		-------
		None
		"""
		epitope_position = int(immunome["position"]) - 1  # 0-index
		epitope_lower_bound = float(immunome["lowerbound"])
		self.proteins[protein_id].letter_annotations["lower_bound_rf"][epitope_position] += epitope_lower_bound

	def compute_sliding_avg(self, protein_id, window_size, save_plot):
		"""
		Average over sliding window to define contiguous regions for lower bound CI of RF
		use window/2-1 from the left and the right of the current base RF
		mean of current base = |--(window/2)-1----|current base|----(window/2)-1--|
		Reference: Supplementary of Grifoni et al.

		Credits: https://stackoverflow.com/questions/20249149/rolling-mean-with-customized-window-with-pandas#20249295
		Parameters
		----------
		protein_id : str
			id of protein to compute average sliding window
		window_size : int
			sliding window size
		save_plot : bool
			create and save plot of sliding window values

		Returns
		-------
		list of float
			sliding_average
		"""
		print("Compute sliding average RF")
		lower_bound_rf = np.asarray(self.proteins[protein_id].letter_annotations["lower_bound_rf"])
		sliding_average = Series(lower_bound_rf).rolling(window=int((window_size / 2) - 1), min_periods=1,
		                                                 center=True).mean()
		assert len(lower_bound_rf) == len(sliding_average), "Sliding average values are less the starting RF values"

		if save_plot:
			fig = plt.figure()
			plt.plot(sliding_average, color='k')
			plt.axhline(y=0.3, color='k', linestyle='--')
			plt.ylim(0, 1)
			plt.xlim(0, len(self.proteins[protein_id].seq))
			plt.xlabel("Amino acid position")
			plt.ylabel("Lower bound\nsliding average")
			plt.title(self.proteins[protein_id].id)
			avg_plot_name = "sliding_avg_" + protein_id + ".png"
			fig.savefig(join(self.out_path, avg_plot_name), bbox_inches='tight', dpi=600)
			print("Saved {}".format(join(self.out_path, avg_plot_name)))
		return sliding_average.tolist()

	def set_RF_zero(self, protein_rec):
		"""
		Set protein sequence feature RF to 0
		across all its amino acid

		Parameters
		----------
		protein_rec :

		Returns
		-------
		Bio.SeqIO.SeqRecord
			protein with set lower
		"""
		print("Init lower bound RF to 0")
		protein_rec.letter_annotations["lower_bound_rf"] = len(protein_rec.seq) * [0.0]
		protein_rec.letter_annotations["sliding_avg_lower_bound"] = len(protein_rec.seq) * [0.0]
		return protein_rec

	def process_all_proteins(self, process_tcells, window_size, save_plot, sliding_avg_cutoff, split_at_min_RF,
	                         save_immunodominant_reg):
		"""
		Process
		* Immunodominance regions for each loaded protein (B-cells)
		* IEBD assay csv (T-cells)
		in order to identify immunodominant regions and save them in fasta file

		Parameters
		----------
		process_tcells : bool
			T cells are to be processed (True) otherwise B cells are to be processed (False)
		window_size : int
			sliding window size
		save_plot : bool
			create and save plot of sliding window values (True) otherwise don't save (False)
		sliding_avg_cutoff : float
			cutoff of sliding average of lower bound RF to find immunodominant regions
		split_at_min_RF : bool
			split a region into sub-regions if minimum RF is too lower than max RF (True) otherwise don't split (False)
		save_immunodominant_reg : bool
			save sequence of each immunodominant region (True) otherwise don't save (False)

		Returns
		-------
		dict of str: Bio.SeqIO.SeqRecord
			processed proteins
		"""
		print("RF=0 -> load(RF,immunome.csv) -> slide window")
		self.load_proteins()
		for protein_id, protein_rec in self.proteins.items():
			print("--- ---")
			print("Process protein: {}".format(protein_id))
			if process_tcells:
				# for T-cells, retrieve epitopes from the IEDB csv files
				self.retrieve_Tcells_epitopes(protein_id, save_immunodominant_reg)
			else:
				# for B-cells, parse the immunobrowser lower bound RF value
				# compute the sliding window
				# and identify immunodominant regions using a cut-off value
				self.proteins[protein_id] = self.set_RF_zero(protein_rec)
				# read immuno csv
				immunome_df = self.load_immunome_csv(protein_id)
				# sum lower bound RF per position
				print("Sum lower bound RF per epitope and amino-acid position")
				immunome_df.apply(self.sum_per_position, protein_id=protein_id, axis=1)
				self.proteins[protein_id].letter_annotations["sliding_avg_lower_bound"] = self.compute_sliding_avg(
					protein_id,
					window_size,
					save_plot)

				self.find_immunodominant_regions(protein_id, sliding_avg_cutoff, split_at_min_RF,
				                                 save_immunodominant_reg)
		return self.proteins

	def map_epitope2allele(self, idbe_assay):
		"""
		Map epitopes to allele information

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		save_epitopes : bool
			Save retrieved epitopes sequences (True), don't save (False)

		Returns
		-------
		dict of str : str
			dictionary mapping unique epitopes to their allele info
		"""
		print("Map unique epitope sequence to allele information")
		uniq_epitopes2allele = {}
		unique_epitopes = unique(idbe_assay["Description"])

		for uniq_epi in unique_epitopes:
			if uniq_epi not in uniq_epitopes2allele:
				all_allele_names = []
				for _, row in idbe_assay.loc[idbe_assay["Description"] == uniq_epi, ["Allele Name"]].iterrows():
					allele = str(row["Allele Name"])
					if "HLA" not in allele:
						allele = "unknown"
					else:
						allele = allele.replace(" ", "_")
					if allele not in all_allele_names:
						all_allele_names.append(allele)

				if "unknown" in all_allele_names and len(all_allele_names) > 1:
					all_allele_names.remove("unknown")
				uniq_epitopes2allele[uniq_epi] = ",".join(all_allele_names)

		return uniq_epitopes2allele

	def calculate_RF_score(self, idbe_assay, uniq_epitopes):
		"""
		Calculate RF score for T cell assays, using:
		RF = (r-sqrt(r))/t,
		where r is the # positive responding assays
		and t is the # of the total tested assays

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		uniq_epitopes : list of str
			list of unique epitope for which the RF score will be calculated

		Returns
		-------
		dict of str : float
			dictionary with key the unique epitope and the value the RF score
		"""
		print("Map unique epitope sequence to RF score")
		rf_scores = {}
		for epi in uniq_epitopes:
			# sum r and t for all assays testing the current epitopte

			t, r = 0.0, 0.0
			for _, row in idbe_assay.loc[idbe_assay["Description"] == epi, ["Number of Subjects Tested",
			                                                                "Number of Subjects Responded"]].iterrows():
				if isnan(row["Number of Subjects Tested"]):
					t = t + 1.0
				else:
					t = t + row["Number of Subjects Tested"]
				if isnan(row["Number of Subjects Responded"]):
					r = r + 1.0
				else:
					r = r + row["Number of Subjects Responded"]
			if t == 0:
				rf_score = 0.0
			else:
				rf_score = (r - sqrt(r)) / t
			rf_scores[epi] = rf_score
		return rf_scores

	def find_epitope_regions(self, protein_rec, epitopes):
		"""
		Find the epitopes regions (start-stop) in the protein record
		Credits: Chapter 20.1.8 Biopython (http://biopython.org/DIST/docs/tutorial/Tutorial.html)
		Parameters
		----------
		protein_rec : Bio.SeqIO.SeqRecord
			protein record
		epitopes : list of str
			epitope sequences

		Returns
		-------
			list of list of int
			sorted list of start end position of epitopes regions
		"""
		print("Map unique epitope sequence to start end positions")
		epi_regions = []
		for epi in epitopes:
			start = protein_rec.seq.find(epi)
			if start != -1:
				end = start + len(epi) - 1
				epi_regions.append([start, end])
		epi_regions.sort(key=lambda x: x[0])  # sort ascending by starting position
		return epi_regions

	def retrieve_Tcells_epitopes(self, protein_id, save_epitopes):
		"""
		Retrieve epitopes and their allele info from the saved csv generated by the assay tab of IEDB

		Parameters
		----------
		protein_id : str
			if of protein
		save_epitopes : bool
			save retrieved epitopes in fasta (True) otherwise don't save

		Returns
		-------
		None
		"""
		print("Retrieve and save T-cells epitopes")

		# load assay csv
		tcell_assay = read_csv(join(self.exp_epitopes_path, "tcell_" + protein_id + ".csv"), sep=",", header=1)
		epitope2allele = self.map_epitope2allele(tcell_assay)

		# save unique epitopes into fasta file
		protein_record = self.proteins[protein_id]

		# find epitope regions
		epi_regions = self.find_epitope_regions(protein_record, list(epitope2allele.keys()))

		# calculate RF score
		epi2rf_score = self.calculate_RF_score(tcell_assay, list(epitope2allele.keys()))

		# for each epitope region with RF score > 0
		# collect all info per region and save region
		epi_frags = []
		for i, region in enumerate(epi_regions):
			reg_start, reg_end = region[0], region[-1] + 1
			epi_frag = protein_record.seq[reg_start:reg_end]
			assert epi_frag in epi2rf_score, "Epitope with sequence {} does not have calculated RF score".format(
				epi_frag)
			rf_score = float("{:.4f}".format(epi2rf_score[epi_frag]))
			if rf_score > 0.0:
				frag_record = SeqRecord(epi_frag, id=protein_id + "_immunodom_frag_" + str(i + 1) + "|reg=" + str(
					reg_start + 1) + "-" + str(reg_end) + "|RF_score=" + str(rf_score)
				                                     + "|HLA=" + epitope2allele[protein_record.seq[reg_start:reg_end]],
				                        name="",
				                        description="")
				epi_frags.append(frag_record)

		if save_epitopes:
			SeqIO.write(epi_frags, join(self.out_path, "immunodom_regions_" + protein_id + ".fasta"), "fasta")
			print("Saved regions at: {}".format(join(self.out_path, "immunodom_regions_" + protein_id + ".fasta")))

	def find_immunodominant_regions(self, protein_id, sliding_avg_cutoff, split_at_min_RF, save_regions):
		"""
		Find immunodominant regions in input protein

		Parameters
		----------
		protein_id : str
			id of protein
		sliding_avg_cutoff : float
			cutoff of sliding average of lower bound RF to find immunodominant regions
		split_at_min_RF : bool
			split a region into sub-regions if minimum RF exists
		save_regions : bool
			save identified regions in fasta file (True), otherwise don't save (False)

		Returns
		-------
		None
		"""
		print("Identify and save immunodominant regions")
		protein_record = self.proteins[protein_id]
		# get all positions with average >= cutoff
		immunodominant_positions = [aa_pos for aa_pos, lower_bound_avg in
		                            enumerate(protein_record.letter_annotations["sliding_avg_lower_bound"]) if
		                            lower_bound_avg >= sliding_avg_cutoff]

		# find consecutive position to define a region
		immunodominant_regions, current_region = [], []
		previous_pos = immunodominant_positions[0]
		for aa_pos in immunodominant_positions:
			if aa_pos - previous_pos > 1:  # create a new subregion
				immunodominant_regions.append(current_region)
				current_region = [aa_pos]
			else:
				current_region.append(aa_pos)
			previous_pos = aa_pos
		immunodominant_regions.append(current_region)

		# split immunodominant regions if a very minimum RF exists
		if split_at_min_RF:
			splitted_immunodominant_regions = []
			for region in immunodominant_regions:
				reg_start, reg_end = region[0], region[-1] + 1
				max_RF = max(protein_record.letter_annotations["sliding_avg_lower_bound"][reg_start:reg_end])
				min_RF = min(protein_record.letter_annotations["sliding_avg_lower_bound"][reg_start:reg_end])
				min_RF_pos = reg_start + protein_record.letter_annotations["sliding_avg_lower_bound"][
				                         reg_start:reg_end].index(min_RF)
				if split_at_min_RF and max_RF - min_RF >= sliding_avg_cutoff and min_RF_pos + 1 < reg_end:
					# split if the min is much lower than max, but also min is NOT in the end of the region
					print("Splitting immunodominant regions to shorter as minimum value too lower than maximum value")
					min_RF_pos = reg_start + protein_record.letter_annotations["sliding_avg_lower_bound"][
					                         reg_start:reg_end].index(min_RF)
					sub_region1, sub_region2 = [reg_start, min_RF_pos], [min_RF_pos + 1, reg_end]
					splitted_immunodominant_regions.append(sub_region1)
					splitted_immunodominant_regions.append(sub_region2)
				else:
					splitted_immunodominant_regions.append(region)
			immunodominant_regions = splitted_immunodominant_regions

		# save protein fragment for each immunodominant region
		immunodom_frags = []
		for i, region in enumerate(immunodominant_regions):
			reg_start, reg_end = region[0], region[-1] + 1
			if reg_end - reg_start + 1 >= 10:
				immunodom_frag = protein_record.seq[reg_start:reg_end]
				max_lower_bound = float("{:.4f}".format(
					max(protein_record.letter_annotations["sliding_avg_lower_bound"][reg_start:reg_end])))
				frag_record = SeqRecord(immunodom_frag, id=protein_id + "_immunodom_frag_" + str(i + 1) + "|reg=" + str(
					reg_start + 1) + "-" + str(reg_end) + "|max_lower_bound=" + str(max_lower_bound), name="",
				                        description="")
				immunodom_frags.append(frag_record)
			else:
				print("Skip region with less than 10 amino-acids")
				print("Skipped region start={} & end={}".format(reg_start, reg_end))
		if save_regions:
			SeqIO.write(immunodom_frags, join(self.out_path, "immunodom_regions_" + protein_id + ".fasta"), "fasta")
			print("Saved regions at: {}".format(join(self.out_path, "immunodom_regions_" + protein_id + ".fasta")))
