from code.utils import is_fasta_file_extension
from code.epitopes.Protein import Protein
from code.epitopes.Epitope import Epitope
from code.epitopes.EpitopeFragment import EpitopeFragment
from os.path import join, splitext, isfile, exists, isdir
from pathlib import Path
from os import scandir
from pandas import read_csv, unique
from math import isnan, sqrt


class IEDBEpitopes:
	"""
	Class to extract B and T-cells epitopes for
	each virus taxon id found in virus_protein folder

	Input epitopes per cell type are downloaded from http://www.iedb.org/database_export_v3.php
	IEDB paper reference:
	Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, Wheeler DK, Sette A, Peters B.
	The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res. 2019 Jan 8;47(D1):D339-D343.
	doi: 10.1093/nar/gky1006. PMID: 30357391
	"""

	def __init__(self, data_path, cell_epitopes_folder, viruses_folder, host_taxon_id, host_name, assay_type,
	             output_path):
		"""
		IEDBEpitopes constructor

		Parameters
		----------
		data_path : str
			input data root
		cell_epitopes_folder : str
			B and T cell IEDB epitopes folder
		viruses_folder : str
			viruses protein folder
		host_taxon_id : str
			NCBI taxon id of host organism whose cells used for the experimental identification of epitopes
		host_name : str
			name of host organism whose cells used for the experimental identification of epitopes
		assay_type : str
			IEDB experimental assay type
		output_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.cell_epitopes_path = join(data_path, cell_epitopes_folder)
		self.viruses_path = join(data_path, viruses_folder)

		self.tcell_iedb_assays, self.bcell_iedb_assays = None, None
		self.host_taxon_id = host_taxon_id.split("_")[1]
		self.host_name = str(host_name)
		self.assay_type = assay_type
		self.url_prefixes = {"taxid": "http://purl.obolibrary.org/obo/NCBITaxon_",
		                     "uniprot": "http://www.uniprot.org/uniprot/"}
		self.output_path = output_path
		self.current_virus_proteins = {}
		self.current_virus_taxon_id = None
		self.current_virus_epitopes = []  # all Epitope objects for current virus
		self.current_virus_epi_fragments = []  # all EpitopeFragment objects for current virus
		self.ncbi_iedb_not_equal = []  # save all epitope sequence information that IEDB seq is not equal to NCBI
		self.epitope_id = 0
		self.epitope_fragment_id = 0
		self.load_iedb_csvs()

	def get_virus_epi_fragments(self):
		"""
		Get current virus epitope fragments in list of tuples

		Returns
		-------
		list of tuple
			all epitope fragments of virus
		"""
		print("Return all epitope fragments of virus taxid={}".format(self.current_virus_taxon_id))
		epi_fragments = []
		for epi_fragment in self.current_virus_epi_fragments:
			epi_frag_attributes = epi_fragment.get_all_attributes()
			fragment = tuple([int(epi_frag_attributes["fragment_id"]), int(epi_frag_attributes["parent_epi_id"]),
			                  epi_frag_attributes["fragment_seq"],
			                  int(epi_frag_attributes["fragment_start"]),
			                  int(epi_frag_attributes["fragment_stop"])])
			epi_fragments.append(fragment)
		return epi_fragments

	def virus_epi_fragments2tsv(self):
		"""
		Write fragments of discontinuous epitopes for current virus to tsv file

		Returns
		-------
		None
		"""
		print("Save current virus epitope fragments to csv")
		epitopes_name = "imported_iedb_epitope_fragment.tsv"
		if not exists(join(self.output_path, epitopes_name)):
			print("Create file: {}".format(epitopes_name))
			with open(join(self.output_path, epitopes_name), "w") as epitopes_out:
				epitopes_out.write(
					"epi_fragment_id\tparent_epitope_id\tepi_fragment_sequence\tepi_fragment_start\tepi_fragment_stop\n")
		with open(join(self.output_path, epitopes_name), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			for epi_fragment in self.current_virus_epi_fragments:
				print("Write IEDB imported epitope fragment")
				epi_frag_attributes = epi_fragment.get_all_attributes()
				epi_frag_row = "\t".join(
					[str(epi_frag_attributes["fragment_id"]), str(epi_frag_attributes["parent_epi_id"]),
					 epi_frag_attributes["fragment_seq"],
					 str(epi_frag_attributes["fragment_start"]),
					 str(epi_frag_attributes["fragment_stop"])])
				epitopes_out.write(epi_frag_row + "\n")
		print("====")

	def get_virus_epitopes(self):
		"""
		Get current virus epitopes in list of tuples

		Returns
		-------
		list of tuple
			all epitopes of virus
		"""
		print("Return all epitopes of virus taxid={}".format(self.current_virus_taxon_id))
		epitopes = []
		for epi in self.current_virus_epitopes:
			epi_attributes = epi.get_all_attributes()
			if epi_attributes["is_linear"]:
				epi_seq = epi_attributes["region_seq"]
			else:
				epi_seq = None
			epitope = tuple([int(epi_attributes["epitope_id"]), int(epi_attributes["virus_taxid"]),
			                 str(epi_attributes["host_taxid"]),
			                 str(epi_attributes["protein_ncbi_id"]), str(epi_attributes["cell_type"]),
			                 epi_attributes["hla_restriction"],
			                 float(epi_attributes["response_frequency"]), epi_seq,
			                 int(epi_attributes["region_start"]), int(epi_attributes["region_stop"]),
			                 epi_attributes["is_imported"], ",".join(epi_attributes["external_links"]),
			                 epi_attributes["prediction_process"], epi_attributes["is_linear"]])
			epitopes.append(epitope)
		return epitopes

	def virus_epitopes2tsv(self):
		"""
		Write epitopes for current virus to tsv file

		Returns
		-------
		None
		"""
		print("Save current virus epitopes to csv")

		epitopes_name = "imported_iedb_epitopes.tsv"
		if not exists(join(self.output_path, epitopes_name)):
			print("Create file: {}".format(epitopes_name))
			with open(join(self.output_path, epitopes_name), "w") as epitopes_out:
				epitopes_out.write(
					"epitope_id\tvirus_taxid\thost_taxid\tprotein_ncbi_id\tcell_type\thla_restriction\tresponse_frequency\tepitope_sequence\tepitope_start\tepitope_stop\tis_imported\texternal_links\tprediction_process\tis_linear\n")

		with open(join(self.output_path, epitopes_name), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			for epitope in self.current_virus_epitopes:
				print("Write IEDB imported epitope")
				epi_attributes = epitope.get_all_attributes()

				if epi_attributes["is_linear"]:
					epi_seq = epi_attributes["region_seq"]
				else:
					epi_seq = None
				epitope_row = "\t".join(
					[str(epi_attributes["epitope_id"]), epi_attributes["virus_taxid"], epi_attributes["host_taxid"],
					 epi_attributes["protein_ncbi_id"], epi_attributes["cell_type"], epi_attributes["hla_restriction"],
					 str(epi_attributes["response_frequency"]), epi_seq,
					 str(epi_attributes["region_start"]), str(epi_attributes["region_stop"]),
					 str(epi_attributes["is_imported"]), ",".join(epi_attributes["external_links"]),
					 str(epi_attributes["prediction_process"]), str(epi_attributes["is_linear"])])
				epitopes_out.write(epitope_row + "\n")
		print("====")

	def load_iedb_csvs(self):
		"""
		Load IEDB B and T cells assays csvs

		Returns
		-------
		None
		"""
		assert isfile(join(self.cell_epitopes_path,
		                   "tcell_full_v3.csv.gz")), "AssertionError: IEDB Tcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		self.tcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "tcell_full_v3.csv.gz"), sep=",", header=1,
		                                  compression='gzip')
		assert isfile(join(self.cell_epitopes_path,
		                   "bcell_full_v3.csv.gz")), "AssertionError: IEDB Bcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		self.bcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "bcell_full_v3.csv.gz"), sep=",", header=1,
		                                  compression='gzip')

	def subset_iedb_by_host_assay_type(self):
		"""
		Subset IEDB records by host taxon id and assay type

		Returns
		-------
		None
		"""
		print("Get working subset of iedb assays")
		if self.assay_type == "positive":
			print("Select positive assays")
			self.tcell_iedb_assays = self.tcell_iedb_assays.loc[
				self.tcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]

			self.bcell_iedb_assays = self.bcell_iedb_assays.loc[
				self.bcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]

		print("Select epitopes experimental identified in host id: {} = {}".format(self.host_taxon_id, self.host_name))
		host_taxid = self.url_prefixes["taxid"] + self.host_taxon_id
		self.tcell_iedb_assays = self.tcell_iedb_assays.loc[self.tcell_iedb_assays["Host IRI"] == host_taxid]
		self.bcell_iedb_assays = self.bcell_iedb_assays.loc[self.bcell_iedb_assays["Host IRI"] == host_taxid]
		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))
		print("T cell subset number of non-unique epitopes: {}".format(self.tcell_iedb_assays.shape[0]))

	def subset_Bcells_by_epi_type(self):
		"""
		Subset B cells to keep only B-cells linear and discontinuous epitopes

		Returns
		-------
		None
		"""
		print("Select only linear or discontinuous epitopes")
		self.bcell_iedb_assays = self.bcell_iedb_assays.loc[
			(self.bcell_iedb_assays["Object Type"] == "Discontinuous peptide") | (
					self.bcell_iedb_assays["Object Type"] == "Linear peptide")]
		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))

	def subset_iedb_by_virus_id(self):
		"""
		Subset iedb record to get the ones related to virus id

		Returns
		-------
		Pandas.DataFrame, Pandas.DataFrame
			B-cell IEDB assays related only to virus id, T-cell IEDB assay only to virus id
		"""
		print("Get iedb only for taxon id={}".format(self.current_virus_taxon_id))

		tcell_iedb_virus = self.tcell_iedb_assays.loc[
			self.tcell_iedb_assays["Parent Species IRI"].str.split("NCBITaxon_").str[
				-1].str.strip() == self.current_virus_taxon_id]
		if tcell_iedb_virus.shape[0] == 0:
			print("2nd attempt: Match with Organism IRI")
			tcell_iedb_virus = self.tcell_iedb_assays.loc[
				self.tcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[
					-1].str.strip() == self.current_virus_taxon_id]
		bcell_iedb_virus = self.bcell_iedb_assays.loc[
			self.bcell_iedb_assays["Parent Species IRI"].str.split("NCBITaxon_").str[
				-1].str.strip() == self.current_virus_taxon_id]
		if bcell_iedb_virus.shape[0] == 0:
			print("2nd attempt: Match with Organism IRI")
			bcell_iedb_virus = self.bcell_iedb_assays.loc[
				self.bcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[
					-1].str.strip() == self.current_virus_taxon_id]
		print("B cell subset for virus contains {} non unique epitopes".format(bcell_iedb_virus.shape[0]))
		print("B cell head: {}".format(bcell_iedb_virus.head()))
		print("T cell subset for virus contains {} non unique epitopes".format(tcell_iedb_virus.shape[0]))
		print("T cell head: {}".format(tcell_iedb_virus.head()))
		return bcell_iedb_virus, tcell_iedb_virus

	def process_all_viruses(self):
		"""
		Process and get IEDB epitopes for each protein of each virus

		Returns
		-------
		None
		"""
		print("Process all viruses found in {}".format(self.viruses_path))
		self.subset_iedb_by_host_assay_type()
		self.subset_Bcells_by_epi_type()
		with scandir(self.viruses_path) as viruses_dir:
			for content in viruses_dir:
				if content.is_dir() and "taxon_" in content.name:
					self.process_virus_proteins(content)
					self.virus_epitopes2tsv()
					self.virus_epi_fragments2tsv()
		self.write_ncbi_iedb_discordant_sequences()
		print("=== ~ ===")

	def process_virus(self, virus_taxid):
		"""
		Extract virus epitopes from IEDB

		Parameters
		----------
		virus_taxid : int
		NCBI taxonomy id of virus

		Returns
		-------
		None
		"""
		print("Process virus with ncbi taxid={}".format(virus_taxid))
		virus_path = Path(join(self.viruses_path, "taxon_" + str(virus_taxid)))
		assert isdir(virus_path), "AssertionError: virus folder was not found in viruses_proteins folder"
		self.subset_iedb_by_host_assay_type()
		self.subset_Bcells_by_epi_type()
		self.process_virus_proteins(virus_path)

	def write_ncbi_iedb_discordant_sequences(self):
		"""
		Write all epitopes that their NCBI and IEDB sequences do not match

		Returns
		-------
		None
		"""
		print("Write epitope information for sequences discordant between NCBI and IEDB at {}".format(
			join(self.output_path, "ncbi_iedb_discordant_epitopes.txt")))
		with open(join(self.output_path, "ncbi_iedb_discordant_epitopes.txt"), 'w') as discordant_out:
			discordant_out.writelines("\n".join(self.ncbi_iedb_not_equal))

	def load_virus_proteins(self, virus_proteins_path):
		"""
		Load virus proteins from argument path

		Parameters
		----------
		virus_proteins_path : str
			virus proteins path

		Returns
		-------
		None
		"""
		print("Load virus proteins")
		self.current_virus_proteins = {}
		with scandir(virus_proteins_path) as proteins_dir:
			for proteins_content in proteins_dir:
				if is_fasta_file_extension(proteins_content.name) and proteins_content.is_file():
					protein = Protein(join(virus_proteins_path, proteins_content.name))
					protein.load()
					protein_id = splitext(proteins_content.name)[0]
					if protein_id not in self.current_virus_proteins:
						self.current_virus_proteins[protein_id] = protein
		print("====")

	def process_virus_proteins(self, virus_proteins_path):
		"""
		Process information for IEDB epitopes for each protein of the virus

		Parameters
		----------
		virus_proteins_path : str
			virus protein folder path

		Returns
		-------
		None
		"""
		self.current_virus_taxon_id = splitext(virus_proteins_path.name)[0].split("_")[1]
		print("=== Virus ===")
		print("Process proteins of virus with taxon id: {}".format(self.current_virus_taxon_id))
		self.load_virus_proteins(virus_proteins_path)
		bcells_current_virus, tcells_current_virus = self.subset_iedb_by_virus_id()
		self.current_virus_epitopes = []  # clear current virus epitopes
		self.current_virus_epi_fragments = []  # clear current virus epitope fragments
		self.ncbi_iedb_not_equal.append("=== Virus taxid={} ===".format(self.current_virus_taxon_id))

		self.process_Bcells(bcells_current_virus)
		self.process_Tcells(tcells_current_virus)

	def process_Bcells(self, bcells_current_virus):
		"""
		Process B-cells assays for current virus

		Parameters
		----------
		bcells_current_virus : Pandas.DataFrame
			B-cells assays for current virus

		Returns
		-------
		None
		"""
		print("\nProcess B-cells")
		self.ncbi_iedb_not_equal.append("B-cells:")
		for protein_id, protein in self.current_virus_proteins.items():
			print("---")
			print("Process protein with uniprot id: {}".format(protein_id))
			iedb_uniprot_id = self.url_prefixes["uniprot"] + protein.get_uniprot_id()
			bcells_current_protein = bcells_current_virus.loc[
				bcells_current_virus["Parent Protein IRI"] == iedb_uniprot_id]
			if bcells_current_protein.shape[0] == 0:
				print("Could not match using uniprot")
				print("2nd attempt: Match with protein name")
				bcells_current_protein = bcells_current_virus.loc[
					bcells_current_virus["Parent Protein"].str.split("[").str[
						0].str.strip().str.lower() == protein.get_name().lower()]
			print("Number of non unique epitopes for protein = {}".format(bcells_current_protein.shape[0]))

			if bcells_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# normalize epitopes
				normalized2unique = self.normalize_unique_epitope_sequences(bcells_current_protein)
				# map epitope to allele
				normalized2allele = self.map_epitope2allele(bcells_current_protein, normalized2unique, False)
				# find epitope regions
				epi_regions, normalized2regions = self.find_epitope_regions(bcells_current_protein, normalized2unique)
				# calculate RF score
				normalized2rf_score = self.calculate_RF_score(bcells_current_protein, normalized2unique)
				# extract external links per unique epitope
				normalized2external_links = self.find_epitope_external_links(bcells_current_protein, normalized2unique)

				# create an Epitope object for each identified IEDB epitope
				host_taxon_id = self.host_taxon_id
				is_imported = True
				prediction_process = "IEDB_import"
				# get protein record
				protein_record = protein.get_record()
				for normalized, region in normalized2regions.items():
					if region[0] == -1 and region[-1] == -1:  # discontinuous epitope
						is_linear = False
						reg_start, reg_end = self.get_discontinous_epi_start_stop(normalized)
					else:  # linear epitope
						is_linear = True
						# decrease by one the region start to convert 1-index to 0-index
						reg_start, reg_end = region[0] - 1, region[-1]
						ncbi_prot_epi = protein_record.seq[reg_start:reg_end]

						if ncbi_prot_epi != normalized:  # append discordant epitopes
							self.ncbi_iedb_not_equal.append(
								"iedb frag:{}, ncbi prot:{}, iedb epitope link(s): {}".format(normalized, ncbi_prot_epi,
								                                                              " ".join(
									                                                              normalized2external_links[
										                                                              normalized])))
					external_links = normalized2external_links[normalized]
					if normalized2rf_score[
						normalized] > -1:  # if tested subjects information is given, export epitope
						epitope = Epitope(self.current_virus_taxon_id, protein.get_ncbi_id(), host_taxon_id, "B cell",
						                  normalized2allele[normalized], normalized2rf_score[normalized],
						                  str(normalized), reg_start,
						                  reg_end, is_imported, external_links, prediction_process, is_linear)
					self.current_virus_epitopes.append(epitope)
					if not is_linear:
						self.current_virus_epi_fragments = self.current_virus_epi_fragments + epitope.get_fragments()
				# print inspection: for each protein print last epitope
				all_attributes = self.current_virus_epitopes[-1].get_all_attributes()
				print(all_attributes)
		print("====")

	def get_discontinous_epi_start_stop(self, epi_description):
		"""
		Get discontinuous epitope start and stop position

		Parameters
		----------
		epi_description : str
			discontinuous epitope description

		Returns
		-------
		int, int
			discontinuous epitope start, discontinuous epitope end
		"""
		discontinuous_epi_aa_pos = epi_description.split(",")
		if epi_description[-1] != ",":
			aa_pos_start, aa_pos_end = discontinuous_epi_aa_pos[0].strip(), discontinuous_epi_aa_pos[-1].strip()
			pos_start = int(aa_pos_start[1:len(aa_pos_start)])
			pos_end = int(aa_pos_end[1:len(aa_pos_end)])
			assert pos_end - pos_start >= 1, "AssertionError: current discontinuous epitope: {} has a fragment with larger pos_start than pos_end".format(
				epi_description)
		else:  # only one position discontinuous epitope contain the comma at the end
			aa_pos = discontinuous_epi_aa_pos[0].strip()
			pos_start = int(aa_pos[1:len(aa_pos)])
			pos_end = int(aa_pos[1:len(aa_pos)])
		return pos_start, pos_end

	def normalize_unique_epitope_sequences(self, iedb_assay):
		"""
		Normalize unique epitope sequences
		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab

		Returns
		-------
		list of str
			normalized unique epitopes
		"""
		uniq_epitopes = unique(iedb_assay["Description"])
		normalized2unique = {}
		for uniq_epi in uniq_epitopes:
			if "," in uniq_epi:
				# case: Discontinuous peptide	T359, T363, K365, K390, G391, D392, R395,
				normalized = uniq_epi
			elif "+" in uniq_epi:
				# case: Linear peptide	ALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFEL + ACET(A1)
				# remove the + ACET() information from the sequence
				normalized = uniq_epi.split("+")[0].strip()
			elif uniq_epi[1:].isdigit():
				# case: Discontinuous peptide (of one base)	P462
				normalized = uniq_epi + ","
			else:
				# case: Linear peptide	APRITFGGPTDSTDNNQN
				assert not uniq_epi.isdigit(), "AssertionError: epitope={} should have been only aminoacid sequence".format(
					uniq_epi)
				normalized = uniq_epi
			normalized2unique[normalized] = uniq_epi
		return normalized2unique

	def map_epitope2allele(self, iedb_assay, normalized2unique, is_tcell):
		"""
		Map epitopes to allele information

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe created from csv of IEDB assay tab
		normalized2unique : dict of str : str
			dictionary mapping normalized epitope to unique epitope
		is_tcell : bool
			Processing T-cells (True), otherwise it is B-cells (False)

		Returns
		-------
		dict of str : str
			dictionary mapping normalzed epitopes to their allele info
		"""
		print("Map unique epitope sequence to allele information")
		normalized2allele = {}
		# for each pair normalized,unique find the allele information based on the unique allele
		# then save it for the corresponding normalized epitope
		for normalized_epi, uniq_epi in normalized2unique.items():
			if is_tcell:
				if uniq_epi not in normalized2allele:
					all_allele_names = []
					for _, row in iedb_assay.loc[iedb_assay["Description"] == uniq_epi, ["Allele Name"]].iterrows():
						allele = str(row["Allele Name"])
						if "HLA" not in allele:
							allele = "unknown"
						else:
							allele = allele.replace(" ", "_")
						if allele not in all_allele_names:
							all_allele_names.append(allele)

					if "unknown" in all_allele_names and len(all_allele_names) > 1:
						all_allele_names.remove("unknown")
					normalized2allele[normalized_epi] = ",".join(all_allele_names)
			else:
				normalized2allele[normalized_epi] = None
		return normalized2allele

	def find_epitope_regions(self, iedb_assay, normalized2unique):
		"""
		Find the epitopes regions using IEDB start end coordinates
		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe created from csv of IEDB assay tab
		normalized2unique : dict of str : str
			normalized epitopes to unique epitopes

		Returns
		-------
			list of list of int
			sorted list of start end position of epitopes regions
			dict of str : list of int
			dictionary to map normalized epitope sequence to start end position
		"""
		print("Map unique epitope sequence to start end positions")
		epi_regions = []
		epi2regions = {}
		for normalized, unique in normalized2unique.items():
			if "," in normalized:  # discontinuous epitope
				start, end = -1, -1
				epi_regions.append([start, end])
			else:
				for _, row in iedb_assay.loc[iedb_assay["Description"] == unique].iterrows():
					if not isnan(row["Starting Position"]):
						start = int(row["Starting Position"])
						end = int(row["Ending Position"])
						break
			epi2regions[normalized] = [start, end]
		epi_regions.sort(key=lambda x: x[0])  # sort ascending by starting position
		return epi_regions, epi2regions

	def find_epitope_external_links(self, iedb_assay, normalized2unique):
		"""
		Find external links for each unique epitope

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		normalized2unique : dict of str : str
			dictionary mapping normalized epitopes to unique sequence

		Returns
		-------
		dict of str : list of str
			dictionary with key the normalized epitope and the value the list of external links
		"""
		print("Extract external links for each unique epitope")
		epi2external_links = {}
		for normalized, unique in normalized2unique.items():
			epi2external_links[normalized] = []
			for _, row in iedb_assay.loc[iedb_assay["Description"] == unique, ["Reference IRI"]].iterrows():
				if str(row["Reference IRI"]) not in epi2external_links[normalized]:
					epi2external_links[normalized].append(str(row["Reference IRI"]))
		return epi2external_links

	def calculate_RF_score(self, idbe_assay, normalized2unique):
		"""
		Calculate RF score for T cell assays, using:
		RF = (r-sqrt(r))/t,
		where r is the # positive responding assays
		and t is the # of the total tested assays

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		uniq_epitopes : dict of str : str
			dictionary mapping normalized epitope sequence to unique sequence

		Returns
		-------
		dict of str : float
			dictionary mapping the normalized epitope to RF score value
		"""
		print("Map unique epitope sequence to RF score")
		rf_scores = {}
		for normalized, unique in normalized2unique.items():
			# sum r and t for all assays testing the current epitope
			t, r = 0.0, 0.0
			for _, row in idbe_assay.loc[idbe_assay["Description"] == unique, ["Number of Subjects Tested",
			                                                                   "Number of Subjects Responded"]].iterrows():
				if isnan(row["Number of Subjects Tested"]):
					t = t + 1.0
				else:
					t = t + row["Number of Subjects Tested"]
				if isnan(row["Number of Subjects Responded"]):
					r = r + 1.0
				else:
					r = r + row["Number of Subjects Responded"]
			if t == 0:  # no information of tested subjects from all occurences of the epitope, assign -1
				rf_score = -1.0
			else:
				rf_score = (r - sqrt(r)) / t
			rf_scores[normalized] = rf_score
		assert len(rf_scores) == len(list(
			normalized2unique.keys())), "AssertionError: number of epitopes is not equal to the number of RF scores"
		return rf_scores

	def process_Tcells(self, tcells_current_virus):
		"""
		Process T-cells assays for current virus

		Parameters
		----------
		tcells_current_virus : Pandas.DataFrame
			T-cells assays for current virus

		Returns
		-------
		None
		"""
		print("\nProcess T-cells")
		self.ncbi_iedb_not_equal.append("T-cells:")
		for protein_id, protein in self.current_virus_proteins.items():
			print("---")
			print("Process protein with uniprot id: {}".format(protein_id))
			iedb_uniprot_id = self.url_prefixes["uniprot"] + protein.get_uniprot_id()
			tcells_current_protein = tcells_current_virus.loc[
				tcells_current_virus["Parent Protein IRI"] == iedb_uniprot_id]

			if tcells_current_protein.shape[0] == 0:
				print("Could not match using uniprot id")
				print("2nd attempt: Match with protein name")
				tcells_current_protein = tcells_current_virus.loc[
					tcells_current_virus["Parent Protein"].str.split("[").str[
						0].str.strip().str.lower() == protein.get_name().lower()]
			print("Number of non unique epitopes for protein = {}".format(tcells_current_protein.shape[0]))

			if tcells_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# normalized epitopes
				normalized2unique = self.normalize_unique_epitope_sequences(tcells_current_protein)
				# map epitope to allele
				normalized2allele = self.map_epitope2allele(tcells_current_protein, normalized2unique, True)
				# find epitope regions
				epi_regions, normalized2regions = self.find_epitope_regions(tcells_current_protein, normalized2unique)
				# calculate RF score
				normalized2rf_score = self.calculate_RF_score(tcells_current_protein, normalized2unique)
				# extract external links per unique epitope
				normalized2external_links = self.find_epitope_external_links(tcells_current_protein, normalized2unique)

				# create an Epitope object for each identified IEDB epitope
				host_taxon_id = self.host_taxon_id
				is_imported = True
				prediction_process = "IEDB_import"
				is_linear = True  # default for T cells
				# get protein record
				protein_record = protein.get_record()
				for normalized, region in normalized2regions.items():
					external_links = normalized2external_links[normalized]
					# decrease by one the region start to convert 1-index to 0-index
					reg_start, reg_end = region[0] - 1, region[-1]
					ncbi_prot_epi = protein_record.seq[reg_start:reg_end]
					if ncbi_prot_epi != normalized:  # append discordant epitopes
						self.ncbi_iedb_not_equal.append(
							"iedb frag:{}, ncbi prot:{}, iedb epitope link(s): {}".format(normalized, ncbi_prot_epi,
							                                                              " ".join(external_links)))
					if normalized2rf_score[
						normalized] > -1:  # if tested subjects information is given, export epitope
						epitope = Epitope(self.current_virus_taxon_id, protein.get_ncbi_id(), host_taxon_id, "T cell",
						                  normalized2allele[normalized], normalized2rf_score[normalized],
						                  str(normalized),
						                  reg_start, reg_end,
						                  is_imported,
						                  external_links, prediction_process, is_linear)
						self.current_virus_epitopes.append(epitope)
				all_attributes = self.current_virus_epitopes[-1].get_all_attributes()
				print(all_attributes)
		print("====")
