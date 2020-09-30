from code.utils import is_fasta_file_extension, create_dir
from code.epitopes.Protein import Protein
from code.epitopes.Epitope import Epitope
from code.epitopes.EpitopeFragment import EpitopeFragment
from os.path import join, splitext, isfile, exists, isdir
from pathlib import Path
from os import scandir
from pandas import read_csv, unique, concat
from math import isnan, sqrt
import requests
import time
from Bio import Entrez


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

	def __init__(self, data_path, cell_epitopes_folder, viruses_folder, ontie_download_folder, host_taxon_id, host_name,
	             assay_type,
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
		ontie_download_folder : str
			ontie downloads folder containing ttls file
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

		self.tcell_iedb_assays, self.bcell_iedb_assays, self.mhc_iedb_assays = None, None, None
		self.host_taxon_id = host_taxon_id.split("_")[1]
		self.host_name = str(host_name)
		self.hosts_info = {}  # keep all host information in a dictionary: {host_iri: {'name': Homo Sapiens,'ncbi_id': 9606}, ..,}
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

		# create directory to save ONTIE.ttl files to map ONTIE ids to ncbi taxon ids
		self.ontie_downloads_path = join(data_path, ontie_download_folder)
		create_dir(self.ontie_downloads_path)  # create if you can access ONTIE efficiently

		# load email for Entrez usage
		self.entrez_email = open(join(data_path, "email.txt")).read().strip()
		Entrez.email = self.entrez_email
		print("Load email for Entrez access")

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
			if not epi_seq:
				epi_seq = ''
			if not epi_attributes["mhc_class"]:
				epi_attributes["mhc_class"] = ''
			if not epi_attributes["mhc_allele"]:
				epi_attributes["mhc_allele"] = ''

			epitope = tuple([int(epi_attributes["epitope_id"]), int(epi_attributes["virus_taxid"]),
			                 str(epi_attributes["host_iri"]), str(epi_attributes['host_name']), int(epi_attributes['host_ncbi_id']),
			                 str(epi_attributes["protein_ncbi_id"]), str(epi_attributes["cell_type"]),
			                 str(epi_attributes["mhc_class"]), str(epi_attributes["mhc_allele"]),
			                 str(epi_attributes["response_frequency_positive"]),
			                 str(epi_attributes["assay_types"]),
			                 epi_seq, int(epi_attributes["region_start"]), int(epi_attributes["region_stop"]),
			                 ",".join(epi_attributes["external_links"]),
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
					"epitope_id\tvirus_taxid\tsource_host_iri\tsource_host_name\thost_ncbi_id\tprotein_ncbi_id\tcell_type\tmhc_class\tmhc_restriction\tresponse_frequency_pos\tassay_types\tepitope_sequence\tepitope_start\tepitope_stop\texternal_links\tprediction_process\tis_linear\n")

		with open(join(self.output_path, epitopes_name), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			for epitope in self.current_virus_epitopes:
				print("Write IEDB imported epitope")
				epi_attributes = epitope.get_all_attributes()

				if epi_attributes["is_linear"]:
					epi_seq = epi_attributes["region_seq"]
				else:
					epi_seq = None
				if not epi_seq:
					epi_seq = ''
				if not epi_attributes["mhc_class"]:
					epi_attributes["mhc_class"] = ''
				if not epi_attributes["mhc_allele"]:
					epi_attributes["mhc_allele"] = ''

				epitope_row = "\t".join(
					[str(epi_attributes["epitope_id"]), epi_attributes["virus_taxid"],
					 epi_attributes["host_iri"], epi_attributes['host_name'], epi_attributes['host_ncbi_id'],
					 epi_attributes["protein_ncbi_id"], epi_attributes["cell_type"],
					 epi_attributes["mhc_class"], epi_attributes["mhc_allele"],
					 str(epi_attributes["response_frequency_positive"]),
					 str(epi_attributes["assay_types"]), epi_seq,
					 str(epi_attributes["region_start"]), str(epi_attributes["region_stop"]),
					 ",".join(epi_attributes["external_links"]),
					 str(epi_attributes["prediction_process"]), str(epi_attributes["is_linear"])])
				epitopes_out.write(epitope_row + "\n")
		print("====")

	def load_iedb_csvs(self):
		"""
		Load IEDB for B, T cells and MHC ligands assays csvs
		Load files in chunks, credits: https://stackoverflow.com/questions/11622652/large-persistent-dataframe-in-pandas
		Returns
		-------
		None
		"""
		assert isfile(join(self.cell_epitopes_path,
		                   "tcell_full_v3.csv.gz")), "AssertionError: IEDB Tcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		# tcell_text_file_reader = read_csv(join(self.cell_epitopes_path, "tcell_full_v3.csv.gz"), sep=",", header=1,
		#                                   compression='gzip', iterator=True, chunksize=1000)
		# self.tcell_iedb_assays = concat(tcell_text_file_reader,ignore_index=True)
		self.tcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "tcell_full_v3.csv.gz"), sep=",", header=1,
		                                  compression='gzip')
		assert isfile(join(self.cell_epitopes_path,
		                   "bcell_full_v3.csv.gz")), "AssertionError: IEDB Bcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		# bcell_text_file_reader = read_csv(join(self.cell_epitopes_path, "bcell_full_v3.csv.gz"), sep=",", header=1,
		#                                   compression='gzip', iterator=True, chunksize=1000)
		# self.bcell_iedb_assays = concat(bcell_text_file_reader,ignore_index=True)
		self.bcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "bcell_full_v3.csv.gz"), sep=",", header=1,
		                                  compression='gzip')
		assert isfile(join(self.cell_epitopes_path,
		                   "mhc_ligand_full.csv.gz")), "AssertionError: IEDB MHC ligand assays csv was not found in {}".format(
			self.cell_epitopes_path)
		# mhc_iedb_assays = read_csv(join(self.cell_epitopes_path, "mhc_ligand_full.csv.gz"), sep=",", header=1,
		#                            compression='gzip', iterator=True, chunksize=10000)
		# self.mhc_iedb_assays = concat(mhc_iedb_assays, ignore_index=True)
		self.mhc_iedb_assays = read_csv(join(self.cell_epitopes_path, "mhc_ligand_full.csv.gz"), sep=",", header=1,
		                                compression='gzip', nrows=1000)
		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))
		print("T cell subset number of non-unique epitopes: {}".format(self.tcell_iedb_assays.shape[0]))
		print("MHC ligand subset number of non-unique epitopes: {}".format(self.mhc_iedb_assays.shape[0]))
		print("---")

	def subset_iedb_by_host_iri(self, iedb_assay, host_iri):
		"""
		Subset IEBD assay dataframe by host IRI

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			iedb assay data
		host_iri : str
			host IRI used to subset the iedb data
		Returns
		-------
		Pandas.DataFrame
			subset of IEDB assay
		"""
		print("\n ### ### New host ### ### ")
		print("Select epitopes found only for host iri:{}".format(host_iri))
		return iedb_assay.loc[iedb_assay["Host IRI"] == host_iri]

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
			self.mhc_iedb_assays = self.mhc_iedb_assays.loc[
				self.mhc_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]
		else:
			print("Select all assay type: positive and negative")
		if self.host_taxon_id != "all":
			print("Select epitopes experimental identified in host id: {} = {}".format(self.host_taxon_id,
			                                                                           self.host_name))
			host_taxid = self.url_prefixes["taxid"] + self.host_taxon_id
			self.tcell_iedb_assays = self.tcell_iedb_assays.loc[self.tcell_iedb_assays["Host IRI"] == host_taxid]
			self.bcell_iedb_assays = self.bcell_iedb_assays.loc[self.bcell_iedb_assays["Host IRI"] == host_taxid]
			self.mhc_iedb_assays = self.mhc_iedb_assays.loc[self.mhc_iedb_assays["Host IRI"] == host_taxid]
		else:
			print("Select all available hosts")

		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))
		print("T cell subset number of non-unique epitopes: {}".format(self.tcell_iedb_assays.shape[0]))
		print("MHC ligand subset number of non-unique epitopes: {}".format(self.mhc_iedb_assays.shape[0]))

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
		Pandas.DataFrame
			B-cell IEDB assays only for virus id
		Pandas.DataFrame
			T-cell IEDB assay only for virus id
		Pandas.DataFrame
			MHC ligand assay only for virus id
		"""
		print("Get iedb only for taxon id={}".format(self.current_virus_taxon_id))
		# select T cell assay by virus id
		tcell_iedb_virus = self.tcell_iedb_assays.loc[
			self.tcell_iedb_assays["Parent Species IRI"].str.split("NCBITaxon_").str[
				-1].str.strip() == self.current_virus_taxon_id]
		if tcell_iedb_virus.shape[0] == 0:
			print("2nd attempt: Match with Organism IRI")
			tcell_iedb_virus = self.tcell_iedb_assays.loc[
				self.tcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[
					-1].str.strip() == self.current_virus_taxon_id]

		# select B cell assays by virus id
		bcell_iedb_virus = self.bcell_iedb_assays.loc[
			self.bcell_iedb_assays["Parent Species IRI"].str.split("NCBITaxon_").str[
				-1].str.strip() == self.current_virus_taxon_id]
		if bcell_iedb_virus.shape[0] == 0:
			print("2nd attempt: Match with Organism IRI")
			bcell_iedb_virus = self.bcell_iedb_assays.loc[
				self.bcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[
					-1].str.strip() == self.current_virus_taxon_id]

		# select MHC ligand assays by virus id
		mhc_iedb_virus = self.mhc_iedb_assays.loc[
			self.mhc_iedb_assays["Parent Species IRI"].str.split("NCBITaxon_").str[
				-1].str.strip() == self.current_virus_taxon_id]
		if mhc_iedb_virus.shape[0] == 0:
			print("2nd attempt: Match with Organism IRI")
			mhc_iedb_virus = self.mhc_iedb_assays.loc[
				self.mhc_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[
					-1].str.strip() == self.current_virus_taxon_id]
		print("B cell subset for virus contains {} non unique epitopes".format(bcell_iedb_virus.shape[0]))
		print("B cell head: {}".format(bcell_iedb_virus.head()))
		print("T cell subset for virus contains {} non unique epitopes".format(tcell_iedb_virus.shape[0]))
		print("T cell head: {}".format(tcell_iedb_virus.head()))
		print("MHC ligand subset for virus contains {} non unique epitopes".format(mhc_iedb_virus.shape[0]))
		print("MHC ligand head: {}".format(mhc_iedb_virus.head()))
		return bcell_iedb_virus, tcell_iedb_virus, mhc_iedb_virus

	@staticmethod
	def concatenate_host_name_ncbi(host_names, host_ncbis):
		"""
		Concatenate host name and ncbi ids into one dictionary

		Parameters
		----------
		host_names : dict of str: str
			map host iri to name
		host_ncbis : dict of str: str
			map host iri to ncbi

		Returns
		-------
		dict of dict of str: str
			Concatenated dictionary of host information
			example: {'http://purl.obolibrary.org/obo/NCBITaxon_9606'{'name':'Homo Sapiens', 'ncbi_id':'9606'}}
		"""
		assert set(list(host_names.keys())) == set(list(host_ncbis.keys())), "AssertionError: host iris are not equal"
		host_info = {}
		for iri, name in host_names.items():
			host_info[iri] = {'name': name, 'ncbi_id': host_ncbis[iri]}
		return host_info

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
		if self.host_taxon_id == "all":
			host_name_iri_tcells = self.tcell_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_bcells = self.bcell_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_mhc = self.mhc_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_all_cell_types = concat([host_name_iri_tcells, host_name_iri_bcells, host_name_iri_mhc])
			# host_name_iri_all_cell_types.to_csv(join(self.output_path, "all_cells_host_info.csv"))
			self.extract_host_info(host_name_iri_all_cell_types)

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
		if self.host_taxon_id == "all":
			host_name_iri_tcells = self.tcell_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_bcells = self.bcell_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_mhc = self.mhc_iedb_assays.loc[:, ['Name', 'Host IRI']]
			host_name_iri_all_cell_types = concat([host_name_iri_tcells, host_name_iri_bcells, host_name_iri_mhc])
			# host_name_iri_all_cell_types.to_csv(join(self.output_path, "all_cells_host_info.csv"))
			self.extract_host_info(host_name_iri_all_cell_types)

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

	def retrieve_ncbi4ontie(self, ontie_id):
		"""
		Recursive function to retrieve ncbi id from ontie id using the IEDB ontology API:
		https://ontology.iedb.org/doc/api.html

		Credits: Tommaso for his help!
		Parameters
		----------
		ontie_id : str
			ontie id in the form: "ONTIE_000001"

		Returns
		-------
		ncbi_id : str
			ncbi id relating to the very first input ontie id
		"""
		print("Retrieve ncbi id for ontie id")
		if not isfile(join(self.ontie_downloads_path, ontie_id + ".ttl")):
			print("ttl is not found, retrieve it from ONTIE")
			# if ontie file is not already downloaded,
			# download through the API
			print("https://ontology.iedb.org/ontology/" + ontie_id + "?format=ttl")
			for i in range(0, 10):
				response = requests.get("https://ontology.iedb.org/ontology/" + ontie_id + "?format=ttl",
				                        time.sleep(20))
				print("response code: {}".format(response.status_code))
				if response.status_code == 200:
					print("Success response code!")
					break
			with open(join(self.ontie_downloads_path, ontie_id + ".ttl"), 'wb') as f:
				f.write(response.content)
		else:  # if ontie file is downloaded, then go on processing file
			print('Already downloaded ONTIE file')

		# process ontie file to get subclass information
		with open(join(self.ontie_downloads_path, ontie_id + ".ttl"), 'r') as f:
			subclass_lines = []
			for line in f.readlines():
				clean_line = line.strip(' ;\n')
				if clean_line[0:15] == 'rdfs:subClassOf':
					subclass_lines.append(clean_line)

		# apply depth-first search: work only with the one branch in ontie ontology graph (subclass edge):
		subclass_info = subclass_lines[0]
		if 'obo:NCBITaxon_' in subclass_info:  # if ncbi id was found, terminate recursion
			return subclass_info.split('obo:NCBITaxon_')[1]
		else:  # if ncbi id is not found in the subclass, continue the recursion using as ontie the current subclass id
			ontie_id = "ONTIE_" + subclass_info.split('ONTIE:')[1]
			return self.retrieve_ncbi4ontie(ontie_id)

	# TODO: function to make ncbi taxon -> ncbi id, ONTIE id -> ncbi id
	def map_host_iri2ncbi(self, iedb_assay):
		"""
		Map each host iri to ncbi id

		Parameters
		----------
		iedb_assay : Pandas.DataFrame

		Returns
		-------
		dict of str: str
			map host iri to ncbi id
		"""
		host_iri2ncbi = {}
		for host_iri in unique(iedb_assay["Host IRI"]):
			print("new host iri")
			ncbi_id = "unknown"
			if "NCBITaxon" in str(host_iri):  # extract existing ncbi taxon id
				ncbi_id = str(host_iri).split("NCBITaxon_")[1]
			elif "ONTIE" in str(host_iri):
				ontie_id = "ONTIE_" + str(host_iri).split("ONTIE_")[1]
				ncbi_id = self.retrieve_ncbi4ontie(ontie_id)
			host_iri2ncbi[host_iri] = ncbi_id
		assert "unknown" not in list(
			host_iri2ncbi.values()), "AssertionError: exists host iri which neither ncbi, nor ontie"
		return host_iri2ncbi

	def find_unique_host_iris(self, iedb_assay):
		"""
		Find all unique host IRIs from the input data frame

		Parameters
		----------
		iedb_assay: Pandas.DataFrame
			dataframe to extract all unique host taxon id

		Returns
		-------
		list of str
			list of all unique host IRIs
		"""
		all_unique_host_iris = []
		for host_iri in unique(iedb_assay["Host IRI"]):
			all_unique_host_iris.append(host_iri)
		return all_unique_host_iris

	def map_host_iri2name(self, iedb_assay):
		"""
		Extract host names for each IRI, IRI -> name

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe to extract host information from

		Returns
		-------
		dict of str: str
			host information
		"""
		host_iri2names = {}
		for host_iri in unique(iedb_assay["Host IRI"]):
			for _, name in iedb_assay.loc[iedb_assay["Host IRI"] == host_iri, "Name"].iteritems():
				if host_iri not in host_iri2names:
					host_iri2names[host_iri] = [str(name)]
				else:
					if str(name) not in host_iri2names[host_iri]:
						host_iri2names[host_iri].append(str(name))
		# sort and get the first name out of all unique available ones
		for host_iri in list(host_iri2names.keys()):
			host_iri2names[host_iri] = sorted(host_iri2names[host_iri])[0]
		return host_iri2names

	def extract_host_info(self, iedb_assay):
		"""
		Extract host information for each IRI, IRI -> name, ncbi id

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe to extract host information from

		"""
		print("Extract host info")
		for host_iri in unique(iedb_assay["Host IRI"]):
			for _, name in iedb_assay.loc[iedb_assay["Host IRI"] == host_iri, "Name"].iteritems():
				if host_iri not in self.hosts_info:
					self.hosts_info[host_iri] = {'name': [], 'ncbi_id': []}
					self.hosts_info[host_iri]['name'].append(str(name))
				else:
					if str(name) not in self.hosts_info[host_iri]['name']:
						self.hosts_info[host_iri]['name'].append(str(name))
		# sort and get the first name
		# for this first name get the taxonomy ncbi id
		for host_iri in list(self.hosts_info.keys()):
			host_name = sorted(self.hosts_info[host_iri]['name'])[0]
			host_id = self.extract_host_ncbi_id(host_iri, host_name)
			self.hosts_info[host_iri]['name'] = host_name
			self.hosts_info[host_iri]['ncbi_id'] = host_id

	def retrieve_ncbi_id_from_name(self, host_name):
		"""
		Retrieve ncbi id from host name

		Parameters
		----------
		host_name : str
			host name

		Returns
		-------
		str
			host ncbi id
		"""
		print("Retrieve ncbi id from host name: {}".format(host_name))
		name_split = host_name.strip().split()
		if len(name_split) > 1:
			species_name = "+".join([name_split[0], name_split[1]]).strip()
		else:
			species_name = name_split[0]
		search = Entrez.esearch(term=species_name, db="taxonomy", retmode="xml")
		record = Entrez.read(search)
		if len(record['IdList']) > 0:
			return record['IdList'][0]
		else:
			return "unknown"

	def extract_host_ncbi_id(self, host_iri, host_name):
		"""
		Extract host ncbi id using the name of the host iri
		Parameters
		----------
		host_iri : str
			host iri
		host_name : str
			host name

		Returns
		-------
		str
			host ncbi id
		"""
		print("Extract host ncbi id")
		# first try to see if you have ncbi id in the host iri
		if "NCBITaxon" in str(host_iri):  # extract existing ncbi taxon id
			ncbi_id = str(host_iri).split("NCBITaxon_")[1]
		else:
			# second try retrieving ncbi from the host name
			response_id = self.retrieve_ncbi_id_from_name(host_name)
			if response_id != "unknown":
				ncbi_id = response_id
			else:  # third try: find ONTIE from web API
				assert "ONTIE" in str(host_iri), "AssertError: ONTIE not in iri"
				ontie_id = "ONTIE_" + str(host_iri).split("ONTIE_")[1]
				ncbi_id = self.retrieve_ncbi4ontie(ontie_id)
		return ncbi_id

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
		bcells_current_virus, tcells_current_virus, mhc_current_virus = self.subset_iedb_by_virus_id()
		self.current_virus_epitopes = []  # clear current virus epitopes
		self.current_virus_epi_fragments = []  # clear current virus epitope fragments
		self.ncbi_iedb_not_equal.append("=== Virus taxid={} ===".format(self.current_virus_taxon_id))
		# bcells_current_virus.to_csv(join(self.output_path, "bcells_virus_" + self.current_virus_taxon_id + ".csv"))
		# tcells_current_virus.to_csv(join(self.output_path, "tcells_virus_" + self.current_virus_taxon_id + ".csv"))
		# mhc_current_virus.to_csv(join(self.output_path, "mhc_virus_" + self.current_virus_taxon_id + ".csv"))

		if self.host_taxon_id == "all":
			available_hosts_bcells = self.find_unique_host_iris(bcells_current_virus)
			available_hosts_tcells = self.find_unique_host_iris(tcells_current_virus)
			available_hosts_mhc = self.find_unique_host_iris(mhc_current_virus)
			for unique_host_iri in list(self.hosts_info.keys()):
				if unique_host_iri in available_hosts_bcells:
					bcells_current_virus_host = self.subset_iedb_by_host_iri(bcells_current_virus, unique_host_iri)
					self.process_Bcells(bcells_current_virus_host)
				if unique_host_iri in available_hosts_tcells:
					tcells_current_virus_host = self.subset_iedb_by_host_iri(tcells_current_virus, unique_host_iri)
					self.process_Tcells(tcells_current_virus_host)
				if unique_host_iri in available_hosts_mhc:
					mhc_current_virus_host = self.subset_iedb_by_host_iri(mhc_current_virus, unique_host_iri)
					self.process_MHC(mhc_current_virus_host)
		else:  # for only one host id the assay is already subset
			self.process_Bcells(bcells_current_virus)
			self.process_Tcells(tcells_current_virus)
			self.process_MHC(mhc_current_virus)

	def process_MHC(self, mhc_current_virus):
		"""
		Process MHC ligand assays for current virus

		Parameters
		----------
		mhc_current_virus : Pandas.DataFrame
			MHC ligand assays for current virus

		Returns
		-------
		None
		"""
		print("\nProcess MHC ligands")
		self.ncbi_iedb_not_equal.append("MHC:")
		for protein_id, protein in self.current_virus_proteins.items():
			print("---")
			print("Process protein with uniprot id: {}".format(protein_id))
			iedb_uniprot_id = self.url_prefixes["uniprot"] + protein.get_uniprot_id()
			mhc_current_protein = mhc_current_virus.loc[
				mhc_current_virus["Parent Protein IRI"] == iedb_uniprot_id]

			if mhc_current_protein.shape[0] == 0:
				print("Could not match using uniprot id")
				print("2nd attempt: Match with protein name")
				mhc_current_protein = mhc_current_virus.loc[
					mhc_current_virus["Parent Protein"].str.split("[").str[
						0].str.strip().str.lower() == protein.get_name().lower()]

			if mhc_current_protein.shape[0] == 0:
				print("Could not match using protein name")
				print("3rd attempt: Match with antigen name")
				mhc_current_protein = mhc_current_virus.loc[
					mhc_current_virus["Antigen Name"].str.strip().str.lower()
					== protein.get_name().lower()]
			print("Number of non unique epitopes for protein = {}".format(mhc_current_protein.shape[0]))

			if mhc_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# normalized epitopes
				normalized2unique = self.normalize_unique_epitope_sequences(mhc_current_protein)
				# map epitope to allele
				normalized2allele = self.map_epitope2MCH_info(mhc_current_protein, normalized2unique, False)
				# find epitope regions
				epi_regions, normalized2regions = self.find_epitope_regions(mhc_current_protein, normalized2unique)
				# calculate RF info
				normalized2rf_info = self.calculate_RF_info(mhc_current_protein, normalized2unique)
				# extract external links per unique epitope
				normalized2external_links = self.find_epitope_external_links(mhc_current_protein, normalized2unique)

				# get current host info: IRI, name and ncbi id
				host_iri = unique(mhc_current_protein["Host IRI"])[0]
				host_name = self.hosts_info[host_iri]['name']
				host_ncbi_id = self.hosts_info[host_iri]['ncbi_id']

				# create an Epitope object for each identified IEDB epitope
				host_taxon_id = self.host_taxon_id
				prediction_process = "IEDB_import"
				is_linear = True  # default for T cells
				# get protein record
				protein_record = protein.get_record()
				for normalized, region in normalized2regions.items():
					external_links = normalized2external_links[normalized]
					reg_start, reg_end = region[0], region[-1]
					if not (
							reg_start == reg_end and reg_start == 0):  # if linear epitopes has known start and end positions
						# decrease by one the region start to convert 1-index to 0-index
						ncbi_prot_epi = protein_record.seq[region[0] - 1:region[-1]]
						if ncbi_prot_epi != normalized:  # append discordant epitopes
							self.ncbi_iedb_not_equal.append(
								"iedb frag:{}, ncbi prot:{}, iedb epitope link(s): {}".format(normalized, ncbi_prot_epi,
								                                                              " ".join(external_links)))

						# create epitope object for the two possible types of assay, if are found in the data
						epitope = Epitope(self.current_virus_taxon_id, protein.get_ncbi_id(), host_iri, host_name,
						                  host_ncbi_id, "MHC Ligand",
						                  normalized2allele[normalized], normalized2rf_info[normalized],
						                  str(normalized),
						                  reg_start, reg_end,
						                  external_links, prediction_process, is_linear)
						self.current_virus_epitopes.append(epitope)
					all_attributes = self.current_virus_epitopes[-1].get_all_attributes()
					print(all_attributes)
		print("====")

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
				print(protein.get_name().lower())
				bcells_current_protein = bcells_current_virus.loc[
					bcells_current_virus["Parent Protein"].str.split("[").str[
						0].str.strip().str.lower() == protein.get_name().lower()]
			if bcells_current_protein.shape[0] == 0:
				print("Could not match using parent protein name")
				print("3rd attempt: Match with antigen name")
				bcells_current_protein = bcells_current_virus.loc[
					bcells_current_virus["Antigen Name"].str.strip().str.lower() == protein.get_name().lower()]
			print("Number of non unique epitopes for protein = {}".format(bcells_current_protein.shape[0]))

			if bcells_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# normalize epitopes
				normalized2unique = self.normalize_unique_epitope_sequences(bcells_current_protein)
				# map epitope to allele
				normalized2allele = self.map_epitope2MCH_info(bcells_current_protein, normalized2unique, True)
				# find epitope regions
				epi_regions, normalized2regions = self.find_epitope_regions(bcells_current_protein, normalized2unique)
				# calculate RF score
				normalized2rf_info = self.calculate_RF_info(bcells_current_protein, normalized2unique)
				# extract external links per unique epitope
				normalized2external_links = self.find_epitope_external_links(bcells_current_protein, normalized2unique)

				# get current host info: IRI, name and ncbi id
				host_iri = unique(bcells_current_protein["Host IRI"])[0]
				host_name = self.hosts_info[host_iri]['name']
				host_ncbi_id = self.hosts_info[host_iri]['ncbi_id']

				# create an Epitope object for each identified IEDB epitope
				prediction_process = "IEDB_import"
				# get protein record
				protein_record = protein.get_record()
				for normalized, region in normalized2regions.items():
					if region[0] == -1 and region[-1] == -1:  # discontinuous epitope
						is_linear = False
						reg_start, reg_end = self.get_discontinous_epi_start_stop(normalized)
					elif not (region[0] == region[-1] and region[
						0] == 0):  # linear epitope with known start and end positions
						is_linear = True
						reg_start, reg_end = region[0], region[-1]
						# decrease by one the region start to convert 1-index to 0-index
						ncbi_prot_epi = protein_record.seq[region[0] - 1: region[-1]]

						if ncbi_prot_epi != normalized:  # append discordant epitopes
							self.ncbi_iedb_not_equal.append(
								"iedb frag:{}, ncbi prot:{}, iedb epitope link(s): {}".format(normalized, ncbi_prot_epi,
								                                                              " ".join(
									                                                              normalized2external_links[
										                                                              normalized])))
					else:
						print("{} is a linear epitope without start and end".format(normalized))
						continue
					external_links = normalized2external_links[normalized]
					epitope = Epitope(self.current_virus_taxon_id, protein.get_ncbi_id(), host_iri,
					                  host_name, host_ncbi_id, "B cell",
					                  normalized2allele[normalized], normalized2rf_info[normalized],
					                  str(normalized),
					                  reg_start, reg_end,
					                  external_links, prediction_process, is_linear)
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
		print("Get discontinuous epitope start stop")
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

	def map_epitope2MCH_info(self, iedb_assay, normalized2unique, is_bcell):
		"""
		Map epitopes to allele information (MHC class and allele name)

		Parameters
		----------
		iedb_assay : Pandas.DataFrame
			dataframe created from csv of IEDB assay tab
		normalized2unique : dict of str : str
			dictionary mapping normalized epitope to unique epitope
		is_bcell : bool
			Processing B-cells (True), otherwise processing B-cells or MHC ligand assays (False)

		Returns
		-------
		dict of str: dict of str : str
			dictionary mapping normalzed epitopes to their allele info
			example: {.."MNTP..FGRQW":{'MHC_class':'I','MHC_allele':'H2-Kd'} ..}

		Defaults:
		variable        default_value       meaning
		class               None            not applicable (B-cell are not related to MHC class)
		class               unknown         not reported MHC class (Tcell or MHC ligand)
		allele              None            not applicable (B-cell are not related to MHC allele)
		allele              unknown         not reported allele name (Tcell or MHC ligand)
		"""
		print("Map unique epitope sequence to allele information")
		normalized2mhc = {}
		# for each pair normalized,unique find the allele information based on the unique allele
		# then save it for the corresponding normalized epitope
		for normalized_epi, uniq_epi in normalized2unique.items():
			if is_bcell:
				normalized2mhc[normalized_epi] = {"class": None, "allele": None}
			else:
				mhc_class_name = None
				if "Class" in list(iedb_assay.columns):  # T-cells
					mhc_class_name = "Class"
				else:  # MHC ligand
					mhc_class_name = "MHC allele class"
				if uniq_epi not in normalized2mhc:
					all_mhc_allele_names = []
					all_mhc_classes = []
					for _, row in iedb_assay.loc[
						iedb_assay["Description"] == uniq_epi, ["Allele Name", mhc_class_name]].iterrows():
						mhc_class = str(row[mhc_class_name])
						mhc_allele = str(row["Allele Name"])
						if mhc_allele != 'nan':
							if mhc_allele.replace(" ", "_") not in all_mhc_allele_names:
								all_mhc_allele_names.append(mhc_allele.replace(" ", "_"))
						if mhc_class != 'nan' and mhc_class not in all_mhc_classes:
							all_mhc_classes.append(mhc_class)
					# after loop, add default values for class and allele
					normalized2mhc[normalized_epi] = {"class": "unknown", "allele": "unknown"}
					# update defaults, if information is present
					if len(all_mhc_classes) > 0:
						normalized2mhc[normalized_epi]["class"] = ",".join(all_mhc_classes)
					if len(all_mhc_allele_names) > 0:
						normalized2mhc[normalized_epi]["allele"] = ",".join(all_mhc_allele_names)

		return normalized2mhc

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
				start, end = 0, 0  # default start end for linear epitopes with unknown start and end positions in the antigenic sequence
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

	def calculate_RF_info(self, idbe_assay, normalized2unique):
		"""
		Calculate RF score for assays, using:
		RF = (r-sqrt(r))/t,
		where r is the # positive responding assays
		and t is the # of the total tested assays
		and also extract if the assay is positive or negative

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		uniq_epitopes : dict of str : str
			dictionary mapping normalized epitope sequence to unique sequence

		Returns
		-------
		dict of dict of
			dictionary mapping the normalized epitope to RF score value and if the assay is positive or negative
			example: {.. {"MNTP..AAQYF":{ 'positive':{'rf_score':0.5,'exists_pos_assay':True},'negative':{'rf_score':0.0,'exists_neg_assay':False} } }..}

		Defaults:
		variable        default_value       meaning
		rf_score            0               total rf_score is 0
		rf_score            -1              no reported necessary information to compute rf_score
		=> -1 will be converted to None in the Epitope table attribute value
		"""
		print("Map unique epitope sequence to RF score")
		rf_info = {}
		for normalized, unique in normalized2unique.items():
			# sum r and t for all assays testing the current epitope
			t_positive, r_positive = 0.0, 0.0
			t_negative = 0.0
			epi_rf_info = {'positive': {'rf_score': -1.0, 'exists_pos_assay': False},
			               'negative': {'rf_score': -1.0, 'exists_neg_assay': False}}
			for _, row in idbe_assay.loc[
				idbe_assay["Description"] == unique, ["Qualitative Measure", "Number of Subjects Tested",
				                                      "Number of Subjects Responded"]].iterrows():
				if row["Qualitative Measure"].lower() == 'negative':
					epi_rf_info['negative']['exists_neg_assay'] = True
					t_negative = 0
				else:  # positive assays
					epi_rf_info['positive']['exists_pos_assay'] = True
					if isnan(row["Number of Subjects Tested"]):
						t_positive = t_positive + 1.0
					else:
						t_positive = t_positive + row["Number of Subjects Tested"]
					if isnan(row["Number of Subjects Responded"]):
						r_positive = r_positive + 1.0
					else:
						r_positive = r_positive + row["Number of Subjects Responded"]
			# after passing all occurences of an unique epitope
			# aggregate the RF for positive and examine if it was found in a negative assay
			if epi_rf_info['positive']['exists_pos_assay']:
				if t_positive == 0:  # no information of tested subjects from all occurrences of the epitope, assign -1
					rf_score = -1.0
				else:
					rf_score = (r_positive - sqrt(r_positive)) / t_positive
				epi_rf_info['positive']['rf_score'] = rf_score
			if epi_rf_info['negative'][
				'exists_neg_assay']:  # for negative assay, RF = -1 (RF is calculated for positive assays)
				# here use t_negative to compute rf_score, if needed
				epi_rf_info['negative']['rf_score'] = -1.0
			rf_info[normalized] = epi_rf_info
		assert len(rf_info) == len(list(
			normalized2unique.keys())), "AssertionError: number of epitopes is not equal to the number of RF scores"
		return rf_info

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

			if tcells_current_protein.shape[0] == 0:
				print("Could not match using protein name")
				print("3rd attempt: Match with antigen name")
				tcells_current_protein = tcells_current_virus.loc[
					tcells_current_virus["Antigen Name"].str.strip().str.lower()
					== protein.get_name().lower()]
			print("Number of non unique epitopes for protein = {}".format(tcells_current_protein.shape[0]))

			if tcells_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# normalized epitopes
				normalized2unique = self.normalize_unique_epitope_sequences(tcells_current_protein)
				# map epitope to allele
				normalized2allele = self.map_epitope2MCH_info(tcells_current_protein, normalized2unique, False)
				# find epitope regions
				epi_regions, normalized2regions = self.find_epitope_regions(tcells_current_protein, normalized2unique)
				# calculate RF info
				normalized2rf_info = self.calculate_RF_info(tcells_current_protein, normalized2unique)
				# extract external links per unique epitope
				normalized2external_links = self.find_epitope_external_links(tcells_current_protein, normalized2unique)

				# get current host info: IRI, name and ncbi id
				host_iri = unique(tcells_current_protein["Host IRI"])[0]
				host_name = self.hosts_info[host_iri]['name']
				host_ncbi_id = self.hosts_info[host_iri]['ncbi_id']

				# create an Epitope object for each identified IEDB epitope
				host_taxon_id = self.host_taxon_id
				prediction_process = "IEDB_import"
				is_linear = True  # default for T cells
				# get protein record
				protein_record = protein.get_record()
				for normalized, region in normalized2regions.items():
					external_links = normalized2external_links[normalized]
					reg_start, reg_end = region[0], region[-1]

					if not (
							reg_start == reg_end and reg_start == 0):  # if linear epitopes has known start and end positions
						# decrease by one the region start to convert 1-index to 0-index
						ncbi_prot_epi = protein_record.seq[region[0] - 1:region[-1]]
						if ncbi_prot_epi != normalized:  # append discordant epitopes
							self.ncbi_iedb_not_equal.append(
								"iedb frag:{}, ncbi prot:{}, iedb epitope link(s): {}".format(normalized, ncbi_prot_epi,
								                                                              " ".join(external_links)))

						# create epitope object for the two possible types of assay, if are found in the data
						epitope = Epitope(self.current_virus_taxon_id, protein.get_ncbi_id(), host_iri,
						                  host_name, host_ncbi_id, "T cell",
						                  normalized2allele[normalized], normalized2rf_info[normalized],
						                  str(normalized),
						                  reg_start, reg_end,
						                  external_links, prediction_process, is_linear)
						self.current_virus_epitopes.append(epitope)
				all_attributes = self.current_virus_epitopes[-1].get_all_attributes()
				print(all_attributes)
		print("====")
