import itertools
from VirusGenoUtil.code.epitopes.EpitopeFragment import EpitopeFragment


class Epitope:
	"""
	Epitope to store the most important information of an epitope
	"""
	# credits: https://stackoverflow.com/questions/1045344/how-do-you-create-an-incremental-id-in-a-python-class/54318273#54318273
	new_id = itertools.count()

	def __init__(self, virus_taxid, protein_ncbi_id, host_iri, host_name, host_ncbi_id, cell_type, mhc_restriction,
	             response_frequency_info,
	             region_seq, region_start, region_stop, external_links, prediction_process, is_linear):
		"""
		Epitope construstor

		Parameters
		----------
		virus_taxid  : str
			virus NCBI taxonomy id
		protein_ncbi_id : str
			protein (antigen) NCBI id
		host_iri : str
			host taxon IRI
		host_name : str
			host name
		host_ncbi_id : str
			host ncbi id
		cell_type : str
			cell type
		mhc_restriction : str
			MHC restriction info (for T-cells)
		response_frequency_info : dict of str : str, float
			response frequency of epitope and if it was identified by a positive or negative assay
		region_seq : str
			epitope region sequence
		region_start : int
			epitope region start
		region_stop : int
			epitope region stop
		external_links : list of str
			IEDB reference links
		prediction_process : str
			epitope identification process
		is_linear : bool
			is a linear epitope (True) otherwise is discontinuous (False)
		"""
		self.id = next(Epitope.new_id)
		self.virus_taxid = virus_taxid
		self.protein_ncbi_id = protein_ncbi_id
		self.host_iri = host_iri
		self.host_name = host_name
		self.host_ncbi_id = host_ncbi_id
		self.cell_type = cell_type
		self.mhc_class = mhc_restriction['class']
		self.mhc_allele = mhc_restriction['allele']
		if response_frequency_info["positive"]["rf_score"] == -1:
			self.response_frequency_positive = None
		else:
			self.response_frequency_positive = str(response_frequency_info["positive"]["rf_score"])

		# convert the two boolean to one categorical
		if response_frequency_info["positive"]["exists_pos_assay"] and response_frequency_info["negative"]["exists_neg_assay"]:
			self.assay_types = "both"
		elif response_frequency_info["positive"]:
			self.assay_types = "positive"
		elif response_frequency_info["negative"]:
			self.assay_type = "negative"

		self.seq = region_seq
		self.region_start = region_start
		self.region_stop = region_stop
		self.prediction_process = prediction_process
		self.external_links = external_links
		self.is_linear = is_linear
		self.epitope_fragments = []
		if not is_linear:
			self.fragment()

	def fragment(self):
		"""
		Fragment discontinuous epitope into fragments

		Returns
		-------
		None
		"""
		print("Fragment discontinuous epitope: {}".format(self.seq))
		current_fragment_aa, current_fragment_pos = [], []
		# get as very first previous position the first position of the discontinuous epitope
		previous_position = int(self.seq.split(",")[0].strip()[1:]) - 1
		for aa_pos in self.seq.strip().split(","):
			aa_pos = aa_pos.strip()
			if aa_pos != "":  # if pos contains amino-acid, process it
				aa, pos = aa_pos[0], int(aa_pos[1:len(aa_pos)])
				if pos == previous_position + 1:  # continue current fragment
					current_fragment_aa.append(aa)
					current_fragment_pos.append(pos)
				elif self.seq.count(",") > 1:
					# current fragment is bigger than one amino-acid
					# current fragment has just finished
					assert len(current_fragment_pos) == len(
						current_fragment_aa), "AssertionError: identified fragment does not contain equal number of amino-acids and amino-acids positions"
					epi_fragment = EpitopeFragment(self.id, ''.join(current_fragment_aa), current_fragment_pos[0],
					                               current_fragment_pos[-1])
					self.epitope_fragments.append(epi_fragment)
					# start up the new fragment
					current_fragment_aa = [aa]
					current_fragment_pos = [pos]
				previous_position = pos

		# create the last fragment
		assert len(current_fragment_pos) == len(
			current_fragment_aa), "AssertionError: identified fragment does not contain equal number of amino-acids and amino-acids positions"
		epi_fragment = EpitopeFragment(self.id, ''.join(current_fragment_aa), current_fragment_pos[0],
		                               current_fragment_pos[-1])
		self.epitope_fragments.append(epi_fragment)

	def get_fragments(self):
		"""
		Get epitope fragments
		Returns
		-------
		list of EpitopeFragment
			list of epitope fragments for discontinuous epitope
		"""
		return self.epitope_fragments

	def get_all_attributes(self):
		"""
		Return all the attributes of the epitope

		Returns
		-------
		dict of str
			all epitope attributes returned in a dictionary
		"""
		return {"epitope_id": self.id,
		        "virus_taxid": self.virus_taxid,
		        "protein_ncbi_id": self.protein_ncbi_id,
		        "host_iri": self.host_iri,
		        "host_name": self.host_name,
		        "host_ncbi_id": self.host_ncbi_id,
		        "cell_type": self.cell_type,
		        "mhc_class": self.mhc_class,
		        "mhc_allele": self.mhc_allele,
		        "response_frequency_positive": self.response_frequency_positive,
		        "assay_types": self.assay_types,
		        "region_seq": self.seq,
		        "region_start": self.region_start,
		        "region_stop": self.region_stop,
		        "external_links": self.external_links,
		        "prediction_process": self.prediction_process,
		        "fragments": self.epitope_fragments,
		        "is_linear": self.is_linear
		        }
