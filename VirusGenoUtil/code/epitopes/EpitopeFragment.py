import itertools


class EpitopeFragment:
	"""
	Class to store information for a fragment of a discontinuous epitope
	"""
	new_id = itertools.count()

	def __init__(self, parent_epi_id, fragm_seq, fragm_start, fragm_stop):
		"""
		EpitopeFragment constructor

		Parameters
		----------
		parent_epi_id : str
			id of the parent epitope (whole epitope that includes the fragment)
		fragm_seq : str
			fragment sequence
		fragm_start : int
			fragment start
		fragm_stop : int
			fragment stop
		"""
		self.fragment_id = next(EpitopeFragment.new_id)
		self.parent_epi_id = parent_epi_id
		self.fragment_seq = fragm_seq
		self.fragment_start = fragm_start
		self.fragment_stop = fragm_stop

	def get_all_attributes(self):
		"""
		Return all the attributes of the epitope fragment

		Returns
		-------
		dict of str
			all epitope fragment attributes returned in a dictionary
		"""
		return {"fragment_id": self.fragment_id, "parent_epi_id": self.parent_epi_id,
		        "fragment_seq": self.fragment_seq, "fragment_start": self.fragment_start,
		        "fragment_stop": self.fragment_stop}
