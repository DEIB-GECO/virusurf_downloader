from os.path import join, splitext, exists


class BepiPredEpitopes:
	"""
	Class to extract epitopes found by BepiPred 2.0
	Following parameters from the reference paper

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: bepipred output tsv file of continous epitopes for the structural proteins of SARS-CoV-2
	(please update) output: file containing the aligned sequences, response frequency and identity
	"""

	def __init__(self, data_path, bepipred_folder, virus_ncbi_id, out_path):
		"""
		BepiPredEpitopes constructor

		Parameters
		----------
		data_path : str
			input data root
		bepipred_folder : str
			bepipred folder
		virus_ncbi_id : str
			virus NCBI id
		out_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.bepipred_folder = join(self.data_path, bepipred_folder)
		self.virus_ncbi_id = virus_ncbi_id
		self.out_path = out_path
		self.prediction_method = "Bepipred_2"

	@staticmethod
	def bepipred_tabs2info(bepipred_line):
		"""
		Convert bepipred prediction line to tabs

		Parameters
		----------
		bepipred_line : str
			bepipred prediction line

		Returns
		-------
		dict of str: str
			bepipred prediction information
		"""
		bepipred_tabs = bepipred_line.split("\t")
		assert len(
			bepipred_tabs) == 9, "AssertionError: Excepting 9 columns for bepipred prediciton line, got {}".format(
			len(bepipred_tabs))
		bepipred_info = {"entry": bepipred_tabs[0], "pos": bepipred_tabs[1], "AA": bepipred_tabs[2],
		                 "RSA": bepipred_tabs[3], "helix": bepipred_tabs[4], "sheet": bepipred_tabs[5],
		                 "coil": bepipred_tabs[6], "pred": bepipred_tabs[7], "epi": bepipred_tabs[8]}
		return bepipred_info

	def write_epitope_row(self, current_epitope, epitopes_filename):
		"""
		Write epitope row in output tsv file

		Parameters
		----------
		current_epitope : dict of str: str
			current epitope information
		epitopes_filename :
			epitope output tsv file

		Returns
		-------
		None
		"""
		with open(join(self.out_path, epitopes_filename), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_filename))
			epitope_row = "\t".join(
				[self.virus_ncbi_id, current_epitope["protein_id"], current_epitope["start"], current_epitope["end"],
				 current_epitope["seq"], "not applicable", self.prediction_method])
			epitopes_out.write(epitope_row + "\n")

	def save2tsv(self, bepipred_file):
		"""
		Save epitopes predicted by BepiPred to tsv file

		Parameters
		----------
		bepipred_file : str
			Bepipred prediction output file

		Returns
		-------
		None
		"""
		print(bepipred_file[-4:])
		assert bepipred_file[-4:] == ".tsv", "Currently support parsing only tsv file"

		# initialize epitopes csv file
		epitopes_filename = splitext(bepipred_file)[0] + "_epi" + ".tsv"
		if not exists(join(self.out_path, epitopes_filename)):
			print("Create file: {}".format(epitopes_filename))
			with open(join(self.out_path, epitopes_filename), "w") as epitopes_out:
				epitopes_out.write(
					"target_organism\ttarget_prot_id\ttarget_epi_start\ttarget_epi_end\ttarget_epi_sequence\tHLA restriction\tprediction_method\n")

		with open(join(self.bepipred_folder, bepipred_file)) as bepipred_out:
			current_epitope = {"start": 0, "end": 0, "seq": None, "protein_id": None}
			running_over_epi = False
			for bepipred_line in bepipred_out.readlines():
				bepipred_line = bepipred_line.strip()

				# determine the line type
				if bepipred_line[0:5] == "#####":  # header protein line
					continue
				elif bepipred_line[0:6] == "#Entry":  # header column line
					continue
				else:  # amino acid lines
					bepipred_info = BepiPredEpitopes.bepipred_tabs2info(bepipred_line)
					if running_over_epi:
						if bepipred_info["epi"] == "E":  # continue reading running epitope
							current_epitope["seq"] = current_epitope["seq"] + bepipred_info["AA"]
							current_epitope["end"] = bepipred_info["pos"]
						else:  # running over epitope is finished
							running_over_epi = False
							if len(current_epitope["seq"]) >= 9:
								self.write_epitope_row(current_epitope, epitopes_filename)
							else:
								print("Exclude current epitope: {} from output tsv as its length is lower that 9 amino-acids".format(current_epitope["seq"]))
							current_epitope = {"start": 0, "end": 0, "seq": None, "protein_id": None}
					else:
						if bepipred_info["epi"] == "E":  # the start of a new epitope is found
							running_over_epi = True
							current_epitope["protein_id"] = bepipred_info["entry"].split("|")[3]
							current_epitope["start"] = bepipred_info["pos"]
							current_epitope["seq"] = bepipred_info["AA"]
						else:  # prediction line not containing epitope between other not containing as well
							continue
