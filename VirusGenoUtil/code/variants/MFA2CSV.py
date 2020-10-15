from Bio import AlignIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
from os.path import join, exists
from os import scandir
from pandas import read_csv
import urllib
import time
from VirusGenoUtil.code.utils import create_dir, is_fasta_file_extension, extract_orf_name


class MFA2CSV:
	"""
	Class to convert a multi-fasta alignment file (MFA)
	to variants csv per target sequence

	Example:
		pos:  012345678
		ref:  ATGGAGAGC
	target_1: ATGGT--GC

	then an example of produced variant csv for target1 will be:
	ref_id   |   target_id   |   ref_location  | ref | alt | type
	ncbi_ref | ncbi_target   |        4        |  A  |  T  | mismatch
	ncbi_ref | ncbi_target   |        5        |  G  |  -  | deletion
	"""

	def __init__(self, data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path):
		"""
		MFA2CSV constructor

		Parameters
		----------
		data_path : str
			absolute path where data lie
		alignments_folder : str
			alignments folder name
		xmls_folder : str
			xmls folder name
		ncbi_ref_id : str
			NCBI reference id
		out_path : str
			absolute path where output will be placed

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.out_path = out_path
		self.alignments_path = join(data_path, alignments_folder)
		self.xmls_path = join(data_path, xmls_folder)
		self.ref_id = ncbi_ref_id
		self.ref_record = None
		self.ref_seq, self.ref_aa = None, None
		self.entrez_email = open(join(data_path, "email.txt")).read().strip()
		Entrez.email = self.entrez_email
		print("Load email for Enterz access")
		self.gene_annots = {}
		create_dir(out_path)

	def parse_ref_xml(self, xml_file):
		"""
		Parse reference sequence from xml file (created by virulign)

		Parameters
		----------
		xml_file : str
			xml file name

		Returns
		-------
		None
		"""
		print("Parse ORF reference sequence from {}".format(xml_file))
		root = ET.parse(join(self.xmls_path, xml_file)).getroot()
		assert "referenceSequence" in root.attrib.keys(), "AssertionError: {} does not contain reference sequence".format(
			xml_file)
		self.ref_seq = Seq(root.attrib["referenceSequence"])
		self.ref_aa = None
		ref_name, ref_description = None, None
		if "name" in root.attrib.keys():
			ref_name = root.attrib["name"].strip()
			if "orf" in ref_name:
				ref_name = "ORF" + ref_name.split("orf")[1]
		if "description" in root.attrib.keys():
			ref_description = root.attrib["description"].strip()

		self.ref_record = SeqRecord(Seq(root.attrib["referenceSequence"]),
		                            id=self.ref_id, name=ref_name,
		                            description=ref_description)

		print("Parsed reference record:\n {}".format(self.ref_record))
		print("---")

	def fetch_genbank(self, seq_id):
		"""
		Fetch target sequence genbank from NCBI Entrez

		Parameters
		----------
		seq_id : str
			sequence NCBI id

		Returns
		-------
		file handle
			fetched genbank handle
		"""
		max_trials = 3
		while max_trials > 0:
			try:
				print("Fetch genbank")
				genbank_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="xml")
				retry = False
			except urllib.error.HTTPError as err:
				if err.code == 500 or err.code == 504:
					# server error or proxy error
					print("re-try")
					retry = True
				else:
					raise
			if retry:
				max_trials -= 1
				time.sleep(20)
			else:
				max_trials = 0
		return genbank_handle

	def parse_genbank(self, seq_id, genbank_handle):
		"""
		Parse fetched genbank handle

		Parameters
		----------
		seq_id : str
			sequence NCBI id
		genbank_handle: file handle
			fetched genbank handle

		Returns
		-------
		None
		"""
		print("Parse fetched genbank")
		self.gene_annots[seq_id] = {}
		for rec in Entrez.read(genbank_handle):
			print("for record")
			annotations = rec["GBSeq_feature-table"]
			for annotation in annotations:
				if annotation["GBFeature_key"] == "gene":
					gene_name = annotation["GBFeature_quals"][0]["GBQualifier_value"]
					location = annotation["GBFeature_location"].split("..")
					if gene_name not in self.gene_annots[seq_id]:
						self.gene_annots[seq_id][gene_name] = {"start": int(location[0]), "end": int(location[1])}

		print("Successfully loaded the following gene annotations: ")
		print(self.gene_annots[seq_id])

	def alignment2csv_nt(self, alignment):
		"""
		Convert MSA nucleic acid alignment to csv
		Credits: following loosy the VCF format tutorial at:
		https://faculty.washington.edu/browning/beagle/intro-to-vcf.html

		Parameters
		----------
		alignment : Bio.Align.MultipleSeqAlignment
			aligned target sequence

		Returns
		-------
		None
		"""
		print("Convert nucleic acid alignment of target id {} to variants csv".format(alignment.id))

		# parse the genbank to get the ORF start end positions
		if self.ref_record.id not in self.gene_annots:
			print("Fetch and then parse genbank annotations for target id: {}".format(alignment.id))
			fetched_genbank = self.fetch_genbank(self.ref_record.id)
			self.parse_genbank(self.ref_record.id, fetched_genbank)

		variants_name = alignment.id + "_variants.csv"
		if not exists(join(self.out_path, variants_name)):
			print("Create file: {}".format(variants_name))
			with open(join(self.out_path, variants_name), "w") as variants_out:
				variants_out.write("ORF,POS,REF_ID,TARGET_ID,REF,ALT\n")
		with open(join(self.out_path, variants_name), "a") as variants_out:
			print("Update file: {}".format(variants_name))
			orf_start_pos = self.gene_annots[self.ref_record.id][self.ref_record.name]["start"]
			for ref_pos, ref_base in enumerate(self.ref_record.seq):
				if ref_base != alignment.seq[ref_pos]:  # variant detected
					print("Write variant")
					variant_row = ",".join(
						[self.ref_record.name, str(orf_start_pos + ref_pos), self.ref_record.id, alignment.id, ref_base,
						 alignment.seq[ref_pos]])
					variants_out.write(variant_row + "\n")

	def alignment2csv_aa(self, alignment):
		"""
		Convert MSA amino acid alignment to csv
		Credits: following loosy the VCF format tutorial at:
		https://faculty.washington.edu/browning/beagle/intro-to-vcf.html

		Parameters
		----------
		alignment : Bio.Align.MultipleSeqAlignment
			aligned target sequence

		Returns
		-------
		None
		"""
		print("Convert amino acid alignment of target id {} to variants csv".format(alignment.id))
		# read the genbank to get the start end in order to sort the variants
		if self.ref_record.id not in self.gene_annots:
			print("Fetch and then parse genbank annotations for target id: {}".format(alignment.id))
			fetched_genbank = self.fetch_genbank(self.ref_record.id)
			self.parse_genbank(self.ref_record.id, fetched_genbank)

		variants_name = alignment.id + "_variants.csv"
		if not exists(join(self.out_path, variants_name)):
			print("Create file: {}".format(variants_name))
			with open(join(self.out_path, variants_name), "w") as variants_out:
				variants_out.write("ORF,ORF_POS,POS,REF_ID,TARGET_ID,REF,ALT\n")

		with open(join(self.out_path, variants_name), "a") as variants_out:
			print("Update file: {}".format(variants_name))
			orf_start_pos = self.gene_annots[self.ref_record.id][self.ref_record.name]["start"]
			orf_end_pos = self.gene_annots[self.ref_record.id][self.ref_record.name]["end"]
			for ref_pos, ref_base in enumerate(self.ref_aa):
				if ref_base != alignment.seq[ref_pos]:  # variant detected
					print("Write variant")
					variant_row = ",".join(
						[self.ref_record.name, str(orf_start_pos) + str(orf_end_pos), str(ref_pos + 1),
						 self.ref_record.id,
						 alignment.id, ref_base,
						 alignment.seq[ref_pos]])
					variants_out.write(variant_row + "\n")

	def sort_variants(self, is_amino_seq):
		"""
		Sort variants csv

		Parameters
		----------
		is_amino_seq : bool
			The MFA contains amino acid sequence (True),
			otherwise it contains nucleic acid sequence (False)

		Returns
		-------
		None
		"""
		print("Sort variants csvs")
		with scandir(self.out_path) as out_dir:
			for out_content in out_dir:
				if out_content.name.endswith(".csv") and out_content.is_file():
					variants_df = read_csv(join(self.out_path, out_content.name), sep=",", header=0)
					if is_amino_seq:
						variants_df.sort_values(by=["ORF_POS", "POS"], ascending=True, inplace=True)
						variants_df.drop("ORF_POS", axis=1, inplace=True)
					else:
						variants_df.sort_values(by="POS", ascending=True, inplace=True)
					variants_df.to_csv(join(self.out_path, out_content.name), sep=",", index=False)

	def run(self, xml_file, mfa_file):
		"""
		Convert MFA file to variants csv per target sequence

		Parameters
		----------
		xml_file : str
			XML file containing the sequence of the ORF element
		mfa_file : str
			MFA file containing the multiple-fasta alignment

		Returns
		-------
		None
		"""
		# parse ref
		self.parse_ref_xml(xml_file)
		# parse multi-fasta alignment
		assert is_fasta_file_extension(
			mfa_file), "AssertionError: Currently supporting only multi-fasta alignment files."
		print("Parse MAF")
		multiple_alignment = AlignIO.read(join(self.alignments_path, mfa_file), "fasta")
		for alignment in multiple_alignment:
			self.alignment2csv(alignment)
		print("---")

	def run_multiple_orfs(self, is_amino_seq):
		"""
		Create multiple MFA for ORFs in csv per target sequence
		using nucleotide-acid bases

		Parameters
		----------
		is_amino_seq : bool
			the MFA contains amino acid sequence (True),
			otherwise it contains nucleic acid sequence(False)

		Returns
		-------
		None
		"""
		with scandir(self.alignments_path) as alignments_dir:
			for alignments_content in alignments_dir:
				if is_fasta_file_extension(alignments_content.name) and alignments_content.is_file():
					# parse reference ORF
					self.parse_ref_xml(extract_orf_name(alignments_content.name) + ".xml")
					# convert MAF --> variants csv
					print("Parse MAF for {}".format(alignments_content.name))
					multiple_alignment = AlignIO.read(join(self.alignments_path, alignments_content.name), "fasta")
					if is_amino_seq:
						# get ORF amino acid sequence
						assert multiple_alignment[
							       0].id.upper() == self.ref_record.name.upper(), "AssertError: Code expects first ORF sequence to be the reference"
						self.ref_aa = multiple_alignment[0].seq
						for alignment in multiple_alignment[1::]:
							self.alignment2csv_aa(alignment)
					else:
						for alignment in multiple_alignment:
							self.alignment2csv_nt(alignment)
					print("---")
		self.sort_variants(is_amino_seq)
