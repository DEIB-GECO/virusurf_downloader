from os import makedirs


def create_dir(base_path):
	"""
	Create directory

	Parameters
	----------
	base_path : str
		full path of directory to be created

	Returns
	-------
	None
	"""
	makedirs(base_path, exist_ok=True)


def is_fasta_file_extension(file_name):
	"""
	Check if file has fasta extension

	Parameters
	----------
	file_name : str
		file name

	Returns
	-------
	bool
		True if file has fasta extension, False otherwise
	"""
	if file_name[-4:] == ".fna":
		return True
	elif file_name[-4:] == ".faa":
		return True
	elif file_name[-6:] == ".fasta":
		return True
	elif file_name[-6:] == ".fastq":
		return True
	elif file_name[-3:] == ".fa":
		return True
	else:
		return False


def extract_orf_name(file_name):
	"""
	Extract ORF name from MFA fasta file

	Parameters
	----------
	file_name : str
		file name
	Returns
	-------
	str
		extracted ORF name
	"""
	assert is_fasta_file_extension(file_name), "AssertionError: extract ORF name from the MFA file"
	if "_" in file_name:
		return file_name.split("_")[0]
