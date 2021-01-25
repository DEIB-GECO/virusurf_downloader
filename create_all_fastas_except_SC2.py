from generate_fasta import generate_fasta
from data_sources.ncbi_any_virus.settings import known_settings
from db_config.read_db_import_configuration import set_db_name
from os import mkdir
from os.path import sep
from datetime import date
import sys
from time import sleep


print("This script expects <db_name> as argument")
try:
    db_name = sys.argv[1]
    print(f"Fasta files will be generated from DB {db_name} for each virus except sars_cov_2 and betacoronavirus_england_1")
    print("Program will resume in 10 s.")
    sleep(10)
    set_db_name(db_name)
except:
    print("Mandatory argument <db_name> is missing.")
    sys.exit(-1)


for _fasta_target, virus_import_parameters in known_settings.items():
    if not virus_import_parameters:
        raise ValueError(f'{_fasta_target} is not recognised as an importable virus')
    if _fasta_target == 'sars_cov_2' or _fasta_target == 'betacoronavirus_england_1':
        continue

    virus_txid = virus_import_parameters["virus_taxon_id"]
    virus_folder = virus_import_parameters["generated_dir_name"]
    fasta_name = f'{_fasta_target}_all_{date.today().strftime("%Y-%b-%d")}.fasta'
    fasta_path = generate_fasta(virus_txid, virus_folder, fasta_name)
    print(f'fasta generated for {_fasta_target} at path\n '
          f'{fasta_path}')
