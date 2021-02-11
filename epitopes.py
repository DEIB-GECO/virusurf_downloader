from typing import Optional
from os.path import sep
from loguru import logger
from vcm import vcm as vcm
from VirusGenoUtil.code.integrate_iedb_epitopes import epitopes_for_virus_taxon
from data_sources.ncbi_services import host_taxon_name_from_ncbi_taxon_id
from db_config import read_db_import_configuration as import_config, database
from time import sleep
import wget
from zipfile import ZipFile
import gzip, shutil
from tqdm import tqdm
from locations import remove_file

epitope_id_mappings = dict()


def import_epitopes(virus_taxon_id: int):
    # remember to update the data source
    logger.warning("Keep in mind to update the data source with 'python main.py download epitopes' before performing"
                   " the current action. This program will resume in 10 seconds.")
    try:
        sleep(10)
    except KeyboardInterrupt:
        return

    # begin
    db_params: dict = import_config.get_database_config_params()
    database.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])
    if virus_taxon_id in [2010960, 186539]:     # == bombali or reston ebolavirus
        logger.info(f'No epitopes available for virus {virus_taxon_id}.')
        return

    virus_db_id = virus_database_id(virus_taxon_id)

    if virus_db_id is None:
        raise Exception('Epitopes must be associated to a Virus DB entity. Before running epitopes, create the '
                        f'virus associated with taxon {virus_taxon_id}')
    if epitopes_already_imported(virus_db_id):
        logger.info('Epitopes for this virus are already imported into the DB. Aborting import.')
        return

    # run epitopes' code
    logger.debug(
        f'calling epitopes for virus taxon {virus_taxon_id} associated to DB virud_id {virus_db_id}')
    epitopes, fragments = epitopes_for_virus_taxon(virus_taxon_id)

    # write to file
    # epitopes_file = open(f'.{sep}epitopes.csv', mode='w')
    # epitopes_file.write('virus_db_id\thost_specie_db_id\thost_name\thost_iri\tprotein_ncbi_id\tcell_type\tmhc_class\t'
    #                     'mhc_allele\tresponse_frequency_positive\tassay_type\tseq\tstart\tstop\text_links\t'
    #                     'prediction_process\tis_linear\tepitope_iri\tiedb_epitope_id\n')
    # epitopes_fragm_file = open(f'.{sep}epitopes_fragments.csv', mode='w')
    # epitopes_fragm_file.write('damianos_epitope_id\tseq\tstart\tstop\n')

    def do(session: database.Session):
        global epitope_id_mappings
        try:
            for epitope in epitopes:
                # get contained values
                damianos_epitope_id, virus_taxon_id, host_iri, host_name, host_taxon_id, protein_ncbi_id, cell_type, \
                mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links, \
                prediction_process, is_linear, epitope_iri, iedb_epitope_id = epitope

                # put host specie foreign key
                host_specie_db_id = create_or_get_host_specie_db_id(session, host_taxon_id)

                # insert epitope in the DB
                epitope = (virus_db_id, host_specie_db_id, host_name, host_iri, protein_ncbi_id, cell_type,
                           mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links,
                           prediction_process, is_linear, epitope_iri, iedb_epitope_id)

                # write to file
                # types = (str(type(i)) for i in epitope)
                # items = (str(i) for i in epitope)
                # for i in zip(items, types):
                #     epitopes_file.write(f'{i[0], i[1]}\t')
                # epitopes_file.write('\n')

                epitope_db_id = vcm.create_epitope(session, epitope)
                # bind epitope ids from Damianos with the ones returned from database
                epitope_id_mappings[damianos_epitope_id] = epitope_db_id

            for fragment in fragments:
                _, damianos_epitope_id, seq, start, stop = fragment

                # bind epitope ids from Damianos with the ones returned from database
                try:
                    epitope_db_id = epitope_id_mappings[damianos_epitope_id]
                except KeyError as e:
                    raise KeyError(
                        f'the epitope fragment ID {damianos_epitope_id} does not appear in the epitope IDs. This epitope fragment'
                        f' will be not inserted into the DB.'
                    )

                fragment = (epitope_db_id, seq, start, stop)

                # write to file
                # types = (str(type(i)) for i in fragment)
                # items = (str(i) for i in fragment)
                # for i in zip(items, types):
                #     epitopes_fragm_file.write(f'{i[0], i[1]}\t')
                # epitopes_fragm_file.write('\n')

                vcm.create_epitope_fragment(session, fragment)

        except Exception as e:
            logger.exception('Exception occurred while computing and importing epitopes. Epitopes won\'t be inserted into the DB.')
            raise database.RollbackAndRaise(e)
        # finally:
            # epitopes_file.close()
            # epitopes_fragm_file.close()

    database.try_py_function(
        do
    )

    # insert one row for each linear epitope into epitope_fragment table
    database.run_script(f".{sep}sql_scripts{sep}insert_linear_epitopes_into_epi_fragments.sql")


def create_or_get_host_specie_db_id(session: database.Session, organism_ncbi_taxon_id):
    specie_db_id = vcm.get_specie_id(session, organism_ncbi_taxon_id)
    if not specie_db_id:
        organism_name_from_ncbi = host_taxon_name_from_ncbi_taxon_id(organism_ncbi_taxon_id)
        if organism_name_from_ncbi is not None:
            organism_name_from_ncbi = organism_name_from_ncbi.lower()
        specie_db_id = vcm.create_or_get_host_specie_alt(session, organism_name_from_ncbi, organism_ncbi_taxon_id)
    return specie_db_id


def virus_database_id(virus_taxon_id) -> Optional[int]:
    class Virus:
        @staticmethod
        def taxon_id():
            return virus_taxon_id

    def get_virus_db_id(session) -> Optional[int]:
        db_virus = vcm.get_virus(session, Virus())
        return db_virus.virus_id if db_virus else None

    return database.try_py_function(
        get_virus_db_id
    )


def epitopes_already_imported(virus_db_id):
    return database.try_py_function(
        vcm.check_existence_epitopes, virus_db_id)


def download_epitope_data() -> (str, str, str):
    """
    :return: the local file path of the updated files from IEDB. In order, the files are:
    - bcell_full_v3
    - tcell_full_v3
    - mhc_ligand_full_v3
    """
    tcell_url = "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
    bcell_url = "http://www.iedb.org/downloader.php?file_name=doc/bcell_full_v3.zip"
    mhc_ligand_url = "http://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip"
    download_dir = f".{sep}VirusGenoUtil{sep}data{sep}iedb_input{sep}cell_epitopes{sep}"

    tcell_file_name = "tcell_full_v3"
    bcell_file_name = "bcell_full_v3"
    mhc_ligand_file_name = "mhc_ligand_full"

    final_tcell_local_file_path = download_dir + tcell_file_name + ".csv.gz"
    final_bcell_local_file_path = download_dir + bcell_file_name + ".csv.gz"
    final_mhc_ligand_local_file_path = download_dir + mhc_ligand_file_name + ".csv.gz"

    download_tcell_local_file_path = download_dir + tcell_file_name + ".zip"
    download_bcell_local_file_path = download_dir + bcell_file_name + ".zip"
    download_mhc_ligand_local_file_path = download_dir + mhc_ligand_file_name + ".zip"

    def download_epitopes_data():
        # make sure the output path does not exist already, or wget assigns a trailing number to it
        remove_file(download_bcell_local_file_path)
        remove_file(download_tcell_local_file_path)
        remove_file(download_mhc_ligand_local_file_path)
        logger.info(f'downloading tcell_full from {tcell_url} ...')
        wget.download(tcell_url, download_tcell_local_file_path)
        logger.info(f'downloading bcell_full from {bcell_url} ...')
        wget.download(bcell_url, download_bcell_local_file_path)
        logger.info(f'downloading mhc_ligand_full from {mhc_ligand_url} ...')
        wget.download(mhc_ligand_url, download_mhc_ligand_local_file_path)
        logger.info('\n')

    def extract_epitopes_data():
        # unzip and gzip as requested by VirusGenoUtil library
        logger.info("transforming downloaded files as required by VirusGenoUtil...")
        io_list = zip(
            (download_tcell_local_file_path, download_bcell_local_file_path, download_mhc_ligand_local_file_path),
            (tcell_file_name, bcell_file_name, mhc_ligand_file_name),
            (final_tcell_local_file_path, final_bcell_local_file_path, final_mhc_ligand_local_file_path))
        for downloaded_file_path, file_name, output_file_path in tqdm(io_list):
            inner_file_name = file_name + ".csv"
            # unzip downloaded file into inner_file_name
            with ZipFile(file=downloaded_file_path, mode='r') as zipped_file:
                zipped_file.extract(member=inner_file_name, path=download_dir)
            # gzip extracted file
            with open(file=download_dir + inner_file_name, mode="rb") as inner_file:
                with gzip.open(output_file_path, mode='wb') as output_file:
                    shutil.copyfileobj(inner_file, output_file)
            # remove inner file as it is only an intermediate product
            remove_file(download_dir + inner_file_name)

    # download/update epitopes data from IEDB
    # IEDB does not send Content-Length HTTP header, so we cannot decide whether the source data is more recent than
    # the local data. Therefore we just download and replace the local files.
    download_epitopes_data()
    extract_epitopes_data()

    return final_bcell_local_file_path, final_tcell_local_file_path, final_mhc_ligand_local_file_path
