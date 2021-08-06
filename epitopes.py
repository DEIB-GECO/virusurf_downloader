import socket
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
from data_sources.ncbi_any_virus.settings import known_settings

epitope_id_mappings = dict()
mat_view_with_data_template_path = f'.{sep}sql_scripts{sep}epitope_views_n_indexes{sep}template{sep}epitope_4_virus+protein_mat_view_with_data_template.sql'
mat_view_with_no_data_template_path = f'.{sep}sql_scripts{sep}epitope_views_n_indexes{sep}template{sep}epitope_4_virus+protein_mat_view_with_no_data_template.sql'


def import_epitopes(virus_taxon_id: int):
    socket.setdefaulttimeout(120)  # sets a timeout on Biopython.Entrez requests (default limit is None = infinite)
                                   # importing ncbi_services involves the timeout to be set to 6
                                   # affects all socket created from now on

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
    logger.trace(f"socket timeout before epitopes is {socket.getdefaulttimeout()}")
    epitopes, fragments = epitopes_for_virus_taxon(virus_taxon_id)
    logger.trace(f"socket timeout after epitopes is {socket.getdefaulttimeout()}")
    map_protein_id_2_name = protein_id_2_name(virus_taxon_id)

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

                # add protein_name to epitope columns:
                prot_name = map_protein_id_2_name.get(protein_ncbi_id)

                # insert epitope in the DB
                epitope = (virus_db_id, host_specie_db_id, host_name, host_iri, prot_name, cell_type,
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
            vcm.DBCache.commit_changes()
        except Exception as e:
            logger.exception('Exception occurred while computing and importing epitopes. Epitopes won\'t be inserted into the DB.')
            vcm.DBCache.rollback_changes()
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


def protein_id_2_name(virus_taxon_id):
    this_virus_settings = [v for k, v in known_settings.items() if v['virus_taxon_id'] == virus_taxon_id]
    # it is necessary to unfold it because it is a list of dictionaries
    try:
        this_virus_settings = this_virus_settings[0]
    except IndexError:
        raise AssertionError(f"Couldn't find the annotation file path for virus with taxon_id"
                             f" {virus_taxon_id} (should be only numeric).")
    product_ncbi_id_2_name = dict()
    with open(this_virus_settings["annotation_file_path"], mode='r') as annotations_file:
        for line in annotations_file.readlines():
            _, _, _, _, _, product_name, product_id, sequence = line.rstrip().split('\t')
            if sequence and product_id != '.' and product_name != '.':
                product_ncbi_id_2_name[product_id] = product_name
    return product_ncbi_id_2_name


def protein_names_of_virus(virus_taxon_id):
    this_virus_settings = [v for k, v in known_settings.items() if v['virus_taxon_id'] == virus_taxon_id]
    # it is necessary to unfold it because it is a list of dictionaries
    try:
        this_virus_settings = this_virus_settings[0]
    except IndexError:
        raise AssertionError(f"Couldn't find the annotation file path for virus with taxon_id"
                             f" {virus_taxon_id} (should be only numeric).")
    protein_names = []
    with open(this_virus_settings["annotation_file_path"], mode='r') as annotations_file:
        for line in annotations_file.readlines():
            _, _, _type, coordinates, gene, product_name, code, sequence = line.rstrip().split('\t')
            if sequence != '.' and product_name != '.':
                protein_names.append(product_name)
    return protein_names


def generate_epitope_mat_view_n_indexes(virus_txid: int, ddl_with_data: bool, output_file_path: Optional[str] = None):
    max_db_object_length = 0                # NEEDED TO CHECK THAT MAXIMUM NAME LENGTH IS < 63
    all_proteins_materialized_views = []

    for protein in protein_names_of_virus(virus_txid):
        # generate protein short name
        short_protein_name = protein.replace('-', '_').replace('(', '').replace(')', '').replace(' ', '_')\
            .replace("'", '').replace('/', '_').replace('\\', '_').lower()
        short_protein_name = short_protein_name[:min(11, len(short_protein_name))]

        # generate materialized view and indexes for current virus and protein
        template_path = mat_view_with_data_template_path if ddl_with_data else mat_view_with_no_data_template_path
        with open(template_path, mode='r') as mat_view_template:
            mat_view_of_protein = f"-- MATERIALIZED VIEW AND INDEXES FOR VIRUS {virus_txid} AND PROTEIN {protein}\n"
            for line in mat_view_template.readlines():
                replaced_line = line\
                    .replace('$prot_name', protein.replace("'", "''"))\
                    .replace('$virus_id', str(virus_txid))\
                    .replace('$short_prot_name', short_protein_name)
                mat_view_of_protein += replaced_line

                # test db object name length (must be < 63 chars)
                object_to_test = None
                if replaced_line.startswith('CREATE INDEX '):
                    object_to_test = replaced_line.lstrip('CREATE INDEX ')
                elif replaced_line.startswith('CREATE MATERIALIZED VIEW '):
                    object_to_test = replaced_line.lstrip('CREATE MATERIALIZED VIEW ')
                if object_to_test:
                    object_to_test = object_to_test.lower().lstrip('public.')
                    if len(object_to_test) > 63:
                        print(replaced_line)
                        raise AssertionError(f'DATABASE OBJECT {object_to_test} ({len(object_to_test)} CHARS)'
                                             f' IN LINE {replaced_line} EXCEEDS '
                                             f'THE 63 CHARACTERS LENGTH LIMIT.')
                    max_db_object_length = max(max_db_object_length, len(object_to_test))

        # collect materialized views for all the proteins of this virus
        all_proteins_materialized_views.append(mat_view_of_protein)

    print(f"length of longest db_object registered for virus {virus_txid}: {max_db_object_length}")

    if output_file_path is not None:
        with open(output_file_path, mode='w') as out:
            for view_n_indexes in all_proteins_materialized_views:
                out.write(view_n_indexes)
                out.write('\n\n')
        print(f"epitope materialized views and indexes for virus  {virus_txid} created at path {output_file_path}")
    else:
        for view_n_indexes in all_proteins_materialized_views:
            print(view_n_indexes)
            print('\n\n')


def generate_epitope_mat_view_n_indexes_4_all_viruses(generate_into_directory_path: str, ddl_with_data: bool) -> None:
    for virus_short_name in known_settings.keys():
        virus_settings = known_settings[virus_short_name]
        virus_txid = virus_settings["virus_taxon_id"]
        if generate_into_directory_path:
            output_file_path = f'{generate_into_directory_path}{sep}{virus_short_name}_epitope_mat_views_n_indexes.sql'
        else:
            output_file_path = None
        generate_epitope_mat_view_n_indexes(virus_txid, ddl_with_data, output_file_path)


def generate_refresh_epitope_mat_view(virus_txid: int, output_file_path: Optional[str] = None):
    # find epitope mat view name template
    epitope_mat_view_template = ""
    with open(mat_view_with_data_template_path, mode='r') as mat_view_template:
        for line in mat_view_template.readlines():
            line = line.strip()
            if line.startswith("CREATE MATERIALIZED VIEW "):
                epitope_mat_view_template = "REFRESH MATERIALIZED VIEW " \
                                            + line.lstrip("CREATE MATERIALIZED VIEW public.") \
                                            + ";"
                break
    if not epitope_mat_view_template:
        raise AssertionError(f"Couldn't find the template name for the epitope materialized views in file {mat_view_with_data_template_path}")

    all_protein_refresh_commands = []
    for protein in protein_names_of_virus(virus_txid):
        # generate protein short name
        short_protein_name = protein.replace('-', '_').replace('(', '').replace(')', '').replace(' ', '_')\
            .replace("'", '').replace('/', '_').replace('\\', '_').lower()
        short_protein_name = short_protein_name[:min(11, len(short_protein_name))]

        all_protein_refresh_commands.append(
            epitope_mat_view_template
                .replace('$virus_id', str(virus_txid)) \
                .replace('$short_prot_name', short_protein_name)
        )

    if output_file_path:
        with open(output_file_path, mode="w") as out:
            for line in all_protein_refresh_commands:
                out.write(line + "\n")
    else:
        for line in all_protein_refresh_commands:
            print(line + "\n")


def generate_refresh_epitope_mat_view_4_all_viruses(generate_into_directory_path: str) -> None:
    for virus_short_name in known_settings.keys():
        virus_settings = known_settings[virus_short_name]
        virus_txid = virus_settings["virus_taxon_id"]
        if generate_into_directory_path:
            output_file_path = f'{generate_into_directory_path}{sep}{virus_short_name}_refresh_epitope_mat_views.sql'
        else:
            output_file_path = None
        generate_refresh_epitope_mat_view(virus_txid, output_file_path)
