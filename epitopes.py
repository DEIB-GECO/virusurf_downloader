from typing import Optional

from loguru import logger
from vcm import vcm as vcm
from VirusGenoUtil.code.integrate_iedb_epitopes import epitopes_for_virus_taxon
from data_sources.ncbi_services import host_taxon_name_from_ncbi_taxon_id
from db_config import read_db_import_configuration as import_config, database

epitope_id_mappings = dict()


def import_epitopes(virus_taxon_id: int):
    db_params:dict = import_config.get_database_config_params()
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
        f'calling epitopes for virus taxon {virus_taxon_id} as bound to DB virus with id {virus_db_id}')
    epitopes, fragments = epitopes_for_virus_taxon(virus_taxon_id)

    # write to file
    epitopes_file = open('./epitopes.csv', mode='w')
    epitopes_file.write('virus_db_id\thost_specie_db_id\thost_name\thost_iri\tprotein_ncbi_id\tcell_type\tmhc_class\tmhc_allele\tresponse_frequency_positive\tassay_type\tseq\tstart\tstop\text_links\tprediction_process\tis_linear\n')
    epitopes_fragm_file = open('./epitopes_fragments.csv', mode='w')
    epitopes_fragm_file.write('damianos_epitope_id\tseq\tstart\tstop\n')

    def do(session: database.Session):
        global epitope_id_mappings
        try:
            for epitope in epitopes:
                # get contained values
                damianos_epitope_id, virus_taxon_id, host_iri, host_name, host_taxon_id, protein_ncbi_id, cell_type, \
                mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links, \
                prediction_process, is_linear = epitope

                # put host specie foreign key
                host_specie_db_id = create_or_get_host_specie_db_id(session, host_taxon_id)

                # insert epitope in the DB
                epitope = (virus_db_id, host_specie_db_id, host_name, host_iri, protein_ncbi_id, cell_type,
                           mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links,
                           prediction_process, is_linear)

                # write to file
                types = (str(type(i)) for i in epitope)
                items = (str(i) for i in epitope)
                for i in zip(items, types):
                    epitopes_file.write(f'{i[0], i[1]}\t')
                epitopes_file.write('\n')

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
                types = (str(type(i)) for i in fragment)
                items = (str(i) for i in fragment)
                for i in zip(items, types):
                    epitopes_fragm_file.write(f'{i[0], i[1]}\t')
                epitopes_fragm_file.write('\n')

                vcm.create_epitope_fragment(session, fragment)

        except Exception as e:
            logger.exception('Exception occurred while computing and importing epitopes. Epitopes won\'t be inserted into the DB.')
            raise database.RollbackAndRaise(e)
        finally:
            epitopes_file.close()
            epitopes_fragm_file.close()

    database.try_py_function(
        do
    )


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
