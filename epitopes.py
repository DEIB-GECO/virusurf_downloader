import sys

from loguru import logger
from data_sources.ncbi_any_virus.ncbi_importer import prepared_parameters
import database_tom
import vcm
from VirusGenoUtil.code.integrate_iedb_epitopes import epitopes_for_virus_taxon


def import_epitopes(parameter_name: str):
    # find necessary arguments
    associated_parameters = prepared_parameters[parameter_name]
    bind_to_organism_taxon_id = associated_parameters[1]

    class Virus:
        @staticmethod
        def taxon_id():
            return bind_to_organism_taxon_id

    db_virus_id = database_tom.try_py_function(
        vcm.get_virus, Virus()
    )
    if db_virus_id is None:
        logger.error('Epitopes must be associated to a Virus DB entity. Before running epitopes, create the '
                     f'virus associated with taxon {bind_to_organism_taxon_id}')
        sys.exit(1)

    # run epitopes

    def create_and_add_to_db(session: database_tom.Session):
        logger.debug(f'calling epitopes for virus taxon {bind_to_organism_taxon_id} as bound to DB virus with id {db_virus_id}')
        try:
            epitopes, fragments = epitopes_for_virus_taxon(bind_to_organism_taxon_id)

            host_specie = session.query(database_tom.HostSpecie).filter(
                database_tom.HostSpecie.host_taxon_id == 9606
            ).first()
            host_specie_id = host_specie.host_id

            vcm.create_epitopes(session, epitopes, db_virus_id, host_specie_id)
            vcm.create_epitopes_fragments(session, fragments)
        except:
            logger.exception('Exception occurred while computing and importing epitopes. Epitopes won\'t be inserted into the DB.')
            raise database_tom.Rollback()

    if not database_tom.try_py_function(vcm.check_existence_epitopes, db_virus_id):
        database_tom.try_py_function(
            create_and_add_to_db
        )





