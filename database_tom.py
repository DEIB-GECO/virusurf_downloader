import sqlalchemy
from sqlalchemy import Column, ForeignKey
from sqlalchemy import String, Integer, Boolean, Float, Date
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.orm.session import Session
from sqlalchemy.engine import Engine
from sqlalchemy.exc import SQLAlchemyError
from loguru import logger

# https://www.compose.com/articles/using-postgresql-through-sqlalchemy/

_db_engine: Engine
_base = declarative_base()
_session_factory: sessionmaker


def config_db_engine(db_name, db_user, db_psw, db_port, recreate_db_from_scratch: bool = False):
    """
    Call this method once to initialize the database session factory and prepare it to execute queries.
    """
    global _db_engine, _session_factory
    logger.info('configuring db... make sure a connection is available')
    _db_engine = sqlalchemy.create_engine(f'postgresql://{db_user}:{db_psw}@localhost:{db_port}/{db_name}')

    try:
        if recreate_db_from_scratch:
            _base.metadata.drop_all(_db_engine, tables=[ExperimentType.__table__, SequencingProject.__table__, Virus.__table__,
                                                        HostSample.__table__, Sequence.__table__, Annotation.__table__,
                                                        Variant.__table__])

        # create tables if not existing
        _base.metadata.create_all(_db_engine)   # throws sqlalchemy.exc.OperationalError if connection is not available
    except sqlalchemy.exc.OperationalError as e:
        logger.error('DB connection not available')
        raise e
    _session_factory = sessionmaker(_db_engine)
    logger.info('db configured')


def try_py_function(func, *args, **kwargs):
    """
    Use this function to perform any action on the database.

    How it works: you first create a function that wraps the actions you want to do with the database. This function
    must have as first argument a session (sqlalchemy.orm.session.Session ) object + other arguments at your choice
    (both positional and key-value).
    Them you call try_py_function( your_function_name, your_arguments ).

    This function takes care of the session lifecycle, executing your function, committing or rolling-back if an exception
    occurs.

    :param func: your function name without parentheses
    :param args: your positional args
    :param kwargs: your kw args
    :return: the return value of your function.
    """
    session: Session = _session_factory()
    try:
        result = func(session, *args, **kwargs)
        session.commit()
        return result
    except SQLAlchemyError as e:
        try:
            session.rollback()
        except SQLAlchemyError:
            logger.exception('An error occurred during DB transaction. Rollback failed')
        raise e
    finally:
        session.close()



class ExperimentType(_base):
    __tablename__ = 'experiment_type'

    experiment_type_id = Column(Integer, primary_key=True, autoincrement=True)

    sequencing_technology = Column(String)
    assembly_method = Column(String)
    coverage = Column(String)


class SequencingProject(_base):
    __tablename__ = 'sequencing_project'

    sequencing_project_id = Column(Integer, primary_key=True, autoincrement=True)

    sequencing_lab = Column(String)
    # submission_date = Column(Date)
    submission_date = Column(String)
    database_source = Column(String)
    bioproject_id = Column(String)


class Virus(_base):
    __tablename__ = 'virus'
    #
    virus_id = Column(Integer, primary_key=True, autoincrement=True)

    taxon_id = Column(Integer)
    taxon_name = Column(String)

    family = Column(String)
    sub_family = Column(String)
    genus = Column(String)
    species = Column(String)
    equivalent_list = Column(String)
    molecule_type = Column(String)
    is_single_stranded = Column(Boolean)
    is_positive_stranded = Column(Boolean)

    def __str__(self):
        return f'ID:{self.virus_id} TAXON_ID:{self.taxon_id} TAXON_NAME:{self.taxon_name}'


class HostSample(_base):
    __tablename__ = 'host_sample'

    host_sample_id = Column(Integer, primary_key=True, autoincrement=True)

    host_taxon_id = Column(Integer)
    host_taxon_name = Column(String)

    collection_date = Column(String)
    isolation_source = Column(String)
    originating_lab = Column(String)
    country = Column(String)
    region = Column(String)
    geo_group = Column(String)

    #     extra
    age = Column(Integer)
    gender = Column(String)


class Sequence(_base):
    __tablename__ = 'sequence'

    sequence_id = Column(Integer, primary_key=True, autoincrement=True)
    # FKs
    experiment_type_id = Column(Integer, ForeignKey(ExperimentType.experiment_type_id), nullable=False)
    virus_id = Column(Integer, ForeignKey(Virus.virus_id), nullable=False)
    host_sample_id = Column(Integer, ForeignKey(HostSample.host_sample_id), nullable=False)
    sequencing_project_id = Column(Integer, ForeignKey(SequencingProject.sequencing_project_id), nullable=False)

    accession_id = Column(String, unique=True, nullable=False)
    alternative_accession_id = Column(String, unique=True, nullable=True)
    strain_name = Column(String)
    is_reference = Column(Boolean, nullable=False)
    is_complete = Column(Boolean)
    nucleotide_sequence = Column(String)
    strand = Column(String)
    length = Column(Integer)
    gc_percentage = Column(Float)
    linage = Column(String)
    clade = Column(String)


class Annotation(_base):
    __tablename__ = 'annotation'

    annotation_id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_id = Column(Integer, ForeignKey(Sequence.sequence_id), nullable=False)

    feature_type = Column(String, nullable=False)
    start = Column(Integer)
    stop = Column(Integer)
    gene_name = Column(String)
    product = Column(String)
    external_reference = Column(String)


class Variant(_base):
    __tablename__ = 'variant'

    variant_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, ForeignKey(Sequence.sequence_id), nullable=False)

    sequence_original = Column(String, nullable=False)
    sequence_alternative = Column(String, nullable=False)
    start_original = Column(Integer)
    start_alternative = Column(Integer)
    variant_length = Column(Integer, nullable=False)
    variant_type = Column(String, nullable=False)

    def get_list(self):
        return [self.start_original, self.variant_length, self.sequence_original, self.sequence_alternative,
                self.variant_type]

    def get_list_columns(self):
        return ['start', 'length', 'sequence_original', 'alt_sequence', 'variant_type']
