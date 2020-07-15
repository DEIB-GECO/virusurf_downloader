import sqlalchemy
from sqlalchemy import Column, ForeignKey, select, join, Index, column
from sqlalchemy import String, Integer, Boolean, Float, Date
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.orm.session import Session
from sqlalchemy.engine import Engine
from sqlalchemy.exc import SQLAlchemyError
from loguru import logger
from time import sleep
from sqlalchemy_utils import create_view

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
            logger.warning(
                'Removal of all table records in 10 seconds. Stop the execution if that\'s not the desired behaviour, and '
                'rerun by setting "recreate_db_from_scratch" to False in module main.py.')
            sleep(10)
            logger.info('removal of all table records in progress...')
            # DROP VIEWS (also removes dependencies on the tables)
            for v in views:
                v.drop()

            # DROP TABLES
            _base.metadata.drop_all(_db_engine, tables=[ExperimentType.__table__, SequencingProject.__table__, Virus.__table__,
                                                        HostSample.__table__, Sequence.__table__, AminoacidVariant.__table__,
                                                        Annotation.__table__, NucleotideVariant.__table__,
                                                        VariantImpact.__table__])

        # CREATE TABLES if not existing
        _base.metadata.create_all(_db_engine)   # throws sqlalchemy.exc.OperationalError if connection is not available
        # CREATE OR REPLACE VIEWS
        for v in views:
            v.create()
    except sqlalchemy.exc.OperationalError as e:
        logger.error('DB connection not available')
        raise e

    _session_factory = sessionmaker(_db_engine)
    logger.info('db configured')


def get_session() -> Session:
    return _session_factory()


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
        logger.trace('rollback of current transaction')
        rollback(session)
        raise e
    except RollbackTransactionAndRaise as e:
        logger.trace('rollback of current transaction')
        rollback(session)
        raise e
    except RollbackTransactionWithoutError as e:
        logger.trace('rollback of current transaction')
        rollback(session)
        if str(e) is not None:
            logger.error(str(e))
    except Exception as e:
        logger.info('Rollback of current transaction.')
        rollback(session)
        raise e
    finally:
        session.close()


def rollback(session):
    try:
        session.rollback()
    except SQLAlchemyError:
        logger.exception('An error occurred during DB transaction. Rollback failed')


class RollbackTransactionWithoutError(Exception):
    pass


class RollbackTransactionAndRaise(Exception):
    pass


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
    submission_date = Column(Date)
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
    n_percentage = Column(Float)
    lineage = Column(String)
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
    aminoacid_sequence = Column(String)
    annotation_nucleotide_sequence = Column(String)


class NucleotideVariant(_base):
    __tablename__ = 'nucleotide_variant'

    nucleotide_variant_id = Column(Integer, primary_key=True)
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


class VariantImpact(_base):
    __tablename__ = 'variant_impact'

    variant_impact_id = Column(Integer, primary_key=True)
    nucleotide_variant_id = Column(Integer, ForeignKey(NucleotideVariant.nucleotide_variant_id), nullable=False)

    effect = Column(String)
    putative_impact = Column(String)
    impact_gene_name = Column(String)


class AminoacidVariant(_base):
    __tablename__ = 'aminoacid_variant'

    aminoacid_variant_id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey(Annotation.annotation_id), nullable=False)

    sequence_aa_original = Column(String, nullable=False)
    sequence_aa_alternative = Column(String, nullable=False)
    start_aa_original = Column(Integer)
    variant_aa_length = Column(Integer, nullable=False)
    variant_aa_type = Column(String, nullable=False)


#   ###################################     VIEWS       ##################################
class View:
    """
    SQLAlchemy gives a create_view function but it doesn't checks
    the existence of the view before creating, so without a prior DROP VIEW, the create_all raises Exception.
    This class offers a workaround to create and drop views when necessary.
    """

    @staticmethod
    def create():
        raise NotImplementedError('Override this method to call _create_view with the correct parameters')

    @staticmethod
    def drop():
        raise NotImplementedError('Override this method to call _drop_view with the correct parameters')

    @staticmethod
    def _create_view(view_name, view_stmt):
        compiled_stmt = view_stmt.compile(compile_kwargs={"literal_binds": True}, dialect=_db_engine.dialect)
        _db_engine.execute(f'CREATE OR REPLACE VIEW {view_name} AS {compiled_stmt}')

    @staticmethod
    def _drop_view(view_name):
        _db_engine.execute(f'DROP VIEW IF EXISTS {view_name}')


class ViewAnnotationCDS(View):
    stmt = select([
        Annotation.annotation_id,
        Annotation.sequence_id,
        Annotation.start,
        Annotation.stop,
        Annotation.gene_name,
        Annotation.product,
        Annotation.external_reference,
        Annotation.aminoacid_sequence
    ]).where(Annotation.feature_type == 'CDS')

    @staticmethod
    def create():
        View._create_view('annotation_cds', ViewAnnotationCDS.stmt)

    @staticmethod
    def drop():
        View._drop_view('annotation_cds')


class ViewNucleotideVariantAnnoatation(View):
    stmt = select([
        NucleotideVariant.nucleotide_variant_id,
        Annotation.feature_type.label('n_feature_type'),
        Annotation.gene_name.label('n_gene_name'),
        Annotation.product.label('n_product')
    ]).select_from(join(Annotation, NucleotideVariant,
                        (NucleotideVariant.start_alternative >= Annotation.start) &
                        (NucleotideVariant.start_alternative <= Annotation.stop) &
                        (NucleotideVariant.sequence_id == Annotation.sequence_id)))

    @staticmethod
    def create():
        View._create_view('nucleotide_variant_annotation', ViewNucleotideVariantAnnoatation.stmt)

    @staticmethod
    def drop():
        View._drop_view('nucleotide_variant_annotation')

    # try:
    #     __table__ = create_view('nucleotide_variant_annotation', stmt, _base.metadata)
    # except sqlalchemy.exc.ProgrammingError:
    #     pass    # view already exists


class ViewNucleotideVariantLimited(View):
    stmt = select([
        NucleotideVariant
    ]).where(NucleotideVariant.variant_length <= 20)

    @staticmethod
    def create():
        View._create_view('nucleotide_variant_limited', ViewNucleotideVariantLimited.stmt)

    @staticmethod
    def drop():
        View._drop_view('nucleotide_variant_limited')


views = [ViewAnnotationCDS, ViewNucleotideVariantAnnoatation, ViewNucleotideVariantLimited]


def create_indexes():
    logger.warning('Generation of indexes: This operation will blindly add new indexes without checking prior existence'
                   'Stop the execution now if that\'s not the desired behaviour')
    sleep(10)
    logger.info('Generating indexes...')

    def column_name(column_obj):
        return str(column_obj).split('.', maxsplit=1)[1]

    _db_engine.execute(f'CREATE INDEX ON {AminoacidVariant.__table__}({column_name(AminoacidVariant.annotation_id)})')
    _db_engine.execute(f'CREATE INDEX ON {AminoacidVariant.__table__}(lower({column_name(AminoacidVariant.variant_aa_type)}))')
    _db_engine.execute(f'CREATE INDEX ON {AminoacidVariant.__table__}({column_name(AminoacidVariant.start_aa_original)})')
    _db_engine.execute(f'CREATE INDEX ON {AminoacidVariant.__table__}({column_name(AminoacidVariant.variant_aa_type)})')

    _db_engine.execute(f'CREATE INDEX ON {Annotation.__table__}({column_name(Annotation.sequence_id)})')
    _db_engine.execute(f'CREATE INDEX ON {Annotation.__table__}({column_name(Annotation.start)})')
    _db_engine.execute(f'CREATE INDEX ON {Annotation.__table__}({column_name(Annotation.stop)})')

    #            for now we'll keep the following index disabled
    # _db_engine.execute(f'CREATE INDEX ON {NucleotideVariant.__table__}(lower({column_name(NucleotideVariant.sequence_alternative)}))')
    _db_engine.execute(f'CREATE INDEX ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.sequence_id)})')    # primary key
    _db_engine.execute(f'CREATE INDEX ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.start_alternative)})')
    _db_engine.execute(f'CREATE INDEX ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.start_original)})')
    _db_engine.execute(f'CREATE INDEX ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.variant_length)})')

    _db_engine.execute(f'CREATE INDEX ON {Sequence.__table__}({column_name(Sequence.experiment_type_id)})')
    _db_engine.execute(f'CREATE INDEX ON {Sequence.__table__}({column_name(Sequence.host_sample_id)})')
    _db_engine.execute(f'CREATE INDEX ON {Sequence.__table__}({column_name(Sequence.sequencing_project_id)})')
    _db_engine.execute(f'CREATE INDEX ON {Sequence.__table__}({column_name(Sequence.virus_id)})')

    _db_engine.execute(f'CREATE INDEX ON {VariantImpact.__table__}({column_name(VariantImpact.nucleotide_variant_id)})')