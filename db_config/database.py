from typing import Optional

import sqlalchemy
from sqlalchemy import Column, select, join, Index, column, REAL
from sqlalchemy import String, Integer, Boolean, Float, Date, BigInteger
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.orm.session import Session
from sqlalchemy.engine import Engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.sql import text
from loguru import logger
from time import sleep

# https://www.compose.com/articles/using-postgresql-through-sqlalchemy/

_db_engine: Engine
_base = declarative_base()
_session_factory: sessionmaker
_last_config_parameters = ()


def config_db_engine(db_name, db_user, db_psw, db_port, recreate_db_from_scratch: bool = False):
    """
    Call this method once to initialize the database session factory and prepare it to execute queries.
    """
    global _db_engine, _session_factory, _last_config_parameters
    _last_config_parameters = (db_name, db_user, db_psw, db_port, recreate_db_from_scratch)
    logger.info('configuring db... make sure a connection is available')
    _db_engine = sqlalchemy.create_engine(f'postgresql://{db_user}:{db_psw}@localhost:{db_port}/{db_name}')\
        .execution_options(schema_translate_map={None: "public"})

    try:
        if recreate_db_from_scratch:
            logger.warning(
                'Removal of all table records in 10 seconds. Stop the execution if that\'s not the desired behaviour, and '
                'rerun by setting "recreate_db_from_scratch" to False in module main.py.')
            logger.warning(
                'Also, check if cache of annotations should be deleted!'
            )
            sleep(10)
            logger.info('removal of all table records in progress...')
            # delete_indexes()

            # DROP VIEWS (also removes dependencies on the tables)
            # for v in views:
            #     v.drop()

            # DROP TABLES
            _base.metadata.drop_all(_db_engine, tables=[ExperimentType.__table__, SequencingProject.__table__, Virus.__table__,
                                                        HostSample.__table__, Sequence.__table__, AminoAcidVariant.__table__,
                                                        AnnotationSequence.__table__, Annotation.__table__,
                                                        NucleotideVariant.__table__, NucleotideSequence.__table__,
                                                        VariantImpact.__table__, Epitope.__table__, EpitopeFragment.__table__,
                                                        DBMeta.__table__])

        # CREATE TABLES if not existing
        _base.metadata.create_all(_db_engine)   # throws sqlalchemy.exc.OperationalError if connection is not available

    except sqlalchemy.exc.OperationalError as e:
        logger.error('DB connection not available')
        raise e

    _session_factory = sessionmaker(bind=_db_engine, autocommit=False)
    logger.info('db configured')


def get_session() -> Session:
    return _session_factory()


def dispose_db_engine():
    _db_engine.dispose()


def re_config_db_engine(recreate_db_from_scratch=False):
    config_db_engine(*_last_config_parameters[0:4], recreate_db_from_scratch)


def rollback(session):
    try:
        session.rollback()
    except SQLAlchemyError:
        logger.exception('An error occurred during DB transaction. Rollback failed')


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
    # when the sessionmaker.autocommit flag is False, there is no need to begin a transaction with session.begin()
    try:
        result = func(session, *args, **kwargs)
        session.commit()
        return result
    except Rollback:
        logger.trace('Rollback of current transaction')
        rollback(session)
    except RollbackAndRaise as e:
        logger.trace('Rollback of current transaction')
        rollback(session)
        raise e.wrapped_exception
    except CommitAndRaise as e:
        session.commit()
        raise e.wrapped_exception
    except Exception as e:
        logger.info('Rollback of current transaction.')
        rollback(session)
        raise e
    finally:
        session.close()


#            ##################  EXCEPTIONS WITH SPECIAL BEHAVIOURS  #####################
class Rollback(Exception):
    pass


class RollbackAndRaise(Exception):
    def __init__(self, wrap_exception: Exception):
        self.wrapped_exception = wrap_exception


class CommitAndRaise(Exception):
    def __init__(self, wrap_exception: Exception):
        self.wrapped_exception = wrap_exception


#           ######################  VCM SCHEMA DEFINITION       #############################
class ExperimentType(_base):
    __tablename__ = 'experiment_type'

    experiment_type_id = Column(Integer, primary_key=True, autoincrement=True)

    sequencing_technology = Column(String)
    assembly_method = Column(String)
    coverage = Column(Integer)


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


class HostSpecie(_base):
    __tablename__ = 'host_specie'

    host_id = Column(Integer, primary_key=True, autoincrement=True)

    host_taxon_id = Column(Integer)
    host_taxon_name = Column(String)


class HostSample(_base):
    __tablename__ = 'host_sample'

    host_sample_id = Column(Integer, primary_key=True, autoincrement=True)
    host_id = Column(Integer)
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
    experiment_type_id = Column(Integer, nullable=False)
    virus_id = Column(Integer, nullable=False)
    host_sample_id = Column(Integer, nullable=False)
    sequencing_project_id = Column(Integer, nullable=False)

    accession_id = Column(String, nullable=False)
    alternative_accession_id = Column(String, nullable=True)
    strain_name = Column(String)
    is_reference = Column(Boolean, nullable=False)
    is_complete = Column(Boolean)
    strand = Column(String)
    length = Column(Integer)
    gc_percentage = Column(Float)
    n_percentage = Column(Float)
    lineage = Column(String)
    clade = Column(String)
    gisaid_only = Column(Boolean, default=True, nullable=False)


class NucleotideSequence(_base):
    __tablename__ = 'nucleotide_sequence'

    sequence_id = Column(Integer, primary_key=True)
    nucleotide_sequence = Column(String)


class Annotation(_base):
    __tablename__ = 'annotation'

    annotation_id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_id = Column(Integer, nullable=False)

    feature_type = Column(String, nullable=False)
    start = Column(Integer)
    stop = Column(Integer)
    gene_name = Column(String)
    product = Column(String)
    external_reference = Column(String)


class AnnotationSequence(_base):
    __tablename__ = 'annotation_sequence'

    annotation_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, nullable=False)

    product = Column(String)
    aminoacid_sequence = Column(String)
    annotation_nucleotide_sequence = Column(String)


class NucleotideVariant(_base):
    __tablename__ = 'nucleotide_variant'

    nucleotide_variant_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, nullable=False)

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
    nucleotide_variant_id = Column(Integer, nullable=False)

    effect = Column(String)
    putative_impact = Column(String)
    impact_gene_name = Column(String)


class AminoAcidVariant(_base):
    __tablename__ = 'aminoacid_variant'

    aminoacid_variant_id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, nullable=False)

    sequence_aa_original = Column(String, nullable=False)
    sequence_aa_alternative = Column(String, nullable=False)
    start_aa_original = Column(Integer)
    variant_aa_length = Column(Integer, nullable=False)
    variant_aa_type = Column(String, nullable=False)


class Epitope(_base):
    __tablename__ = 'epitope'

    epitope_id = Column(Integer, primary_key=True, autoincrement=True)
    epitope_iri = Column(String)
    iedb_epitope_id = Column(Integer)
    virus_id = Column(String, nullable=False)
    host_id = Column(Integer, nullable=False)
    source_host_name = Column(String)
    source_host_iri = Column(String)
    protein_ncbi_id = Column(String)
    cell_type = Column(String)
    mhc_class = Column(String)
    mhc_allele = Column(String)
    response_frequency_pos = Column(REAL)
    epitope_sequence = Column(String)
    epi_annotation_start = Column(Integer)
    epi_annotation_stop = Column(Integer)
    external_link = Column(String)
    prediction_process = Column(String)
    is_linear = Column(Boolean)
    assay_type = Column(String)


class EpitopeFragment(_base):
    __tablename__ = 'epitope_fragment'

    epi_fragment_id = Column(Integer, primary_key=True, autoincrement=True)
    epitope_id = Column(Integer)
    epi_fragment_sequence = Column(String)
    epi_frag_annotation_start = Column(Integer)
    epi_frag_annotation_stop = Column(Integer)


class DBMeta(_base):
    __tablename__ = 'db_meta'

    virus_id = Column(Integer, primary_key=True)
    source = Column(String, primary_key=True)
    date_of_import = Column(Date)


class Overlap(_base):
    __tablename__ = 'overlap'

    sequence_id = Column(Integer, nullable=False, primary_key=True)
    accession_id = Column(String, nullable=False)
    overlapping_accession_id = Column(String, nullable=False, primary_key=True)
    overlapping_source = Column(String, nullable=False)


class PipelineEvent(_base):
    __tablename__ = 'pipeline_event'

    event_id = Column(Integer, primary_key=True)
    event_name = Column(String)
    event_date = Column(String)
    added_items = Column(Integer)
    removed_items = Column(Integer)
    changed_items = Column(Integer)

#   ###################################     VIEWS       ##################################
# class View:
#     """
#     SQLAlchemy gives a create_view function but it doesn't checks
#     the existence of the view before creating, so without a prior DROP VIEW, the create_all raises Exception.
#     This class offers a workaround to create and drop views when necessary.
#     """
#
#     @staticmethod
#     def create():
#         raise NotImplementedError('Override this method to call _create_view with the correct parameters')
#
#     @staticmethod
#     def drop():
#         raise NotImplementedError('Override this method to call _drop_view with the correct parameters')
#
#     @staticmethod
#     def _create_view(view_name, view_stmt):
#         compiled_stmt = view_stmt.compile(compile_kwargs={"literal_binds": True}, dialect=_db_engine.dialect)
#         _db_engine.execute(f'CREATE OR REPLACE VIEW {view_name} AS {compiled_stmt}')
#
#     @staticmethod
#     def _drop_view(view_name):
#         _db_engine.execute(f'DROP VIEW IF EXISTS {view_name}')
#
#
# class MaterializedView:
#     """
#     SQLAlchemy gives a create_view function but it doesn't checks
#     the existence of the view before creating, so without a prior DROP VIEW, the create_all raises Exception.
#     This class offers a workaround to create and drop views when necessary.
#     """
#
#     @staticmethod
#     def create():
#         raise NotImplementedError('Override this method to call _create_view with the correct parameters')
#
#     @staticmethod
#     def drop():
#         raise NotImplementedError('Override this method to call _drop_view with the correct parameters')
#
#     @staticmethod
#     def _create_view(view_name, view_stmt):
#         MaterializedView._drop_view(view_name)
#         compiled_stmt = view_stmt.compile(compile_kwargs={"literal_binds": True}, dialect=_db_engine.dialect)
#         _db_engine.execute(f'CREATE MATERIALIZED VIEW {view_name} AS {compiled_stmt}')
#
#     @staticmethod
#     def _drop_view(view_name):
#         _db_engine.execute(f'DROP MATERIALIZED VIEW IF EXISTS {view_name}')
#
#
# class ViewAnnotationCDS(View):
#     stmt = select([
#         Annotation.annotation_id,
#         Annotation.sequence_id,
#         Annotation.start,
#         Annotation.stop,
#         Annotation.gene_name,
#         Annotation.product,
#         Annotation.external_reference,
#         Annotation.aminoacid_sequence
#     ]).where(Annotation.feature_type == 'CDS')
#
#     @staticmethod
#     def create():
#         ViewAnnotationCDS._create_view('annotation_cds', ViewAnnotationCDS.stmt)
#
#     @staticmethod
#     def drop():
#         ViewAnnotationCDS._drop_view('annotation_cds')
#
#
# class ViewAnnotation(View):
#     stmt = select([
#         Annotation.sequence_id,
#         Annotation.product.label('annotation_view_product'),
#         Annotation.aminoacid_sequence.label('annotation_view_aminoacid_sequence'),
#         Annotation.annotation_nucleotide_sequence.label('annotation_view_nucleotide_sequence')
#     ]).where(
#         # noqa              # == ignore warning on " != None" for this case
#         (Annotation.product != None)
#         & ((Annotation.aminoacid_sequence != None) | (Annotation.annotation_nucleotide_sequence != None))
#     )
#
#     @staticmethod
#     def create():
#         ViewAnnotation._create_view('annotation_view', ViewAnnotation.stmt)
#
#     @staticmethod
#     def drop():
#         ViewAnnotation._drop_view('annotation_view')
#
#
# class ViewNucleotideVariantAnnotation(MaterializedView):
#     stmt = select([
#         NucleotideVariant.nucleotide_variant_id,
#         Annotation.feature_type.label('n_feature_type'),
#         Annotation.gene_name.label('n_gene_name'),
#         Annotation.product.label('n_product')
#     ]).select_from(join(Annotation, NucleotideVariant,
#                         (NucleotideVariant.start_original >= Annotation.start) &
#                         (NucleotideVariant.start_original <= Annotation.stop) &
#                         (NucleotideVariant.sequence_id == Annotation.sequence_id)))
#
#     @staticmethod
#     def create():
#         ViewNucleotideVariantAnnotation._create_view('nucleotide_variant_annotation', ViewNucleotideVariantAnnotation.stmt)
#
#     @staticmethod
#     def drop():
#         ViewNucleotideVariantAnnotation._drop_view('nucleotide_variant_annotation')
#
#     # try:
#     #     __table__ = create_view('nucleotide_variant_annotation', stmt, _base.metadata)
#     # except sqlalchemy.exc.ProgrammingError:
#     #     pass    # view already exists
#
#
# class ViewNucleotideVariantLimited(View):
#     stmt = select([
#         NucleotideVariant
#     ]).where(NucleotideVariant.variant_length <= 20)
#
#     @staticmethod
#     def create():
#         ViewNucleotideVariantLimited._create_view('nucleotide_variant_limited', ViewNucleotideVariantLimited.stmt)
#
#     @staticmethod
#     def drop():
#         ViewNucleotideVariantLimited._drop_view('nucleotide_variant_limited')
#
#
# class HostSampleView(View):
#     stmt = select([
#         HostSample, HostSpecie
#     ]).select_from(join(HostSpecie, HostSpecie, HostSample.host_id == HostSpecie.host_id))
#
#     @staticmethod
#     def create():
#         HostSampleView._create_view('host_sample_view', HostSampleView.stmt)
#
#     @staticmethod
#     def drop():
#         HostSampleView._drop_view('host_sample_view')
#
#
# views = [ViewAnnotationCDS, ViewNucleotideVariantAnnotation, ViewNucleotideVariantLimited, ViewAnnotation, HostSampleView]
#
#
# def create_views():
#     # CREATE OR REPLACE VIEWS
#     for v in views:
#         v.create()


# #                   ##############################      INDEXES     ##################################
# noinspection SqlDialectInspection,SqlNoDataSourceInspection
# def delete_indexes():
#     # the following names must match the ones declared during the generation of the indexes (see code below)
#     indexes_to_drop = ['aa__ann_id', 'aa__var_type_lower', 'aa__start_original', 'aa__var_type_normal',
#                         'ann__seq_id', 'ann__start', 'ann__stop',
#                         'nuc_var__seq_id', 'nuc_var__start_alt', 'nuc_var__start_orig', 'nuc_var__length',
#                         'seq__experiment_id', 'seq__host_id', 'seq__seq_proj_id', 'seq__virus_id',
#                         'impact__var_id', 'nuc_var_ann__var_id'
#                         ]
#     for i in indexes_to_drop:
#         try:
#             _db_engine.execute(f'DROP INDEX {i}')
#         except sqlalchemy.exc.ProgrammingError:
#             pass
#
#
# # noinspection SqlNoDataSourceInspection,SqlDialectInspection
# def create_indexes():
#
#     def column_name(column_obj):
#         return str(column_obj).split('.', maxsplit=1)[1]
#
#     # one of the indexes depends on a materialized view. Check its existence before continuing
#     nuc_var_matview_exists = _db_engine.execute(
#         f"SELECT EXISTS ( SELECT FROM pg_catalog.pg_matviews WHERE matviewname = 'nucleotide_variant_annotation' )"
#         ).first().values()[0] is True
#     if not nuc_var_matview_exists:
#         raise RuntimeError('One of the indexes is based on a materialized view. First create the views, then indexes.\n'
#                            'Creation of indexes aborted. The database has not been changed.')
#
#     logger.info('Deleting previous version of indexes (if present)')
#     delete_indexes()
#     logger.info('Generating indexes...')
#
#     _db_engine.execute(f'CREATE INDEX aa__ann_id ON {AminoAcidVariant.__table__}({column_name(AminoAcidVariant.annotation_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX aa__var_type_lower ON {AminoAcidVariant.__table__}(lower({column_name(AminoAcidVariant.variant_aa_type)})) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX aa__start_original ON {AminoAcidVariant.__table__}({column_name(AminoAcidVariant.start_aa_original)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX aa__var_type_normal ON {AminoAcidVariant.__table__}({column_name(AminoAcidVariant.variant_aa_type)}) TABLESPACE default_ts;')
#
#     _db_engine.execute(f'CREATE INDEX ann__seq_id ON {Annotation.__table__}({column_name(Annotation.sequence_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX ann__start ON {Annotation.__table__}({column_name(Annotation.start)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX ann__stop ON {Annotation.__table__}({column_name(Annotation.stop)}) TABLESPACE default_ts;')
#
#     #            for now we'll keep the following index disabled
#     # _db_engine.execute(f'CREATE INDEX nuc_var__alt ON {NucleotideVariant.__table__}(lower({column_name(NucleotideVariant.sequence_alternative)}))')
#     _db_engine.execute(f'CREATE INDEX nuc_var__seq_id ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.sequence_id)}) TABLESPACE default_ts;')    # primary key
#     _db_engine.execute(f'CREATE INDEX nuc_var__start_alt ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.start_alternative)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX nuc_var__start_orig ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.start_original)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX nuc_var__length ON {NucleotideVariant.__table__}({column_name(NucleotideVariant.variant_length)}) TABLESPACE default_ts;')
#
#     _db_engine.execute(f'CREATE INDEX seq__experiment_id ON {Sequence.__table__}({column_name(Sequence.experiment_type_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX seq__host_id ON {Sequence.__table__}({column_name(Sequence.host_sample_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX seq__seq_proj_id ON {Sequence.__table__}({column_name(Sequence.sequencing_project_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE INDEX seq__virus_id ON {Sequence.__table__}({column_name(Sequence.virus_id)}) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE UNIQUE INDEX seq__accession_id ON {Sequence.__table__}(lower({column_name(Sequence.accession_id)})) TABLESPACE default_ts;')
#     _db_engine.execute(f'CREATE UNIQUE INDEX seq__alternative_accession_id ON {Sequence.__table__}(lower({column_name(Sequence.alternative_accession_id)})) TABLESPACE default_ts;')
#
#     _db_engine.execute(f'CREATE INDEX impact__var_id ON {VariantImpact.__table__}({column_name(VariantImpact.nucleotide_variant_id)}) TABLESPACE default_ts;')
#
#     _db_engine\
#         .execute(f'CREATE INDEX nuc_var_ann__var_id ON nucleotide_variant_annotation USING btree (nucleotide_variant_id) TABLESPACE default_ts;')


#                   ##############################      CHIMERA SEQUENCES       ######################
def disambiguate_chimera_sequences():
    stmt = """
        update sequence 
        set accession_id = CONCAT(accession_id, '_', (select right(virus.taxon_name,1) from virus where virus.virus_id = sequence.virus_id)), 
        alternative_accession_id = CONCAT(alternative_accession_id, '_', (select right(virus.taxon_name,1) from virus where virus.virus_id = sequence.virus_id)) 
        from virus 
        where accession_id in( 
            select distinct accession_id 
            from sequence 
            group by accession_id 
            having count(accession_id)>1
        )
    """
    _db_engine.execute(stmt)


def run_script(path: str, use_session: Optional[Session] = None):
    # Create an empty command string
    sql_command = ''

    if use_session is None:
        session = get_session()
    else:
        session = use_session

    try:
        with open(path, mode='r') as sql_file:
            for line in sql_file:
                # Ignore commented lines
                if not line.startswith('--') and line.strip('\n'):
                    # Append line to the command string
                    sql_command += " " + line.strip('\n')

                    # If the command string ends with ';', it is a full statement
                    if sql_command.endswith(';'):
                        # Try to execute statement and commit it
                        try:
                            result = session.execute(text(sql_command))
                            session.commit()
                            return result
                        # Assert in case of error
                        except Exception as e:
                            logger.exception(f'SQL command:\n{sql_command}\n Failed to execute. ')
                            rollback(session)
    finally:
        if use_session is None:
            session.close()



