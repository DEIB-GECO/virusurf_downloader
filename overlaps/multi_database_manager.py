from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.session import Session
from sqlalchemy.engine import Engine
from loguru import logger
from typing import Optional, List
from sqlalchemy import or_, func
from sqlalchemy import create_engine
from db_config.database_tom import Rollback, RollbackAndRaise, CommitAndRaise, rollback, \
    Sequence, SequencingProject, Virus

_db_engine: Engine
_base = declarative_base()

databases = {}
user_asked_to_commit = False


def config_db_engine(db_name, db_user, db_psw, db_port):
    """
    Call this method once to initialize the database session factory and prepare it to execute queries.
    """
    global databases
    if not databases.get(db_name):
        last_config_parameters = (db_name, db_user, db_psw, db_port)
        logger.info('configuring db... make sure a connection is available')
        db_engine = create_engine(f'postgresql://{db_user}:{db_psw}@localhost:{db_port}/{db_name}')

        session_factory = sessionmaker(db_engine)
        logger.info(f'db {db_name} configured and added to "databases" map')

        databases[db_name] = (db_engine, session_factory, last_config_parameters)
    return databases


def get_session(db_name: str) -> Session:
    try:
        return databases[db_name][1]()
    except KeyError:
        logger.error(f'database {db_name} is not configured. Call {__name__}.config_db_engine with {db_name} before requesting a session.')


def dispose_db_engine(db_name: str):
    databases[db_name][0].dispose()


def re_config_db_engine(db_name: str):
    config_db_engine(*databases[db_name][2])


def try_py_function(db_name: str, func, *args, **kwargs):
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
    session: Session =get_session(db_name)
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


def source_sequences(session, database_source: Optional[List[str]] = None, virus_taxon_name: Optional[str] = None):
    # select
    query = session.query(Sequence)
    # choose tables to join
    tables_in_from = [Sequence]
    if database_source:
        tables_in_from.append(SequencingProject)
    if virus_taxon_name:
        tables_in_from.append(Virus)
    query = query.select_from(*tables_in_from)
    # joins and filters
    if database_source:
        if len(database_source) > 2:
            raise NotImplementedError()
        query = query.filter(SequencingProject.sequencing_project_id == Sequence.sequencing_project_id)
        if len(database_source) == 1:
            query = query.filter(SequencingProject.database_source == database_source[0])
        else:
            query = query.filter(or_(SequencingProject.database_source == database_source[0],
                                     SequencingProject.database_source == database_source[1]))
    if virus_taxon_name:
        query = query.filter(Virus.virus_id == Sequence.virus_id,
                             Virus.taxon_name == virus_taxon_name)
    # ensure to take only sequence with a strain value otherwise we cannot compare anything with the target sequences
    query = query.filter(Sequence.strain_name.isnot(None))
    return query.yield_per(100)


def target_sequences(session, matching_strain:str, database_source: Optional[List[str]] = None, virus_taxon_name: Optional[str] = None) -> List:
    # select
    query = session.query(Sequence)
    # choose tables to join
    tables_in_from = [Sequence]
    if database_source:
        tables_in_from.append(SequencingProject)
    if virus_taxon_name:
        tables_in_from.append(Virus)
    query = query.select_from(*tables_in_from)
    # joins and filter with database_source
    if database_source:
        if len(database_source) > 2:
            raise NotImplementedError()
        # join with sequencing project
        query = query.filter(SequencingProject.sequencing_project_id == Sequence.sequencing_project_id)
        # filter
        if len(database_source) == 1:
            query = query.filter(SequencingProject.database_source == database_source[0])
        else:
            query = query.filter(or_(SequencingProject.database_source == database_source[0],
                                     SequencingProject.database_source == database_source[1]))
    # join and filter with virus taxon name
    if virus_taxon_name:
        query = query.filter(Virus.virus_id == Sequence.virus_id,
                             Virus.taxon_name == virus_taxon_name)
    # filter on strain name
    query = query.filter(Sequence.strain_name.isnot(None),
                         func.lower(Sequence.strain_name).contains(matching_strain.strip().lower()))
    return query.all()