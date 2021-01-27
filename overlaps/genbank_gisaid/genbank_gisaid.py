import sys
from typing import Optional, List
from overlaps.multi_database_manager import config_db_engine, Session, Sequence, SequencingProject, get_session, \
    rollback, HostSample, target_sequences, Virus, Overlap, user_asked_to_commit, insert_overlaps_in_db, \
    cleanup_overlap_tables
from sqlalchemy import func, or_
from loguru import logger
from tqdm import tqdm
from os.path import sep
from datetime import date
import db_config.read_db_overlaps_configuration as db_config

# read only values
source_db_name = ""
source_database_source = ['GenBank', 'RefSeq']
source_name = 'GenBank'
target_db_name = ""
target_name = 'GISAID'


total_only_strain_1_to_n = 0
total_strain_plus_length_1_to_n = 0
total_only_strain_1_to_1 = 0
total_strain_plus_length_1_to_1 = 0
gisaid_only_false_tuples = 0
output_record = []


def source_sequences(session, database_source, virus_taxon_name, count_only: Optional[bool] = False):
    assert database_source is not None and virus_taxon_name is not None, \
        'If one between virus_taxon_name and database_source is None, the query extracting the source' \
        ' sequences is wrong'
    if count_only:
        query = session.query(func.count(Sequence.sequence_id))
    else:
        query = session.query(Sequence, HostSample)
    query = query \
            .select_from(SequencingProject, Virus, HostSample) \
            .filter(Sequence.strain_name.isnot(None),
                    SequencingProject.sequencing_project_id == Sequence.sequencing_project_id,
                    HostSample.host_sample_id == Sequence.host_sample_id,
                    Sequence.virus_id == Virus.virus_id,
                    Virus.taxon_name == virus_taxon_name,
                    Sequence.strain_name != 'NA',
                    func.length(Sequence.strain_name) > 2,
                    Sequence.sequence_id.notin_(
                        session.query(Overlap.sequence_id)
                        .filter(Overlap.overlapping_source == target_name)
                        .distinct())
                    )
    if len(database_source) == 1:
        query = query.filter(SequencingProject.database_source == database_source[0])
    elif len(database_source) == 2:
        or_clause = or_(SequencingProject.database_source == database_source[0],
                        SequencingProject.database_source == database_source[1])
        query = query.filter(or_clause)
    else:
        raise NotImplementedError()
    if count_only:
        return query.scalar()
    else:
        return query.yield_per(100)


def mark_overlaps():
    source_session = get_session(source_db_name)
    target_session = get_session(target_db_name)
    cleanup_overlap_tables(source_session, target_session)
    all_target_ids_changed = set()
    # dest_to_source_matches = {}
    global total_only_strain_1_to_n, total_strain_plus_length_1_to_n, output_record, total_only_strain_1_to_1, total_strain_plus_length_1_to_1, gisaid_only_false_tuples

    try:
        count_source_seq = source_sequences(session=source_session,
                                            database_source=source_database_source,
                                            virus_taxon_name='Severe acute respiratory syndrome coronavirus 2',
                                            count_only=True)
        for source_seq, source_host in tqdm(total=count_source_seq,
                                            iterable=source_sequences(session=source_session,
                                                                      database_source=source_database_source,
                                                                      virus_taxon_name='Severe acute respiratory syndrome coronavirus 2')):

            only_strain = []
            strain_plus_length = []

            target_seq_query = target_sequences(target_session,
                                                matching_strain=source_seq.strain_name)

            for target_seq in target_seq_query:
                if target_seq.length == source_seq.length:
                    strain_plus_length.append(target_seq)
                else:
                    only_strain.append(target_seq)

            if len(strain_plus_length) > 0:
                if len(strain_plus_length) > 1 and source_host.country is not None:
                    strain_plus_length = [a for a in strain_plus_length if (source_host.country.lower() in a.strain_name.strip().lower())]
                if len(strain_plus_length) > 1 and source_seq.strain_name.isnumeric():
                    candidates_from_switzerland = [g for g in strain_plus_length if "Switzerland" in g.strain_name]
                    candidates_from_england = [h for h in strain_plus_length if "England" in h.strain_name]
                    if len(candidates_from_switzerland) == 1:
                        strain_plus_length = candidates_from_switzerland
                        ids = [i.accession_id for i in strain_plus_length]
                        composed_string = f'{source_seq.accession_id} matches with {target_name} strain+length+only_numerical_strain+Switzerland_origin on {ids}.'
                    elif len(candidates_from_england) == len(strain_plus_length):
                        # if all targets are England sequences, it is a fake match
                        strain_plus_length = []
                        composed_string = ''
                    else:
                        strain_plus_length = filter_exact_strain_number_match(source_seq.strain_name, strain_plus_length)
                        ids = [i.accession_id for i in strain_plus_length]
                        composed_string = f'{source_seq.accession_id} matches with {target_name} strain+length on {ids}.'
                else:
                    strain_plus_length = filter_most_exact_strain_match(source_seq.strain_name, strain_plus_length)
                    ids = [i.accession_id for i in strain_plus_length]
                    composed_string = f'{source_seq.accession_id} matches with {target_name} strain+length (possibly also country) on {ids}.'
                if len(strain_plus_length) > 1:
                    total_strain_plus_length_1_to_n += 1
                    composed_string += f' WARNING:  1-N MATCH'
                elif len(strain_plus_length) == 1:
                    total_strain_plus_length_1_to_1 += 1
                if strain_plus_length:
                    output_record.append(composed_string)
                    output_record.append(f'\tstrain values: GENBANK: {source_seq.strain_name} GISAID: {[b.strain_name for b in strain_plus_length]}')

                    if len(set([s.sequence_id for s in strain_plus_length])) != len(strain_plus_length):
                        difference = len(strain_plus_length) - len(set([s.accession_id for s in strain_plus_length]))
                        logger.error(f'mismatch of {difference}')
                    for g in strain_plus_length:
                        gisaid_only_false_tuples += 1
                        g.gisaid_only = False
                        all_target_ids_changed.add(g.accession_id)
                        target_session.merge(g)

                        # UPDATE DEST TO SOURCE RELATION
                        # if not dest_to_source_matches.get(g.accession_id):
                        #     dest_to_source_matches[g.accession_id] = {
                        #         'GISAID_strain': g.strain_name,
                        #         'GenBank-ids': [source_seq.accession_id],
                        #         'GenBank-strains': [source_seq.strain_name]
                        #     }
                        # else:
                        #     dest_elem = dest_to_source_matches[g.accession_id]
                        #     dest_elem['GenBank-ids'].append(source_seq.accession_id)
                        #     dest_elem['GenBank-strains'].append(source_seq.strain_name)
                    insert_overlaps_in_db(source_session, target_session, source_seq, strain_plus_length, source_name,
                                          target_name)
            elif len(only_strain) > 0:
                if len(only_strain) > 1 and source_host.country is not None:
                    only_strain = [c for c in only_strain if (source_host.country.lower() in c.strain_name.strip().lower())]
                if len(only_strain) > 1 and source_seq.strain_name.isnumeric():
                    candidates_from_switzerland = [g for g in only_strain if "Switzerland" in g.strain_name]
                    candidates_from_england = [h for h in only_strain if "England" in h.strain_name]
                    if len(candidates_from_switzerland) == 1:
                        only_strain = candidates_from_switzerland
                        ids = [i.accession_id for i in only_strain]
                        composed_string = f'{source_seq.accession_id} matches with {target_name} strain_only+only_numerical_strain+Switzerland_origin on {ids}.'
                    elif len(candidates_from_england) == len(only_strain):
                        # if all targets are England sequences, it is a fake match
                        only_strain = []
                        composed_string = ''
                    else:
                        only_strain = filter_exact_strain_number_match(source_seq.strain_name, only_strain)
                        ids = [i.accession_id for i in only_strain]
                        composed_string = f'{source_seq.accession_id} matches with {target_name} strain_only on {ids}.'
                else:
                    only_strain = filter_most_exact_strain_match(source_seq.strain_name, only_strain)
                    ids = [i.accession_id for i in only_strain]
                    composed_string = f'{source_seq.accession_id} matches with {target_name} strain_only (possibly also country) on {ids}.'
                if len(only_strain) > 1:
                    total_only_strain_1_to_n += 1
                    composed_string += f' WARNING:  1-N MATCH'
                elif len(only_strain) == 1:
                    total_only_strain_1_to_1 += 1
                if only_strain:
                    output_record.append(composed_string)
                    output_record.append(
                        f'\tstrain values: GENBANK: {source_seq.strain_name} GISAID: {[y.strain_name for y in only_strain]}')

                    if len(set([s.sequence_id for s in only_strain])) != len(only_strain):
                        difference = len(only_strain) - len(set([s.accession_id for s in only_strain]))
                        logger.error(f'mismatch of {difference}')
                    for g in only_strain:
                        gisaid_only_false_tuples += 1
                        all_target_ids_changed.add(g.accession_id)
                        g.gisaid_only = False
                        target_session.merge(g)

                        # UPDATE DEST TO SOURCE RELATION
                        # if not dest_to_source_matches.get(g.accession_id):
                        #     dest_to_source_matches[g.accession_id] = {
                        #         'GISAID_strain': g.strain_name,
                        #         'GenBank-ids': [source_seq.accession_id],
                        #         'GenBank-strains': [source_seq.strain_name]
                        #     }
                        # else:
                        #     dest_elem = dest_to_source_matches[g.accession_id]
                        #     dest_elem['GenBank-ids'].append(source_seq.accession_id)
                        #     dest_elem['GenBank-strains'].append(source_seq.strain_name)
                    insert_overlaps_in_db(source_session, target_session, source_seq, only_strain, source_name,
                                          target_name)

            target_session.flush()
        logger.debug(f'total number of target sequences changed {len(all_target_ids_changed)}')
        if user_asked_to_commit:
            target_session.commit()
            source_session.commit()
    except KeyboardInterrupt:
        rollback(source_session)
        rollback(target_session)
        output_record.append("COMPUTATION INTERRUPTED. TOTALS MAY BE INCOMPLETE !!")
    except Exception as e:
        logger.exception("")
        rollback(source_session)
        rollback(target_session)
        output_record.append("COMPUTATION INTERRUPTED. TOTALS MAY BE INCOMPLETE !!")
    finally:
        source_session.close()
        target_session.close()

        totals_string = f'TOTALS:\n' \
                        f'1-1 MATCHES: strain+length: {total_strain_plus_length_1_to_1} -- only strain {total_only_strain_1_to_1}\n' \
                        f'1-N MATCHES: strain+length: {total_strain_plus_length_1_to_n} -- only strain {total_only_strain_1_to_n} (search "WARN" to find \'em)\n' \
                        f'{target_name} tuples in which to set SEQUENCE.gisaid_only to False {gisaid_only_false_tuples}\n'

        logger.info(totals_string)

        output_path = f'.{sep}overlaps{sep}{source_name}_{target_name}{sep}'
        output_path += f'{source_name}@{source_db_name}_overlapping_{target_name}@{target_db_name}'
        output_path += f'_{date.today().strftime("%Y-%b-%d")}.txt'
        output_path = output_path.lower()
        with open(file=output_path, mode='w') as w:
            for line in output_record:
                w.write(line + "\n")
            w.write(totals_string)
            w.write('\n\n')

            # WRITE DEST TO SOURCE RELATION
            # for key in dest_to_source_matches.keys():
            #     val = dest_to_source_matches[key]
            #     if len(val['GenBank-ids']) > 1:
            #         string = f"GISAID id {key} having strain {val['GISAID_strain']} matches GenBanks ids {val['GenBank-ids']} having strain {val['GenBank-strains']}"
            #         w.write(string+'\n')


def filter_exact_strain_number_match(source_strain_name: str, target_rows) -> []:
    result = []
    for row in target_rows:
        target_strain = row.strain_name
        start_index = target_strain.rfind(source_strain_name)
        preceding_char_index = start_index - 1
        following_char_index = start_index + len(source_strain_name)
        if (preceding_char_index < 0 or not target_strain[preceding_char_index].isnumeric()) and (
                following_char_index == len(target_strain) or not target_strain[following_char_index].isnumeric()):
            result.append(row)
    return result


def filter_most_exact_strain_match(source_strain_name: str, target_rows) -> []:
    result = []
    for row in target_rows:
        target_strain = row.strain_name
        start_index = target_strain.rfind(source_strain_name)
        following_char_index = start_index + len(source_strain_name)
        if following_char_index == len(target_strain) or not target_strain[following_char_index].isnumeric():
            result.append(row)
    return result


def run():
    genbank_db = db_config.get_import_params_for("genbank")
    gisaid_db = db_config.get_import_params_for("gisaid")
    global source_db_name, target_db_name
    source_db_name = genbank_db["db_name"]
    target_db_name = gisaid_db["db_name"]
    config_db_engine(genbank_db["db_name"], genbank_db["db_user"], genbank_db["db_psw"], genbank_db["db_port"])
    config_db_engine(gisaid_db["db_name"], gisaid_db["db_user"], gisaid_db["db_psw"], gisaid_db["db_port"])
    if not user_asked_to_commit:
        logger.warning('OPERATION WON\'T BE COMMITTED TO THE DB')
    mark_overlaps()


