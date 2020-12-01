from typing import Optional
from overlaps.multi_database_manager import config_db_engine, Sequence, SequencingProject, get_session, rollback, \
    source_sequences, target_sequences, HostSample, user_asked_to_commit
from loguru import logger
from tqdm import tqdm
from vcm.vcm import create_or_get_host_sample, create_or_get_sequencing_project
from time import sleep
from os.path import sep

# read only values
source_db_name = 'vcm_11_1'
source_name = 'COGUK'
source_database_source = ['COG-UK']
target_db_name = 'vcm_gisaid_12'
target_name = 'GISAID'
target_database_source = None


total_only_strain_1_to_n = 0
total_strain_plus_length_1_to_n = 0
total_only_strain_1_to_1 = 0
total_strain_plus_length_1_to_1 = 0
total_only_strain_final = 0
total_strain_plus_length_final = 0
output_record = []

class MockedVirusSampleClass:
    # original values
    host_sample_db_object: Optional[HostSample] = None
    sequencing_project_db_object: Optional[SequencingProject] = None

    # values with which the base values should be overwritten
    host_sample_overwritten_values: Optional[dict] = None
    sequencing_project_overwritten_values: Optional[dict] = None

    def assign_base_values(self, host_sample_db_object, sequencing_project_db_object):
        self.host_sample_db_object = host_sample_db_object
        self.sequencing_project_db_object = sequencing_project_db_object

    def overwrite_host_sample_values(self, dictionary):
        self.host_sample_overwritten_values = dictionary

    def overwrite_sequencing_project_values(self, dictionary):
        self.sequencing_project_overwritten_values = dictionary

    # MOCKING VirusSample class
    # SEQUENCING PROJECT FIELDS
    def submission_date(self):
        return self.sequencing_project_overwritten_values.get("submission_date") or \
               self.sequencing_project_db_object.submission_date

    def sequencing_lab(self):
        return self.sequencing_project_overwritten_values.get("sequencing_lab") or \
               self.sequencing_project_db_object.sequencing_lab

    def bioproject_id(self):
        return self.sequencing_project_overwritten_values.get("bioproject_id") or \
               self.sequencing_project_db_object.bioproject_id

    def database_source(self):
        return self.sequencing_project_overwritten_values.get("database_source") or \
               self.sequencing_project_db_object.database_source

    # HOST SAMPLE FIELDS
    def originating_lab(self):
        return self.host_sample_overwritten_values.get("originating_lab") or \
               self.host_sample_db_object.originating_lab

    def collection_date(self):
        return self.host_sample_overwritten_values.get("collection_date") or \
               self.host_sample_db_object.collection_date

    def isolation_source(self):
        return self.host_sample_overwritten_values.get("isolation_source") or \
               self.host_sample_db_object.isolation_source

    def gender(self):
        return self.host_sample_overwritten_values.get("gender") or \
               self.host_sample_db_object.gender

    def age(self):
        return self.host_sample_overwritten_values.get("age") or \
               self.host_sample_db_object.age

    def country__region__geo_group(self):
        return self.host_sample_overwritten_values.get("country__region__geo_group") or \
               (self.host_sample_db_object.country, self.host_sample_db_object.region, self.host_sample_db_object.geo_group)


mocked_virus_sample_obj = MockedVirusSampleClass()


# noinspection PyTypeChecker
def put_gisaid_metadata_into_coguk(coguk_db_session, coguk_seq_id, gisaid_db_session, gisaid_seq_id):
    global mocked_virus_sample_obj
    # get GISAID's metadata
    gisaid_metadata = gisaid_db_session.query(HostSample, SequencingProject, Sequence)\
        .filter(Sequence.accession_id == gisaid_seq_id,
                SequencingProject.sequencing_project_id == Sequence.sequencing_project_id,
                HostSample.host_sample_id == Sequence.host_sample_id)\
        .with_entities(HostSample.isolation_source, HostSample.originating_lab, SequencingProject.sequencing_lab, SequencingProject.submission_date, Sequence.is_complete)\
        .one_or_none()
    if gisaid_metadata:
        isolation_source, originating_lab, sequencing_lab, submission_date, is_complete = gisaid_metadata

        # get COG-UK metadata for this sequence
        coguk_metadata = coguk_db_session.query(Sequence, HostSample, SequencingProject) \
            .filter(Sequence.accession_id == coguk_seq_id,
                    SequencingProject.sequencing_project_id == Sequence.sequencing_project_id,
                    HostSample.host_sample_id == Sequence.host_sample_id).one_or_none()
        sequence, host_sample, sequencing_project = coguk_metadata

        # update associated entities (less than one object per sequence)
        mocked_virus_sample_obj.assign_base_values(host_sample, sequencing_project)
        mocked_virus_sample_obj.overwrite_host_sample_values({
            "isolation_source": isolation_source,
            "originating_lab": originating_lab
        })
        mocked_virus_sample_obj.overwrite_sequencing_project_values({
            "submission_date": submission_date,
            "sequencing_lab": sequencing_lab
        })
        new_host_id = create_or_get_host_sample(coguk_db_session, mocked_virus_sample_obj, host_sample.host_id)
        new_seq_proj_id = create_or_get_sequencing_project(coguk_db_session, mocked_virus_sample_obj)
        # update sequence (unique object per row)
        sequence.is_complete = is_complete
        sequence.host_sample_id = new_host_id
        sequence.sequencing_project_id = new_seq_proj_id


def mark_overlaps():
    source_session = get_session(source_db_name)
    target_session = get_session(target_db_name)
    global total_only_strain_1_to_n, total_strain_plus_length_1_to_n, output_record, total_only_strain_final, \
        total_strain_plus_length_final, total_only_strain_1_to_1, total_strain_plus_length_1_to_1

    try:
        for source_seq in tqdm(source_sequences(source_session, database_source=source_database_source)):
        # for coguk_id in tqdm(read_coguk_overlapping_ids_from_file()):
        #     source_seq = source_session.query(Sequence).filter(Sequence.accession_id == coguk_id).one_or_none()
        #     if source_seq is None:
        #         continue

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
                ids = [i.accession_id for i in strain_plus_length]
                composed_string = f'{source_seq.accession_id} matches with {target_name} strain+length on {ids}.'
                if len(strain_plus_length) > 1:
                    composed_string += f' Only {ids[0]} is used for copying metadata'
                    total_strain_plus_length_1_to_n += 1
                else:
                    total_strain_plus_length_1_to_1 += 1
                output_record.append(composed_string)
                total_strain_plus_length_final += 1
                put_gisaid_metadata_into_coguk(source_session, source_seq.accession_id, target_session, strain_plus_length[0].accession_id)
                for x in strain_plus_length:
                    x.gisaid_only = False
                    target_session.merge(x)
            elif len(only_strain) > 0:
                ids = [i.accession_id for i in only_strain]
                composed_string = f'{source_seq.accession_id} matches with {target_name} strain only on {ids}.'
                if len(only_strain) > 1:
                    composed_string += f' Only {ids[0]} is used for copying metadata'
                    total_only_strain_1_to_n += 1
                else:
                    total_only_strain_1_to_1 += 1
                output_record.append(composed_string)
                total_only_strain_final += 1
                put_gisaid_metadata_into_coguk(source_session, source_seq.accession_id, target_session, only_strain[0].accession_id)
                for x in only_strain:
                    x.gisaid_only = False

        if user_asked_to_commit:
            source_session.commit()
            target_session.commit()
    except KeyboardInterrupt:
        logger.info("rollback of changes")
        rollback(source_session)
        rollback(target_session)
    except Exception:
        logger.exception("")
        logger.info("rollback of changes")
        rollback(source_session)
        rollback(target_session)
    finally:
        source_session.close()
        target_session.close()

        totals_string = f'TOTALS:\n' \
                        f'1-1 PAIRS strain+length: {total_strain_plus_length_1_to_1} -- only strain {total_only_strain_1_to_1}\n' \
                        f'1-N PAIRS strain+length: {total_strain_plus_length_1_to_n} -- only strain {total_only_strain_1_to_n}\n' \
                        f'PAIRS actually used for updating {source_name} metadata:' \
                        f'strain+length: {total_strain_plus_length_final} -- only strain {total_only_strain_final}\n'
        logger.info(totals_string)

        output_path = f'.{sep}overlaps{sep}{source_name}_{target_name}{sep}{source_name}_{target_name}_overlaps.txt'.lower()
        with open(file=output_path, mode='w') as w:
            for line in output_record:
                w.write(line + "\n")
            w.write(totals_string)


def read_coguk_overlapping_ids_from_file():
    logger.warning('PICKING KNOWN OVERLAPPING SEQUENCES BETWEEN vcm_11 AND vcm_gisaid_12 INSTEAD OF SCANNING THE DB\n'
                   'Stop if that\'s not the desired behavior')
    sleep(10)
    coguk_overlapping_vcm_11_ids = f'.{sep}overlaps{sep}coguk_gisaid{sep}coguk_overlapping_vcm_11_ids.txt'

    with open(coguk_overlapping_vcm_11_ids) as overlpping_coguk_ids:
        for line in overlpping_coguk_ids:
            yield line.rstrip()

def run(db_user, db_password, db_port):
    config_db_engine(source_db_name, db_user, db_password, db_port)
    config_db_engine(target_db_name, db_user, db_password, db_port)
    if not user_asked_to_commit:
        logger.warning('OPERATION WON\'T BE COMMITTED TO THE DB')
    mark_overlaps()


