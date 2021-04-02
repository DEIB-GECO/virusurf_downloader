import json
import os
from os.path import exists, sep
from locations import get_local_folder_for, FileType, remove_file
from datetime import date
from json import JSONDecodeError
from typing import Generator, List, Tuple
from loguru import logger
from tqdm import tqdm
import bz2
from data_sources.gisaid_sars_cov_2.sample import GISAIDSarsCov2Sample, round_
from data_sources.virus import VirusSource
from db_config.database import try_py_function, Sequence, HostSample, SequencingProject, Annotation, AminoAcidVariant
from collections import Counter


# noinspection PyMethodMayBeStatic


class GISAIDSarsCov2(VirusSource):

    name = 'gisaid_sars_cov_2'
    data_path = get_local_folder_for(name, FileType.SequenceOrSampleData) + 'export.json'
    credentials_path = f'.{sep}data_sources{sep}gisaid_sars_cov_2{sep}credentials_gisaid.csv'

    def __init__(self):
        super().__init__()
        logger.info(f'importing virus {GISAIDSarsCov2.name}')
        self.update_source_data()

    def taxon_id(self):
        return 2697049

    def taxon_name(self):
        return 'Severe acute respiratory syndrome coronavirus 2'

    def family(self):
        return 'Coronaviridae'

    def sub_family(self):
        return 'Orthocoronavirinae'

    def genus(self):
        return 'Betacoronavirus'

    def species(self):
        return 'Severe acute respiratory syndrome-related coronavirus'

    def equivalent_names(self):
        return 'SARS-CoV-2, 2019-nCoV, COVID-19, COVID-19 virus, COVID19, HCoV-19, Human coronavirus 2019, SARS-2, SARS-CoV2, SARS2, Wuhan coronavirus, Wuhan seafood market pneumonia virus'

    def molecule_type(self):
        return 'RNA'

    def is_single_stranded(self):
        return True

    def is_positive_stranded(self):
        return True

    def count_sequences_in_file(self):
        with open(self.data_path, mode='r') as input_file:
            num_lines = sum(1 for line in input_file)
        return num_lines

    def update_source_data(self):
        logger.info("Downloading updates from source...")
        # read user credentials
        if not os.path.exists(self.credentials_path):
            with open(self.credentials_path, mode='w') as credentials_file:
                credentials_file.write("# Lines starting with # are comments.\n"
                                       "# Write in the following line <username>,<password> to use for downloading "
                                       "updated sequence data from GISAID.")
            raise AssertionError(f"No GISAID credentials provided. Please update the file at path {self.credentials_path}")
        with open(self.credentials_path, mode='r') as credentials_file:
            for line in credentials_file.readlines():
                if line.startswith("#"):
                    continue
                try:
                    username, psw = line.split(",")
                    username = username.strip().rstrip()
                    psw = psw.strip().rstrip()
                except Exception as e:
                    logger.error(f"Error encountered while parsing GISAID credentials file at path {self.credentials_path}")
                    raise e
            if not username or not psw:
                raise AssertionError(f"No GISAID credentials provided. Please update the file at path {self.credentials_path}")
        # download updated data from source
        download_path = get_local_folder_for(self.name, FileType.SequenceOrSampleData)
        download_path += "export_" + date.today().strftime("%Y-%b-%d") + ".json.bz2"
        remove_file(download_path)
        remove_file(self.data_path)
        os.system(f"wget --user {username} --password {psw} -O {download_path} ***REMOVED***")
        if not exists(download_path):
            raise ValueError("download of ***REMOVED*** with username "
                             f"'{username}' and password '{psw}' failed.")
        # extract archive to self.data_path
        with bz2.open(filename=download_path, mode='rt') as compressed_file:
            with open(file=self.data_path, mode="w") as decompressed_file:
                for line in compressed_file:
                    decompressed_file.write(line)

    def get_sequences_in_current_data(self) -> dict:
        def do(session):
            all_sequences = dict()
            for db_items in session.query(Sequence, HostSample, SequencingProject) \
                                      .filter(Sequence.strain_name.isnot(None),
                                              Sequence.host_sample_id == HostSample.host_sample_id,
                                              Sequence.sequencing_project_id == SequencingProject.sequencing_project_id) \
                                      .yield_per(100):
                source_seq, source_host, source_prj = db_items
                all_sequences[source_seq.accession_id] = (source_seq.strain_name, int(source_seq.length),
                                                          source_seq.gc_percentage, source_seq.n_percentage,
                                                          source_host.collection_date, source_host.originating_lab,
                                                          str(source_prj.submission_date) if source_prj.submission_date is not None else None,
                                                          source_prj.sequencing_lab,
                                                          source_host.province, source_host.region, source_host.country,
                                                          source_host.isolation_source)
            return all_sequences

        return try_py_function(do)

    def get_sequences_of_updated_source(self) -> Generator[GISAIDSarsCov2Sample, None, None]:
        with open(self.data_path, mode='r') as input_file:
            for line in input_file:
                try:
                    yield GISAIDSarsCov2Sample(json.loads(line))
                except JSONDecodeError:
                    pass

    # ### COMPARE ANNOTATIONS AND AMINO ACID VARIANTS ### #
    # It's far more efficient to download the annotations of a group of sequences altogether instead of
    # downloading them for each sequence at a time
    def collect_aa_variants_from_db(self, sample_acc_ids: list):
        all_db_annotations = dict()     # acc_id -> [(aa_var1), (aa_var_2), ...]

        def do(session):
            nonlocal all_db_annotations
            all_aa_variants = session \
                .query(Sequence.accession_id, Annotation.start, Annotation.stop, Annotation.gene_name,
                       Annotation.feature_type, Annotation.product,
                       AminoAcidVariant.start_aa_original, AminoAcidVariant.sequence_aa_original,
                       AminoAcidVariant.sequence_aa_alternative) \
                .select_from(Annotation, AminoAcidVariant, Sequence) \
                .filter(Annotation.annotation_id == AminoAcidVariant.annotation_id,
                        Annotation.sequence_id == Sequence.sequence_id,
                        Sequence.accession_id.in_(sample_acc_ids)) \
                .order_by(Sequence.accession_id, Annotation.start, AminoAcidVariant.start_aa_original) \
                .all()
            for acc_id, ann_start, ann_stop, gene_name, feature_type, product, \
                start_aa_original, sequence_aa_original,sequence_aa_alternative in all_aa_variants:
                # group aa_variants by accession id into a map
                try:
                    seq_annotations = all_db_annotations[acc_id]
                except KeyError:
                    seq_annotations = list()
                    all_db_annotations[acc_id] = seq_annotations
                content = (ann_start, ann_stop, gene_name, feature_type, product,
                           start_aa_original, sequence_aa_original, sequence_aa_alternative)
                seq_annotations.append(content)

        try_py_function(do)
        return all_db_annotations

    def annotations_changed(self, local_annotations_list: List[Tuple], sample: GISAIDSarsCov2Sample):
        local_ann = Counter(local_annotations_list)
        updated_ann = Counter(self.collect_aa_variants_from_file(sample))
        annotations_in_db_not_in_file = local_ann - updated_ann
        annotations_in_file_not_in_db = updated_ann - local_ann

        if sum(annotations_in_db_not_in_file.values()) > 0 or sum(annotations_in_file_not_in_db.values()) > 0:
            return True
        else:
            return False

    def collect_aa_variants_from_file(self, sample: GISAIDSarsCov2Sample):
        formatted_annotations = []
        for start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, _mutations in sample.annotations_and_amino_acid_variants():
            for mut in _mutations:
                original_aa, alternative_aa, start_pos, length, _type = mut
                content = (start, stop, gene_name, feature_type, product, start_pos, original_aa, alternative_aa)
                formatted_annotations.append(content)
        return formatted_annotations

    def deltas(self):
        """
        Returns the set of Sequence.accession_id to remove from current database, and the set of Sequence.accession_id
        that should be imported from the source. Note that because sequences in the current data may have been updated,
        the two sets are not necessarily disjoint. So first remove the sequences to be removed and only then
         insert the new ones.
         """
        acc_id_remote = set([x.primary_accession_number() for x in self.get_sequences_of_updated_source()])  # all the sequences from remote
        logger.info('Collecting current metadata...')
        current_data = self.get_sequences_in_current_data()
        acc_id_current = set(current_data.keys())

        acc_id_missing_in_remote = acc_id_current - acc_id_remote
        acc_id_missing_in_current = acc_id_remote - acc_id_current

        acc_id_present_in_current_and_remote = acc_id_current & acc_id_remote
        acc_id_changed_updatable = dict()
        acc_id_changed_not_updatable = []

        changes_distribution = {
            "sequence": 0,
            "host_sample": 0,
            "sequencing_project": 0,
            "annotations": 0
        }

        # collect annotations from DB only for sequences that can have changes
        logger.info('Collecting current annotations...')
        aa_variants_local = self.collect_aa_variants_from_db(list(acc_id_present_in_current_and_remote))
        # to compare sequences present in both, have to scan the file of remote sequences
        logger.info('Comparing current data with updated source data')
        for new_sequence in tqdm(self.get_sequences_of_updated_source(), total=self.count_sequences_in_file()):
            acc_id = new_sequence.primary_accession_number()
            try:
                current_sequence_data = current_data[acc_id]
                # if the gisaid id is the same, compare metadata to detect if the sample changed over time

                # if strain or length changes, overlaps of this sequene become invalid, otherwise they can be kept
                if current_sequence_data[0] != new_sequence.strain() \
                        or current_sequence_data[1] != new_sequence.length():
                    acc_id_changed_not_updatable.append(acc_id)
                else:
                    # detect which tables have changes
                    changes = {
                        "sequence": False,
                        "host_sample": False,
                        "sequencing_project": False,
                        "annotations": False
                    }
                    # changes in sequence table (implies a likely change in the sequence and variants too)
                    if current_sequence_data[2] != new_sequence.gc_percent() \
                            or current_sequence_data[3] != new_sequence.n_percent():
                        changes["sequence"] = True
                        changes_distribution["sequence"] = changes_distribution["sequence"] + 1
                    if self.annotations_changed(aa_variants_local[acc_id], new_sequence):
                        changes["annotations"] = True
                        changes_distribution["annotations"] = changes_distribution["annotations"] + 1
                    # changes in host sample table
                    if current_sequence_data[4] != new_sequence.collection_date()[0] \
                            or current_sequence_data[5] != new_sequence.originating_lab() \
                            or (current_sequence_data[8], current_sequence_data[9], current_sequence_data[10]) != \
                            new_sequence.province__region__country__geo_group()[0:3] \
                            or current_sequence_data[11] != new_sequence.isolation_source():
                        changes["host_sample"] = True
                        changes_distribution["host_sample"] = changes_distribution["host_sample"] + 1
                    # changes in sequencing project table
                    if current_sequence_data[6] != new_sequence.submission_date() \
                            or current_sequence_data[7] != new_sequence.sequencing_lab():
                        changes["sequencing_project"] = True
                        changes_distribution["sequencing_project"] = changes_distribution["sequencing_project"] + 1

                    if True in changes.values():
                        acc_id_changed_updatable[acc_id] = changes
            except KeyError:
                pass  # the accession id is not present in current data. it's a new sequence

        # compute additional sets
        acc_id_changed_not_updatable = set(acc_id_changed_not_updatable)
        acc_id_unchanged = acc_id_present_in_current_and_remote - acc_id_changed_updatable.keys() - acc_id_changed_not_updatable

        acc_id_to_remove = acc_id_missing_in_remote | acc_id_changed_not_updatable
        acc_id_to_import = acc_id_missing_in_current | acc_id_changed_not_updatable

        logger.info(f'\n'
                    f'# Sequences from remote source: {len(acc_id_remote)}. Of which\n'
                    f'# {len(acc_id_unchanged)} present also locally and unchanged\n'
                    f'# {len(acc_id_missing_in_current)} never seen before.\n'
                    f'# {len(acc_id_changed_not_updatable)} have changes in the remote source that impact overlaps (strain/length)\n'  # TODO temporary info
                    f"# {len(acc_id_changed_updatable)} have changes that don't affect overlaps\n"
                    f'# Sequences from local source: {len(acc_id_current)}. Of which\n'
                    f'# {len(acc_id_missing_in_remote)} are missing from remote and must be removed from local.\n'
                    f'# In conclusion: {len(acc_id_to_remove)} sequences will be removed because missing or changed strain/length in remote\n'
                    f'# {len(acc_id_to_import)} sequences will be imported because novel or changed strain/length in remote\n'
                    f'# {len(acc_id_changed_updatable)} sequence will be updated with changes from remote.')

        logger.info(f"distribution of updates: {changes_distribution}")

        return acc_id_to_remove, acc_id_to_import, acc_id_changed_updatable

