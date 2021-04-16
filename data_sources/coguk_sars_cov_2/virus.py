import os
import shutil
from collections import OrderedDict
from data_sources.ncbi_any_virus.ncbi_importer import AnyNCBIVNucSample
from data_sources.ncbi_services import download_ncbi_taxonomy_as_xml, get_samples_accession_ids, download_or_get_ncbi_sample_as_xml
from data_sources.ncbi_any_virus.settings import known_settings as ncbi_known_settings
from lxml import etree
import wget as wget
from tqdm import tqdm
from typing import Generator, Optional, Collection
from loguru import logger
from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from db_config.database import Sequence, HostSample, SequencingProject, try_py_function, NucleotideSequence
from locations import get_local_folder_for, FileType, remove_file
from xml_helper import text_at_node
from urllib.request import Request, urlopen
import nuc_aa_pipeline


class COGUKSarsCov2:

    name = 'coguk_sars_cov_2'
    sequence_file_url = 'https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_all.fasta'
    metadata_file_url = 'https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv'

    def __init__(self):
        logger.info(f'importing virus {self.name} using NCBI SC2 for taxonomy data')
        # fetch taxonomy data from NCBI
        taxonomy_file_path = download_ncbi_taxonomy_as_xml(
            get_local_folder_for(source_name=self.name, _type=FileType.TaxonomyData),
            self.taxon_id())
        try:
            self.tax_tree = etree.parse(taxonomy_file_path, parser=etree.XMLParser(remove_blank_text=True))
        except etree.XMLSyntaxError as e:  # happens on AWS if for some reason the downloaded file is corrupted
            remove_file(taxonomy_file_path)
            ncbi_sc2_taxonomy_dir = get_local_folder_for(
                source_name=ncbi_known_settings["sars_cov_2"]["generated_dir_name"], _type=FileType.TaxonomyData)
            alternative_taxonomy_path = ncbi_sc2_taxonomy_dir + f"{ncbi_known_settings['sars_cov_2']['virus_taxon_id']}.xml"
            if os.path.exists(alternative_taxonomy_path):
                shutil.copyfile(alternative_taxonomy_path, taxonomy_file_path)
                self.tax_tree = etree.parse(taxonomy_file_path, parser=etree.XMLParser(remove_blank_text=True))
            else:
                logger.error(f"Taxonomy file of SARS-CoV-2 was empty. Attempt to use the one from {ncbi_sc2_taxonomy_dir} "
                             f"failed because the filed doesn't exist. Can't proceed.")
                raise e
        # fetch latest source data
        download_dir = get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)
        self.sequence_file_path, self.metadata_file_path = download_or_get_sample_data(download_dir)

    def taxon_id(self):
        return 2697049

    def taxon_name(self):
        return text_at_node(self.tax_tree, './Taxon/ScientificName')

    def family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName')

    def sub_family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName')

    def genus(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName')

    def species(self):
        # species_taxon_id = text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/TaxId')
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName',
                            mandatory=False)

    def equivalent_names(self):
        genbank_acronym = text_at_node(self.tax_tree, './/GenbankAcronym', mandatory=False)
        equivalent_names = self.tax_tree.xpath('.//EquivalentName')
        equivalent_names = [x.text for x in equivalent_names]
        if genbank_acronym:
            equivalent_names.insert(0, genbank_acronym)
        equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
        equivalent_names = ", ".join(equivalent_names)
        return equivalent_names

    def molecule_type(self):
        return 'RNA'

    def is_single_stranded(self):
        return True

    def is_positive_stranded(self):
        return True

    @staticmethod
    def get_sequences_in_current_data() -> dict:
        def do(session):
            all_sequences = dict()
            for db_items in session.query(Sequence, HostSample, NucleotideSequence) \
                    .select_from(Sequence, HostSample, SequencingProject, NucleotideSequence) \
                    .filter(Sequence.strain_name.isnot(None),
                            Sequence.host_sample_id == HostSample.host_sample_id,
                            Sequence.sequence_id == NucleotideSequence.sequence_id,
                            Sequence.sequencing_project_id == SequencingProject.sequencing_project_id,
                            SequencingProject.database_source == 'COG-UK') \
                    .yield_per(100):
                source_seq, source_host, nuc_seq_db_obj = db_items
                all_sequences[source_seq.accession_id] = (source_seq.length,
                                                          source_seq.gc_percentage, source_seq.n_percentage,
                                                          source_host.collection_date,
                                                          source_host.region, source_host.country,
                                                          nuc_seq_db_obj.nucleotide_sequence,
                                                          source_seq.lineage)
            return all_sequences

        return try_py_function(do)

    def get_sequences_of_updated_source(self, filter_accession_ids: Optional[Collection] = None):
        meta = {}
        with open(self.metadata_file_path, mode='r') as metadata_file:
            header = metadata_file.readline()  # skip header
            if not header.startswith("sequence_name,country,adm1,pillar_2,sample_date,epi_week,lineage"):
                raise AssertionError(f"Unexpected metadata columns in file {self.metadata_file_path}.\n"
                                     f"Check the compatibility of this new metadata file before importing.")
            for line in metadata_file:
                try:
                    key, content = line.split(sep=',', maxsplit=1)  # key = <sequence name>, content = <other metadata>
                    meta[key] = content.rstrip()
                except:
                    logger.error(
                        f'Unable to parse the following line from the metadata file {self.metadata_file_path}:\n\t'
                        f'{line}')

        accession_ids_of_interest = filter_accession_ids if filter_accession_ids is not None else meta.keys()

        if len(accession_ids_of_interest) == 0:
            yield from ()
        else:
            # scan the multi-fasta file
            nuc_sequences_with_errors = 0
            with open(self.sequence_file_path, mode='r') as seq_file:
                progress = tqdm(total=len(accession_ids_of_interest))
                while True:
                    sample_key = seq_file.readline()
                    sample_sequence = seq_file.readline()
                    # exit loop when sequences are finished
                    if not sample_sequence or not sample_key:
                        break
                    else:
                        # in the sequence file, the strain name is preceded by a '>', while in the metadata file not
                        sample_key = sample_key[1:].rstrip()
                        # ignore sequences not in the accession_ids_of_interest list
                        if sample_key not in accession_ids_of_interest:
                            continue
                        else:
                            progress.update()
                            sample_sequence = sample_sequence.rstrip().replace("?", "n")
                            if not nuc_aa_pipeline.is_valid_sequence(sample_sequence):
                                nuc_sequences_with_errors += 1
                                continue
                            sample = {
                                COGUKSarsCov2Sample.STRAIN_NAME: sample_key,
                                COGUKSarsCov2Sample.NUC_SEQUENCE: sample_sequence
                            }
                            try:
                                sample[COGUKSarsCov2Sample.METADATA_RAW_STRING] = meta[sample_key]
                            except KeyError:
                                logger.error(
                                    f'Found a sequence without a paired metadata string. Foreign key was {sample_key}')
                            try:
                                sample_obj = COGUKSarsCov2Sample(sample)
                            except Exception:
                                logger.exception(
                                    f"Sample {sample_key} skipped due to an error while parsing the input data.")
                                continue
                            yield sample_obj
                if nuc_sequences_with_errors > 0:
                    logger.warning(f"{nuc_sequences_with_errors} were rejected because they include one of the following "
                                   f"invalid characters: {nuc_aa_pipeline.not_allowed_chars}")
                logger.debug(f"known nucleotides observed in this set of sequences: {nuc_aa_pipeline.used_characters}")

    def deltas(self):
        """
        Returns the set of Sequence.accession_id to remove from current database, and the set of Sequence.accession_id
        that should be imported from the source. Note that because sequences in the current data may have been updated,
        the two sets are not necessarily disjoint. So first remove the sequences to be removed and only then
         insert the new ones.
         """
        logger.info("Reading updated data...")
        acc_id_remote = set([x.primary_accession_number() for x in self.get_sequences_of_updated_source()])  # all the sequences from remote
        logger.info('Reading current data...')
        current_data = self.get_sequences_in_current_data()
        acc_id_current = set(current_data.keys())

        acc_id_missing_in_remote = acc_id_current - acc_id_remote
        acc_id_missing_in_current = acc_id_remote - acc_id_current

        acc_id_present_in_current_and_remote = acc_id_current & acc_id_remote
        acc_id_changed = []

        # to compare sequences present in both, have to scan the file of remote sequences
        logger.info('Comparing current data with that from source...')
        for new_sequence in self.get_sequences_of_updated_source(filter_accession_ids=acc_id_present_in_current_and_remote):
            acc_id = new_sequence.primary_accession_number()
            try:
                current_sequence_data = current_data[acc_id]
                # if the gisaid id is the same, compare metadata to detect if the sample changed over time
                if current_sequence_data[0] != new_sequence.length() \
                        or current_sequence_data[1] != new_sequence.gc_percent() \
                        or current_sequence_data[2] != new_sequence.n_percent() \
                        or current_sequence_data[3] != new_sequence.collection_date()[0] \
                        or (current_sequence_data[4], current_sequence_data[5]) \
                        != new_sequence.province__region__country__geo_group()[1:3] \
                        or current_sequence_data[7] != new_sequence.lineage() \
                        or current_sequence_data[6] != new_sequence.nucleotide_sequence():
                    acc_id_changed.append(acc_id)
            except KeyError:
                pass  # the accession id is not present in current data. it's a new sequence

        # compute additional sets
        acc_id_changed = set(acc_id_changed)
        acc_id_unchanged = acc_id_present_in_current_and_remote - acc_id_changed

        if len(acc_id_missing_in_remote) > 5000:
            logger.warning(f'The number of local sequences missing from remote is really HIGH '
                           f'({len(acc_id_missing_in_remote)}). It has ben cut to 5K max.')
            for i in range(len(acc_id_missing_in_remote) - 5000):
                acc_id_missing_in_remote.pop()

        acc_id_to_remove = acc_id_missing_in_remote | acc_id_changed
        acc_id_to_import = acc_id_missing_in_current | acc_id_changed

        logger.info(f'\n'
                    f'# Sequences from remote source: {len(acc_id_remote)}. Of which\n'
                    f'# {len(acc_id_unchanged)} present also locally and unchanged\n'
                    f'# {len(acc_id_missing_in_current)} never seen before.\n'
                    f'# {len(acc_id_changed)} have different attributes in the remote source\n'
                    f'# Sequences from local source: {len(acc_id_current)}. Of which\n'
                    f'# {len(acc_id_missing_in_remote)} are missing from remote and must be removed from local.\n'
                    f'# In conclusion: {len(acc_id_to_remove)} sequences will be removed because missing or changed in remote\n'
                    f'# and {len(acc_id_to_import)} sequences will be imported because novel or changed in remote.')

        return acc_id_to_remove, acc_id_to_import

    # following function is a copy of deltas with some debug prints
    # def deltas(self):
    #     """
    #     Returns the set of Sequence.accession_id to remove from current database, and the set of Sequence.accession_id
    #     that should be imported from the source. Note that because sequences in the current data may have been updated,
    #     the two sets are not necessarily disjoint. So first remove the sequences to be removed and only then
    #      insert the new ones.
    #      """
    #     logger.info("Reading updated data...")
    #     acc_id_remote = set([x.primary_accession_number() for x in
    #                          self.get_sequences_of_updated_source()])  # all the sequences from remote
    #     logger.info('Reading current data...')
    #     current_data = self.get_sequences_in_current_data()
    #     acc_id_current = set(current_data.keys())
    #
    #     acc_id_missing_in_remote = acc_id_current - acc_id_remote
    #     acc_id_missing_in_current = acc_id_remote - acc_id_current
    #
    #     acc_id_present_in_current_and_remote = acc_id_current & acc_id_remote
    #     acc_id_changed = []
    #     acc_id_changed_sequence = []
    #     counter = 0
    #
    #     # to compare sequences present in both, have to scan the file of remote sequences
    #     logger.info('Comparing current data with that from source...')
    #     for new_sequence in self.get_sequences_of_updated_source(
    #             filter_accession_ids=acc_id_present_in_current_and_remote):
    #         acc_id = new_sequence.primary_accession_number()
    #         if counter > 10:
    #             break
    #         try:
    #             current_sequence_data = current_data[acc_id]
    #             # if the gisaid id is the same, compare metadata to detect if the sample changed over time
    #             if current_sequence_data[0] != new_sequence.length() \
    #                     or current_sequence_data[1] != new_sequence.gc_percent() \
    #                     or current_sequence_data[2] != new_sequence.n_percent() \
    #                     or current_sequence_data[3] != new_sequence.collection_date()[0] \
    #                     or (current_sequence_data[4], current_sequence_data[5]) \
    #                     != new_sequence.province__region__country__geo_group()[1:3]:
    #                 acc_id_changed.append(acc_id)
    #                 counter += 1
    #
    #                 logger.trace(
    #                     f"current data of {acc_id}: {current_sequence_data[0]}, {current_sequence_data[1]}, {current_sequence_data[2]}, {current_sequence_data[3]}, {current_sequence_data[4]}, {current_sequence_data[5]}")
    #                 logger.trace(
    #                     f"updated data of {acc_id}: {new_sequence.length()}, {new_sequence.gc_percent()}, {new_sequence.n_percent()}, {new_sequence.collection_date()[0]}, {new_sequence.province__region__country__geo_group()[1]}, {new_sequence.province__region__country__geo_group()[2]}")
    #
    #                 if current_sequence_data[0] != new_sequence.length():
    #                     logger.trace("is length")
    #                 elif current_sequence_data[1] != new_sequence.gc_percent():
    #                     logger.trace("is gc percent")
    #                 elif current_sequence_data[2] != new_sequence.n_percent():
    #                     logger.trace("is n percent")
    #                 elif current_sequence_data[3] != new_sequence.collection_date()[0]:
    #                     logger.trace("is coll date")
    #                 elif [current_sequence_data[4],
    #                       current_sequence_data[5]] != new_sequence.province__region__country__geo_group()[1:3]:
    #                     logger.trace("is country or region")
    #
    #             elif current_sequence_data[6] != new_sequence.nucleotide_sequence():
    #                 acc_id_changed_sequence.append(acc_id)
    #
    #                 logger.trace(f"{acc_id} has a different sequence")
    #
    #         except KeyError:
    #             pass  # the accession id is not present in current data. it's a new sequence
    #
    #     # compute additional sets
    #     acc_id_changed = set(acc_id_changed)
    #     acc_id_changed_sequence = set(acc_id_changed_sequence)
    #
    #     acc_id_unchanged = acc_id_present_in_current_and_remote - acc_id_changed - acc_id_changed_sequence
    #
    #     acc_id_to_remove = acc_id_missing_in_remote | acc_id_changed | acc_id_changed_sequence
    #     acc_id_to_import = acc_id_missing_in_current | acc_id_changed | acc_id_changed_sequence
    #
    #     logger.info(f'\n'
    #                 f'# Sequences from remote source: {len(acc_id_remote)}. Of which\n'
    #                 f'# {len(acc_id_unchanged)} present also locally and unchanged\n'
    #                 f'# {len(acc_id_missing_in_current)} never seen before.\n'
    #                 f'# {len(acc_id_changed)} have different metadata attributes in the remote source\n'
    #                 f'# {len(acc_id_changed_sequence)} have different sequence in the remote source\n'
    #                 f'# Sequences from local source: {len(acc_id_current)}. Of which\n'
    #                 f'# {len(acc_id_missing_in_remote)} are missing from remote and must be removed from local.\n'
    #                 f'# In conclusion: {len(acc_id_to_remove)} sequences will be removed because missing or changed metdata in remote\n'
    #                 f'# and {len(acc_id_to_import)} sequences will be imported because novel or changed metadtata in remote.')
    #
    #     return acc_id_to_remove, acc_id_to_import

    def reference_sequence(self) -> str:
        if not hasattr(self, 'reference_sample'):
            # import a reference sequence from a different dataset that we'll use to call nucleotide variants
            # get reference sample accession id
            ncbi_reference_sample_query = ncbi_known_settings["sars_cov_2"]["reference_sample_query"]
            reference_accession_id = get_samples_accession_ids(ncbi_reference_sample_query)
            assert len(reference_accession_id) == 1, \
                "no reference sample found or multiple RefSeqs. Please correct the query used on NCBI nuccore"
            # download file as XML
            reference_seq_file_path = download_or_get_ncbi_sample_as_xml(
                get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData),
                reference_accession_id[0]
            )
            # parse file and cache the object
            self.reference_sample = AnyNCBIVNucSample(reference_seq_file_path, reference_accession_id[0])
        return self.reference_sample.nucleotide_sequence()


def download_or_get_sample_data(containing_directory: str) -> (str, str):
    """
    :return: the local file path of the downloaded sequence and metadata files.
    """
    def get_download_size(url) -> Optional[int]:
        """
        Returns the size of the downloadable resource in bytes supported by the protocol
        of the downloadable resource; None otherwise.
        """
        req = Request(url=url, method='HEAD')
        f = urlopen(req)
        if int(f.status) == 200:
            return int(f.headers['Content-Length'])
        else:
            return None

    def get_local_file_size(path) -> Optional[int]:
        """
        Returns the size of the local file in bytes if it exists; None otherwise.
        """
        if os.path.exists(path):
            return os.stat(path).st_size
        else:
            return None

    def download_coguk_data():
        # make sure the output path does not exist already, or wget assigns a trailing number to it
        logger.info(f'downloading sample sequences for COG-UK data from {COGUKSarsCov2.sequence_file_url} ...')
        wget.download(COGUKSarsCov2.sequence_file_url, sequence_local_file_path)
        logger.info(f'\ndownloading sample metadata for COG-UK data from {COGUKSarsCov2.metadata_file_url} ...')
        wget.download(COGUKSarsCov2.metadata_file_url, metadata_local_file_path)
        logger.info('\n')

    sequence_local_file_path = containing_directory + COGUKSarsCov2.sequence_file_url.rsplit('/', maxsplit=1)[1]
    metadata_local_file_path = containing_directory + COGUKSarsCov2.metadata_file_url.rsplit('/', maxsplit=1)[1]
    if not os.path.exists(sequence_local_file_path) or not os.path.exists(metadata_local_file_path):
        download_coguk_data()
    else:
        # compare size of local with remote ones
        if get_download_size(COGUKSarsCov2.sequence_file_url) != get_local_file_size(sequence_local_file_path) or \
                get_download_size(COGUKSarsCov2.metadata_file_url) != get_local_file_size(metadata_local_file_path):
            remove_file(sequence_local_file_path)
            remove_file(metadata_local_file_path)
            download_coguk_data()
    return sequence_local_file_path, metadata_local_file_path
