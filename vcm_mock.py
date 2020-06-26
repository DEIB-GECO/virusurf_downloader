from typing import List
# noinspection PyPackageRequirements
from Bio import Entrez
from lxml import etree
from tqdm import tqdm

from database import ExperimentType, SequencingProject, Virus, HostSample, Sequence, Annotation
from locations import *
from data_sources.ncbi_sars_cov_2.sample import NCBISarsCov2Sample
from xml_helper import *
import string
import random
from loguru import logger

# Set of MOCK methods that create the corresponding rows in the Virus Conceptual Model database tables
logger.warning('You are using a mocked version of the VCM: You can use it to develop the DB and test the flow of operations'
               ' but the data inserted are just meaningless examples.\n'
               'If you want to use real data, import vcm.py instead.')

#   #############################    VCM    ############################
def create_or_get_virus(session, tax_tree):
    virus = session.query(Virus).filter(Virus.taxon_id == 101).one_or_none()

    if not virus:
        virus = Virus(taxon_id=101,
                      taxon_name='a_taxon_name',
                      family='a_family',
                      sub_family='a_subfamily',
                      genus='a_genus',
                      species='a_specie',
                      equivalent_list='equivalent_names',
                      molecule_type='a_molecule_type',
                      is_single_stranded=True,
                      is_positive_stranded=True
                      )
        session.add(virus)
        session.flush()
        print('virus added')
    else:
        print('virus with taxon id 101 already_present')
    return virus


def create_or_get_experiment(session, sample: NCBISarsCov2Sample):
    sequencing_technology = 'a_seq_technology',
    assembly_method = 'an_assembly_method',
    coverage = 'some_coverage'

    experiment = session.query(ExperimentType).filter(ExperimentType.sequencing_technology == sequencing_technology,
                                                      ExperimentType.assembly_method == assembly_method,
                                                      ExperimentType.coverage == coverage).one_or_none()
    if not experiment:
        experiment = ExperimentType(
            sequencing_technology=sequencing_technology,
            assembly_method=assembly_method,
            coverage=coverage)
        session.add(experiment)
        session.flush()
        print('experiment type added')
    else:
        print('experiment type already existing')
    return experiment


def create_or_get_sequencing_project(session, sample: NCBISarsCov2Sample):
    sequencing_lab = 'a_laboratory'
    submission_date = 'a_date'
    bioproject_id = 'a_bioproject_id'
    database_source = 'a_db_source'

    sequencing_project = session.query(SequencingProject).filter(SequencingProject.sequencing_lab == sequencing_lab,
                                                                 SequencingProject.submission_date == submission_date,
                                                                 SequencingProject.bioproject_id == bioproject_id,
                                                                 SequencingProject.database_source == database_source
                                                                 ).one_or_none()

    if not sequencing_project:
        sequencing_project = SequencingProject(sequencing_lab=sequencing_lab,
                                               submission_date=submission_date,
                                               bioproject_id=bioproject_id,
                                               database_source=database_source)
        session.add(sequencing_project)
        session.flush()
        print('sequencing project added')
    else:
        print('sequencing project already existing')
    return sequencing_project


def create_or_get_host_sample(session, sample: NCBISarsCov2Sample):

    host_taxon_id = 222
    host_taxon_name = 'a_name'
    collection_date = 'a_date'
    isolation_source = 'a_source'
    originating_lab = 'a_lab'
    country = 'a_country'
    region = 'a_region'
    geo_group = 'a_geo_group'
    age = 40
    gender = 'a_gender'

    host_sample = session.query(HostSample).filter(HostSample.host_taxon_id == host_taxon_id,
                                                   HostSample.host_taxon_name == host_taxon_name,

                                                   HostSample.collection_date == collection_date,
                                                   HostSample.isolation_source == isolation_source,

                                                   HostSample.originating_lab == originating_lab,
                                                   HostSample.country == country,
                                                   HostSample.region == region,
                                                   HostSample.geo_group == geo_group,
                                                   HostSample.age == age,
                                                   HostSample.gender == gender,
                                                   ).one_or_none()

    if not host_sample:
        host_sample = HostSample(host_taxon_id=host_taxon_id,
                                 host_taxon_name=host_taxon_name,

                                 collection_date=collection_date,
                                 isolation_source=isolation_source,

                                 originating_lab=originating_lab,
                                 country=country,
                                 region=region,
                                 geo_group=geo_group,
                                 age=age,
                                 gender=gender,
                                 )
        session.add(host_sample)
        session.flush()
        print('host sample added')
    else:
        print('host_sample_already existing')
    return host_sample


def create_or_get_sequence(session, virus_sample, virus_id: int, experiment: ExperimentType, host_sample: HostSample, sequencing_project: SequencingProject):
    # data from sample
    accession_id = 'primary_accession_id'
    alternative_accession_id = 'alternativea_accession_id'
    strain_name = 'a_strain'
    is_reference = False
    is_complete = False
    nucleotide_sequence = 'a_nuc_sequence'
    strand = '+'
    length = '10'
    gc_percentage = '0.44'
    # foreign keys
    experiment_type_id = experiment.experiment_type_id
    sequencing_project_id = sequencing_project.sequencing_project_id
    virus_id = virus_id
    host_sample_id = host_sample.host_sample_id

    sequence = session.query(Sequence).filter(Sequence.accession_id == accession_id,
                                              Sequence.alternative_accession_id == alternative_accession_id,
                                              Sequence.strain_name == strain_name,
                                              Sequence.is_reference == is_reference,
                                              Sequence.is_complete == is_complete,
                                              Sequence.nucleotide_sequence == nucleotide_sequence,
                                              Sequence.strand == strand,
                                              Sequence.length == length,
                                              Sequence.gc_percentage == gc_percentage,
                                              Sequence.experiment_type_id == experiment_type_id,
                                              Sequence.sequencing_project_id == sequencing_project_id,
                                              Sequence.virus_id == virus_id,
                                              Sequence.host_sample_id == host_sample_id).one_or_none()
    if not sequence:
        sequence = Sequence(accession_id=accession_id,
                            alternative_accession_id=alternative_accession_id,
                            strain_name=strain_name,
                            is_reference=is_reference,
                            is_complete=is_complete,
                            nucleotide_sequence=nucleotide_sequence,
                            strand=strand,
                            length=length,
                            gc_percentage=gc_percentage,
                            experiment_type_id=experiment_type_id,
                            sequencing_project_id=sequencing_project_id,
                            virus_id=virus_id,
                            host_sample_id=host_sample_id)
        session.add(sequence)
        session.flush()
        print('sequence added')
    else:
        print('sequence already existing')
    return sequence


def create_or_get_annotation(session, sample: NCBISarsCov2Sample, sequence: Sequence):
    # sample data
    feature_type = 'a_type'
    start = 1
    stop = 2
    gene_name = 'a_name'
    product = 'a_product'
    external_reference = 'a_reference'
    # foreign keys
    sequence_id = sequence.sequence_id

    annotation = session.query(Annotation).filter(Annotation.sequence_id == sequence_id,
                                                  Annotation.feature_type == feature_type,
                                                  Annotation.start ==  start,
                                                  Annotation.stop == stop,
                                                  Annotation.gene_name == gene_name,
                                                  Annotation.product == product,
                                                  Annotation.external_reference == external_reference).one_or_none()
    if not annotation:
        annotation = Annotation(sequence_id=sequence_id,
                                feature_type=feature_type,
                                start=start,
                                stop=stop,
                                gene_name=gene_name,
                                product=product,
                                external_reference=external_reference)
        session.add(annotation)
        session.flush()
        print('annotation added')
    else:
        print('annotation already existing')
    return [annotation]


#   ##############################      HELPER METHODS  #################à#
def download_virus_taxonomy_as_xml(taxon_id):
    # write taxonomy tree
    local_file_path = f"{local_folder_taxonomy}/{taxon_id}.xml"
    if not os.path.exists(local_file_path):
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                open(local_file_path, 'w') as f:
            for line in handle:
                f.write(line)
    tax_tree = etree.parse(local_file_path, parser=etree.XMLParser(remove_blank_text=True))
    return tax_tree


def get_virus_sample_accession_ids(virus_specie_taxon_id: int) -> (int, List[int]):
    """
    For a virus specie, it returns the accession ids of the reference sample (refseq) and a list of ids for the non-reference
    samples.
    """
    # get reference sample accession id
    handle = Entrez.esearch(db="nuccore", term=f"(txid{virus_specie_taxon_id}[Organism]) AND srcdb_refseq[Properties]")
    response = Entrez.read(handle)
    handle.close()
    # Example response
    # {
    #     "Count": "1",                       <-- number of total records matching the query
    #     "RetMax": "1",
    #     "RetStart": "0",
    #     "IdList": ["1798174254"],            <-- accession id of refseq
    #     "TranslationSet": [],
    #     "TranslationStack": [
    #         {"Term": "txid2697049[Organism]", "Field": "Organism", "Count": "5511", "Explode": "Y"},
    #         {"Term": "srcdb_refseq[Properties]", "Field": "Properties", "Count": "66995641", "Explode": "N"},
    #         "AND"
    #     ],
    #     "QueryTranslation": "txid2697049[Organism] AND srcdb_refseq[Properties]"
    # }
    assert int(response['Count']) == 1, \
        "no reference sample found or multiple RefSeqs" + "please check: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch"
    refseq_accession_id = response['IdList'][0]

    # get non-refseq accession ids
    #   total number of sequences
    with Entrez.esearch(db="nuccore",
                        term=f"(txid{virus_specie_taxon_id}[Organism]) NOT srcdb_refseq[Properties]",
                        rettype='count') as handle:
        response = Entrez.read(handle)
        total_records = int(response['Count'])
    #   do pagination
    non_refseq_accessions_ids = list()
    RECORDS_PER_PAGE = 1000
    page_number = 0
    with tqdm(total=total_records) as progress_bar:
        while total_records > page_number*RECORDS_PER_PAGE:
            with Entrez.esearch(db="nuccore", term=f"(txid{virus_specie_taxon_id}[Organism]) NOT srcdb_refseq[Properties]",
                                retmax=RECORDS_PER_PAGE, retstart=page_number*RECORDS_PER_PAGE) as handle:
                response = Entrez.read(handle)
            for x in response['IdList']:
                non_refseq_accessions_ids.append(int(x))
            page_number += 1
            progress_bar.update(RECORDS_PER_PAGE if page_number*RECORDS_PER_PAGE < total_records else total_records-((page_number-1)*RECORDS_PER_PAGE))
    assert len(non_refseq_accessions_ids) == total_records, 'Some of the non-refseq accession ids were not correctly downloaded'
    return refseq_accession_id, non_refseq_accessions_ids


# TODO fix statement without effect
def merge_intervals(e):
    intervals = e.xpath(".//INSDInterval")
    intervals2 = []
    for i in intervals:
        start = int(text_at_node(i, './/INSDInterval_from'))
        stop = int(text_at_node(i, './/INSDInterval_to'))
        intervals2.append((start, stop))

    if intervals:
        min_start = min(x[0] for x in intervals2)
        max_stop = max(x[1] for x in intervals2)
        return min_start, max_stop
    else:
        None
    return len(intervals)


def random_string():
    return ''.join(random.choices(string.ascii_uppercase, k=4))
