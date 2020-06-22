from collections import OrderedDict
from typing import List
# noinspection PyPackageRequirements
from Bio import Entrez
from lxml import etree
from database_tom import ExperimentType, SequencingProject, Virus, HostSample, Sequence, Annotation, Variant
from locations import *
from tqdm import tqdm
from virus_sample import VirusSample
from xml_helper import *

# Set of methods that create the corresponding rows in the Virus Conceptual Model database tables


#   #############################    VCM    ############################
def create_or_get_virus(session, tax_tree):
    taxon_id = int(text_at_node(tax_tree, './Taxon/TaxId'))
    scientific_name = text_at_node(tax_tree, './Taxon/ScientificName')

    family = text_at_node(tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName')
    subfamily = text_at_node(tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName')
    genus = text_at_node(tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName')
    species = text_at_node(tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName')
    #     species_taxon_id = text_at_node(tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/TaxId')

    genbank_acronym = text_at_node(tax_tree, './/GenbankAcronym')
    equivalent_names = tax_tree.xpath('.//EquivalentName')
    equivalent_names = [x.text for x in equivalent_names]
    if genbank_acronym:
        equivalent_names.insert(0, genbank_acronym)
    equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
    equivalent_names = ", ".join(equivalent_names)

    print(family, subfamily, genus, species, genbank_acronym, equivalent_names)

    molecule_type = 'RNA'
    is_single_stranded = True
    is_positive_stranded = True

    virus = session.query(Virus).filter(Virus.taxon_id == taxon_id).one_or_none()

    if not virus:
        #         print("not exists")
        virus = Virus(taxon_id=taxon_id,
                      taxon_name=scientific_name,
                      family=family,
                      sub_family=subfamily,
                      genus=genus,
                      species=species,
                      equivalent_list=equivalent_names,
                      molecule_type=molecule_type,
                      is_single_stranded=is_single_stranded,
                      is_positive_stranded=is_positive_stranded
                      )
        session.add(virus)
        session.flush()
    return virus


def create_or_get_experiment(session, sample: VirusSample):
    sequencing_technology = sample.sequencing_technology()
    assembly_method = sample.assembly_method()
    coverage = sample.coverage()

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
    return experiment


def create_or_get_sequencing_project(session, sample: VirusSample):

    submitted = sample.submitted()
    assert submitted == "Submitted ", "Journal value submitted"
    submission_date = sample.submission_date()
    sequencing_lab = sample.sequencing_lab()
    bioproject_id = sample.bioproject_id()
    database_source = sample.database_source()

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
    return sequencing_project


def create_or_get_host_sample(session, sample: VirusSample):

    host_taxon_name = sample.taxon_name()
    host_taxon_id = sample.taxon_id()
    gender = sample.gender()
    age = sample.age()

    originating_lab = sample.originating_lab()
    collection_date = sample.collection_date()
    isolation_source = sample.isolation_source()

    country, region, geo_group = sample.country__region__geo_group()

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
        #         print("not exists")
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
    return host_sample


def create_or_get_sequence(session, virus_sample, virus_id: int, experiment: ExperimentType, host_sample: HostSample, sequencing_project: SequencingProject):
    # data from sample
    accession_id = virus_sample.primary_accession_number()
    alternative_accession_id = str(virus_sample.alternative_accession_number())
    strain_name = virus_sample.strain()
    is_reference = virus_sample.is_reference()
    is_complete = virus_sample.is_complete()
    nucleotide_sequence = virus_sample.nucleotide_sequence()
    strand = virus_sample.strand()
    length = virus_sample.length()
    gc_percentage = virus_sample.gc_percent()
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
    return sequence


def create_or_get_annotation(session, sample: VirusSample, sequence: Sequence):
    annotations = []
    for start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence in sample.annotations():
        if feature_type != 'source':
            annotation = session.query(Annotation).filter(Annotation.feature_type == feature_type,
                                                          Annotation.start == start,
                                                          Annotation.stop == stop,
                                                          Annotation.gene_name == gene_name,
                                                          Annotation.product == product,
                                                          Annotation.external_reference == db_xref_merged,
                                                          Annotation.sequence_id == sequence.sequence_id,
                                                          Annotation.aminoacid_sequence == amino_acid_sequence).one_or_none()
            if not annotation:
                annotation = Annotation(feature_type=feature_type,
                                        start=start,
                                        stop=stop,
                                        gene_name=gene_name,
                                        product=product,
                                        external_reference=db_xref_merged,
                                        sequence_id=sequence.sequence_id,
                                        aminoacid_sequence=amino_acid_sequence)
                session.add(annotation)
            annotations.append(annotation)
    session.flush()
    return annotations


def create_or_get_nucleotide_variants(session, sample: VirusSample, sequence: Sequence, aligner):
    if aligner is None:
        raise Exception("Can't create nucleotide variants because the aligner is None")
    sequence_variants: list = session.query(Variant).filter(Variant.sequence_id == sequence.sequence_id).all()
    if not sequence_variants:
        for (sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type) in sample.call_variants(aligner):
            row = Variant(sequence_id=sequence.sequence_id,
                          sequence_original=sequence_original,
                          sequence_alternative=sequence_alternative,
                          start_original=start_original,
                          start_alternative=start_alternative,
                          variant_length=variant_length,
                          variant_type=variant_type)
            session.add(row)
            sequence_variants.append(row)
        session.flush()
    return sequence_variants


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
