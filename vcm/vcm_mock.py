from typing import List, Tuple

import sqlalchemy
from sqlalchemy import cast
from vcm.vcm import DBCache
from data_sources.virus_sample import VirusSample
from db_config.database import AminoAcidVariant, ExperimentType, SequencingProject, Virus, HostSample, Sequence, Annotation, \
    Epitope, HostSpecie, EpitopeFragment, NucleotideSequence, AnnotationSequence, PipelineEvent
from xml_helper import *
import string
import random
from loguru import logger


logger.warning('You are using a mocked version of the VCM: You can use it to test the access and creation of '
               'SQLAlchemy ORM objects but no data is actually inserted into the configured DB.\n'
               'If you want to insert data, import vcm.py instead.')

cache_host_specie = DBCache.host_specie
cache_host_sample = DBCache.host_sample
cache_experiment_type = DBCache.experiment_type
cache_sequencing_project = DBCache.sequencing_project
cache_virus = DBCache.virus

#   #############################    VCM    ############################
def create_or_get_virus(session, a_virus):
    # virus = session.query(Virus).filter(Virus.taxon_id == a_virus.taxon_id()).one_or_none()
    # if not virus:
    virus = Virus(taxon_id=a_virus.taxon_id(),
                  taxon_name=a_virus.taxon_name(),
                  family=a_virus.family(),
                  sub_family=a_virus.sub_family(),
                  genus=a_virus.genus(),
                  species=a_virus.species(),
                  equivalent_list=a_virus.equivalent_names(),
                  molecule_type=a_virus.molecule_type(),
                  is_single_stranded=a_virus.is_single_stranded(),
                  is_positive_stranded=a_virus.is_positive_stranded()
                  )
    virus.virus_id = 1
    return virus.virus_id


def create_or_get_experiment(session, sample: VirusSample):
    sequencing_technology = sample.sequencing_technology()
    assembly_method = sample.assembly_method()
    coverage = sample.coverage()

    # experiment = session.query(ExperimentType).filter(ExperimentType.sequencing_technology == sequencing_technology,
    #                                                   ExperimentType.assembly_method == assembly_method,
    #                                                   ExperimentType.coverage == coverage).one_or_none()
    # if not experiment:
    experiment = ExperimentType(
        sequencing_technology=sequencing_technology,
        assembly_method=assembly_method,
        coverage=coverage)
    experiment.experiment_type_id = 1
    return experiment.experiment_type_id


def create_or_get_sequencing_project(session, sample: VirusSample):
    submission_date = sample.submission_date()
    sequencing_lab = sample.sequencing_lab()
    bioproject_id = sample.bioproject_id()
    database_source = sample.database_source()

    # sequencing_project = session.query(SequencingProject).filter(SequencingProject.sequencing_lab == sequencing_lab,
    #                                                              SequencingProject.submission_date == submission_date,
    #                                                              SequencingProject.bioproject_id == bioproject_id,
    #                                                              SequencingProject.database_source == database_source
    #                                                              ).one_or_none()
    #
    # if not sequencing_project:
    sequencing_project = SequencingProject(sequencing_lab=sequencing_lab,
                                           submission_date=submission_date,
                                           bioproject_id=bioproject_id,
                                           database_source=database_source)
    sequencing_project.sequencing_project_id = 1
    return sequencing_project.sequencing_project_id


def create_or_get_host_specie(session, sample: VirusSample) -> int:
    host_taxon_name = sample.host_taxon_name()
    if host_taxon_name:
        host_taxon_name = host_taxon_name.lower()
    host_taxon_id = sample.host_taxon_id()
    host_specie = HostSpecie(host_taxon_id=host_taxon_id,
                             host_taxon_name=host_taxon_name)

    # logger.info(f'ACC.ID: {sample.primary_accession_number()} - {sample.alternative_accession_number()} HOST_SPECIE: {host_taxon_name} - {host_taxon_id}')
    return 1


def create_or_get_host_specie_alt(session, organism_name: str, organism_ncbi_id: int):
    return 1


def create_or_get_host_sample(session, sample: VirusSample, host_specie_id: int) -> int:
    gender = sample.gender()
    age = sample.age()

    originating_lab = sample.originating_lab()
    collection_date, precision = sample.collection_date()
    isolation_source = sample.isolation_source()

    province, region, country, geo_group = sample.province__region__country__geo_group()

    # host_sample = session.query(HostSample).filter(HostSample.host_id == host_specie_id,
    #                                                HostSample.collection_date == collection_date,
    #                                                HostSample.isolation_source == isolation_source,
    #                                                HostSample.originating_lab == originating_lab,
    #                                                HostSample.province == province,
    #                                                HostSample.region == region,
    #                                                HostSample.country == country,
    #                                                HostSample.geo_group == geo_group,
    #                                                HostSample.age == age,
    #                                                HostSample.gender == gender,
    #                                                ).one_or_none()
    #
    # if not host_sample:
        #         print("not exists")
    host_sample = HostSample(host_id=host_specie_id,
                             collection_date=collection_date,
                             coll_date_precision=precision,
                             isolation_source=isolation_source,
                             originating_lab=originating_lab,
                             province=province,
                             region=region,
                             country=country,
                             geo_group=geo_group,
                             age=age,
                             gender=gender,
                             )
    host_sample.host_sample_id = 1
    logger.info(f'HOST_SAMPLE: coll_date:{collection_date} - coll-date_prec:{precision} - isol_s:{isolation_source} - '
                f'orig_lab:{originating_lab} - prov:{province} - reg:{region} - country:{country} - geo_g:{geo_group} -'
                f' age:{age} - gend:{gender}')
    return host_sample.host_sample_id


def create_and_get_sequence(session, virus_sample: VirusSample, virus_id: int, experiment_id, host_sample_id, sequencing_project_id):
    # data from sample
    accession_id = virus_sample.primary_accession_number()
    alternative_accession_id = virus_sample.alternative_accession_number()
    strain_name = virus_sample.strain()
    is_reference = virus_sample.is_reference()
    is_complete = virus_sample.is_complete()
    nucleotide_sequence = virus_sample.nucleotide_sequence()
    strand = virus_sample.strand()
    length = virus_sample.length()
    gc_percentage = virus_sample.gc_percent()
    n_percentage = virus_sample.n_percent()
    lineage = virus_sample.lineage()
    clade = virus_sample.clade()

    # sequence = session.query(Sequence).filter(Sequence.accession_id == accession_id,
    #                                           Sequence.alternative_accession_id == alternative_accession_id,
    #                                           Sequence.strain_name == strain_name,
    #                                           Sequence.is_reference == is_reference,
    #                                           Sequence.is_complete == is_complete,
    #                                           Sequence.n_percentage == n_percentage,
    #                                           Sequence.nucleotide_sequence == nucleotide_sequence,
    #                                           Sequence.strand == strand,
    #                                           Sequence.length == length,
    #                                           Sequence.gc_percentage == gc_percentage,
    #                                           Sequence.lineage == lineage,
    #                                           Sequence.clade == clade,
    #                                           Sequence.experiment_type_id == experiment_id,
    #                                           Sequence.sequencing_project_id == sequencing_project_id,
    #                                           Sequence.virus_id == virus_id,
    #                                           Sequence.host_sample_id == host_sample_id).one_or_none()
    # if not sequence:
    sequence = Sequence(accession_id=accession_id,
                        alternative_accession_id=alternative_accession_id,
                        strain_name=strain_name,
                        is_reference=is_reference,
                        is_complete=is_complete,
                        n_percentage=n_percentage,
                        strand=strand,
                        length=length,
                        gc_percentage=gc_percentage,
                        lineage=lineage,
                        clade=clade,
                        experiment_type_id=experiment_id,
                        sequencing_project_id=sequencing_project_id,
                        virus_id=virus_id,
                        host_sample_id=host_sample_id)
    sequence.sequence_id = 1
    nucleotide_sequence_db_obj = None
    if nucleotide_sequence:
        nucleotide_sequence_db_obj = NucleotideSequence(sequence_id=sequence.sequence_id,
                                                        nucleotide_sequence=nucleotide_sequence)
    return sequence, nucleotide_sequence_db_obj


def create_annotation_and_aa_variants(session, sample: VirusSample, sequence: Sequence, reference_sample: VirusSample):
    for start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, aa_variants in sample.annotations_and_amino_acid_variants(reference_sample):
        annotation = Annotation(start=start,
                                stop=stop,
                                gene_name=gene_name,
                                feature_type=feature_type,
                                product=product,
                                external_reference=db_xref_merged,
                                sequence_id=sequence.sequence_id)
        if aa_variants:
            for original, alternative, mutpos, mut_len, mut_type in aa_variants:
                aa_variant = AminoAcidVariant(annotation_id=annotation.annotation_id,
                                              sequence_aa_original=original,
                                              sequence_aa_alternative=alternative,
                                              start_aa_original=mutpos,
                                              variant_aa_length=mut_len,
                                              variant_aa_type=mut_type)


def create_annotation_and_amino_acid_variants(session, sequence_id, *args):
    pass
    # print(args)
    # gene_name, product, protein_id, feature_type, start, stop, nuc_seq, amino_acid_seq, aa_variants = args
    # annotation = session.query(Annotation).filter(Annotation.start == start,
    #                                               Annotation.stop == stop,
    #                                               Annotation.feature_type == feature_type,
    #                                               Annotation.gene_name == gene_name,
    #                                               Annotation.product == product,
    #                                               Annotation.external_reference == protein_id,
    #                                               Annotation.sequence_id == sequence_id,
    #                                               Annotation.annotation_nucleotide_sequence == nuc_seq,
    #                                               Annotation.aminoacid_sequence == amino_acid_seq).one_or_none()
    # if not annotation:
    #     annotation = Annotation(start=start,
    #                             stop=stop,
    #                             gene_name=gene_name,
    #                             feature_type=feature_type,
    #                             product=product,
    #                             sequence_id=sequence_id,
    #                             external_reference=protein_id,
    #                             annotation_nucleotide_sequence=nuc_seq,
    #                             aminoacid_sequence=amino_acid_seq)
    #     if aa_variants:
    #         for gen_name, protein, protein_id, start_pos, sequence_aa_original, sequence_aa_alternative, variant_aa_type in aa_variants:
    #             aa_variant = AminoAcidVariant(annotation_id=annotation.annotation_id,
    #                                           sequence_aa_original=sequence_aa_original,
    #                                           sequence_aa_alternative=sequence_aa_alternative,
    #                                           start_aa_original=start_pos,
    #                                           variant_aa_length=max(len(sequence_aa_original), len(sequence_aa_alternative)),
    #                                           variant_aa_type=variant_aa_type)


def create_nucleotide_variants_and_impacts(session, sample: VirusSample, sequence_id: int, aligner):
    pass
    # for (sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type, variant_impacts) in sample.nucleotide_variants_and_effects(aligner):
    #     nuc_variant_db_row = session.query(NucleotideVariant).filter(NucleotideVariant.sequence_id == sequence_id,
    #                                                                  NucleotideVariant.sequence_original == sequence_original,
    #                                                                  NucleotideVariant.sequence_alternative == sequence_alternative,
    #                                                                  NucleotideVariant.start_original == start_original,
    #                                                                  NucleotideVariant.start_alternative == start_alternative,
    #                                                                  NucleotideVariant.variant_length == variant_length,
    #                                                                  NucleotideVariant.variant_type == variant_type).one_or_none()
    #     if not nuc_variant_db_row:
    #         # insert nucleotide variant
    #         nuc_variant_db_row = NucleotideVariant(sequence_id=sequence_id,
    #                                                sequence_original=sequence_original,
    #                                                sequence_alternative=sequence_alternative,
    #                                                start_original=start_original,
    #                                                start_alternative=start_alternative,
    #                                                variant_length=variant_length,
    #                                                variant_type=variant_type)
    #         # insert related impact
    #         for effect, putative_impact, impact_gene_name in variant_impacts:
    #             impact_db_row = VariantImpact(nucleotide_variant_id=nuc_variant_db_row.nucleotide_variant_id,
    #                                           effect=effect,
    #                                           putative_impact=putative_impact,
    #                                           impact_gene_name=impact_gene_name)


def create_nuc_variants_and_impacts(session, sequence_id, args):
    # seq_original = args['sequence_original']
    # seq_alternative = args['sequence_alternative']
    # start_original = args['start_original']
    # start_alternative = args['start_alternative']
    # variant_length = args['variant_length']
    # variant_type = args['variant_type']
    # impacts = args['annotations']
    # nuc_variant_db_row = session.query(NucleotideVariant).filter(NucleotideVariant.sequence_id == sequence_id,
    #                                                              NucleotideVariant.sequence_original == seq_original,
    #                                                              NucleotideVariant.sequence_alternative == seq_alternative,
    #                                                              NucleotideVariant.start_original == start_original,
    #                                                              NucleotideVariant.start_alternative == start_alternative,
    #                                                              NucleotideVariant.variant_length == variant_length,
    #                                                              NucleotideVariant.variant_type == variant_type).one_or_none()
    # if not nuc_variant_db_row:
    #     # insert nucleotide variant
    #     nuc_variant_db_row = NucleotideVariant(sequence_id=sequence_id,
    #                                            sequence_original=seq_original,
    #                                            sequence_alternative=seq_alternative,
    #                                            start_original=start_original,
    #                                            start_alternative=start_alternative,
    #                                            variant_length=variant_length,
    #                                            variant_type=variant_type)
    #     # insert related impact
    #     for effect, putative_impact, impact_gene_name in impacts:
    #         impact_db_row = VariantImpact(nucleotide_variant_id=nuc_variant_db_row.nucleotide_variant_id,
    #                                       effect=effect,
    #                                       putative_impact=putative_impact,
    #                                       impact_gene_name=impact_gene_name)
    pass


def create_epitope(session, epitope: Tuple):
    db_virus_id, host_specie_db_id, host_name, host_iri, protein_name, cell_type, \
    mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links, \
    prediction_process, is_linear = epitope

    epitope = Epitope(virus_id=db_virus_id,
                      host_id=host_specie_db_id,
                      source_host_name=host_name,
                      source_host_iri=host_iri,
                      protein_name=protein_name,
                      cell_type=cell_type,
                      mhc_class=mhc_class,
                      mhc_allele=mhc_allele,
                      response_frequency_positive=response_frequency_positive,
                      epitope_sequence=seq,
                      epi_annotation_start=start,
                      epi_annotation_stop=stop,
                      external_link=ext_links,
                      prediction_process=prediction_process,
                      is_linear=is_linear,
                      assay_type=assay_type
                      )
    return 1


def create_epitope_fragment(session, epi_fragment: Tuple):
    epitope_id, seq, start, stop = epi_fragment
    fragment = EpitopeFragment(epitope_id=epitope_id,
                               epi_fragment_sequence=seq,
                               epi_frag_annotation_start=start,
                               epi_frag_annotation_stop=stop)


def get_virus(session, a_virus) -> Optional[Virus]:
    # return session.query(Virus).filter(Virus.taxon_id == a_virus.taxon_id()).one_or_none()
    v = Virus(virus_id=1)
    return v


def get_specie_id(session, organism_taxon_id:int):
    return 1


def get_reference_sequence_of_virus(session, a_virus: Virus) -> Optional[Sequence]:
    # return session.query(Sequence).filter(
    #     Sequence.virus_id == a_virus.virus_id,
    #     # noqa              # == ignore warning on " == True" for this case
    #     Sequence.is_reference == True
    # ).one_or_none()
    return None


def random_string():
    return ''.join(random.choices(string.ascii_uppercase, k=4))


def sequence_alternative_accession_ids(session, virus_id: int, sources: Optional[List[str]] = None):
    query = session.query(cast(Sequence.alternative_accession_id, sqlalchemy.Integer)).filter(
        Sequence.virus_id == virus_id
    )
    if sources:
        query = query\
            .join(SequencingProject, Sequence.sequencing_project_id == SequencingProject.sequencing_project_id)\
            .filter(SequencingProject.database_source.in_(sources))
    result = query.all()
    return [_[0] for _ in result]


def sequence_primary_accession_ids(session, virus_id: int, sources: Optional[List[str]] = None):
    query = session.query(Sequence.accession_id).filter(
        Sequence.virus_id == virus_id
    )
    if sources:
        query = query\
            .join(SequencingProject, Sequence.sequencing_project_id == SequencingProject.sequencing_project_id)\
            .filter(SequencingProject.database_source.in_(sources))
    result = query.all()
    return [_[0] for _ in result]


def remove_sequence_and_meta(session, primary_sequence_accession_id: Optional[str], alternative_sequence_accession_id: Optional[str]):
    pass


def remove_sequence_and_meta_list(session, primary_sequence_accession_id: Optional[List[str]] = None,
                                  alternative_sequence_accession_id: Optional[List[str]] = None):
    pass


def check_existence_epitopes(session, virus_id):
    one_epitope = session.query(Epitope).filter(Epitope.virus_id == virus_id).first()
    return one_epitope is not None


def update_db_metadata(session, virus_db_id: int, database_source: str):
    pass


def clean_objects_unreachable_from_sequences(session):
    pass


def insert_data_update_pipeline_event(session, event: PipelineEvent):
    pass
