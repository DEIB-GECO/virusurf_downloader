from collections import OrderedDict
from typing import List
# noinspection PyPackageRequirements
from Bio import Entrez
from lxml import etree

from data_sources.virus import VirusSource
from database_tom import ExperimentType, SequencingProject, Virus, HostSample, Sequence, Annotation, NucleotideVariant, \
    VariantImpact, AminoacidVariant
from locations import *
from tqdm import tqdm
from data_sources.virus_sample import VirusSample
from xml_helper import *

# Set of methods that create the corresponding rows in the Virus Conceptual Model database tables


#   #############################    VCM    ############################
def create_or_get_virus(session, a_virus) -> Virus:
    virus = session.query(Virus).filter(Virus.taxon_id == a_virus.taxon_id()).one_or_none()
    if not virus:
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

    host_taxon_name = sample.host_taxon_name()
    host_taxon_id = sample.host_taxon_id()
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


def create_or_get_sequence(session, virus_sample: VirusSample, virus_id: int, experiment: ExperimentType, host_sample: HostSample, sequencing_project: SequencingProject):
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
                                              Sequence.n_percentage == n_percentage,
                                              Sequence.nucleotide_sequence == nucleotide_sequence,
                                              Sequence.strand == strand,
                                              Sequence.length == length,
                                              Sequence.gc_percentage == gc_percentage,
                                              Sequence.lineage == lineage,
                                              Sequence.clade == clade,
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
                            n_percentage=n_percentage,
                            nucleotide_sequence=nucleotide_sequence,
                            strand=strand,
                            length=length,
                            gc_percentage=gc_percentage,
                            lineage=lineage,
                            clade=clade,
                            experiment_type_id=experiment_type_id,
                            sequencing_project_id=sequencing_project_id,
                            virus_id=virus_id,
                            host_sample_id=host_sample_id)
        session.add(sequence)
        session.flush()
    return sequence


def create_annotation_and_aa_variants(session, sample: VirusSample, sequence: Sequence, reference_sample: VirusSample):
    for start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, aa_variants in sample.annotations_and_amino_acid_variants(reference_sample):
        annotation = session.query(Annotation).filter(Annotation.start == start,
                                                      Annotation.stop == stop,
                                                      Annotation.feature_type == feature_type,
                                                      Annotation.gene_name == gene_name,
                                                      Annotation.product == product,
                                                      Annotation.external_reference == db_xref_merged,
                                                      Annotation.sequence_id == sequence.sequence_id,
                                                      Annotation.aminoacid_sequence == amino_acid_sequence).one_or_none()
        if not annotation:
            annotation = Annotation(start=start,
                                    stop=stop,
                                    gene_name=gene_name,
                                    feature_type=feature_type,
                                    product=product,
                                    external_reference=db_xref_merged,
                                    sequence_id=sequence.sequence_id,
                                    aminoacid_sequence=amino_acid_sequence)
            session.add(annotation)
            session.flush()
            if aa_variants:
                for original, alternative, mutpos, mut_len, mut_type in aa_variants:
                    aa_variant = AminoacidVariant(annotation_id=annotation.annotation_id,
                                                  sequence_aa_original=original,
                                                  sequence_aa_alternative=alternative,
                                                  start_aa_original=mutpos,
                                                  variant_aa_length=mut_len,
                                                  variant_aa_type=mut_type)
                    session.add(aa_variant)
            session.flush()


def create_nucleotide_variants_and_impacts(session, sample: VirusSample, sequence_id: int, aligner):
    for (sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type, variant_impacts) in sample.nucleotide_variants_and_effects(aligner):
        nuc_variant_db_row = session.query(NucleotideVariant).filter(NucleotideVariant.sequence_id == sequence_id,
                                                                     NucleotideVariant.sequence_original == sequence_original,
                                                                     NucleotideVariant.sequence_alternative == sequence_alternative,
                                                                     NucleotideVariant.start_original == start_original,
                                                                     NucleotideVariant.start_alternative == start_alternative,
                                                                     NucleotideVariant.variant_length == variant_length,
                                                                     NucleotideVariant.variant_type == variant_type).one_or_none()
        if not nuc_variant_db_row:
            # insert nucleotide variant
            nuc_variant_db_row = NucleotideVariant(sequence_id=sequence_id,
                                                   sequence_original=sequence_original,
                                                   sequence_alternative=sequence_alternative,
                                                   start_original=start_original,
                                                   start_alternative=start_alternative,
                                                   variant_length=variant_length,
                                                   variant_type=variant_type)
            session.add(nuc_variant_db_row)
            session.flush()
            # insert related impact
            for effect, putative_impact, impact_gene_name in variant_impacts:
                impact_db_row = VariantImpact(nucleotide_variant_id=nuc_variant_db_row.nucleotide_variant_id,
                                              effect=effect,
                                              putative_impact=putative_impact,
                                              impact_gene_name=impact_gene_name)
                session.add(impact_db_row)


def get_virus(session, a_virus) -> Optional[Virus]:
    return session.query(Virus).filter(Virus.taxon_id == a_virus.taxon_id()).one_or_none()


def get_reference_sequence_of_virus(session, a_virus: Virus) -> Optional[Sequence]:
    return session.query(Sequence).filter(
        Sequence.virus_id == a_virus.virus_id,
        # noqa  # ignore warning (something == True) for this case
        Sequence.is_reference == True
    ).one_or_none()
