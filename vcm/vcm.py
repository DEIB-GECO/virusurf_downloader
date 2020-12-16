from typing import List, Tuple
from db_config.database_tom import ExperimentType, SequencingProject, Virus, HostSample, Sequence, Annotation, NucleotideVariant, \
    VariantImpact, AminoAcidVariant, Epitope, EpitopeFragment, HostSpecie, DBMeta
from locations import *
from data_sources.virus_sample import VirusSample
from xml_helper import *
from datetime import datetime

"""
Set of methods that create the corresponding rows in the Virus Conceptual Model database tables
"""

cache_host_specie = dict()
cache_host_sample = dict()
cache_experiment_type = dict()
cache_sequencing_project = dict()
cache_virus = dict()
epitope_id_mappings = dict()

cache_typos_corrections = dict()


#   #############################    VCM    ############################
def create_or_get_virus(session, a_virus):
    global cache_virus

    virus_key = (a_virus.taxon_id(), a_virus.taxon_name(), a_virus.family(), a_virus.sub_family(), a_virus.genus(), a_virus.species(), a_virus.equivalent_names(), a_virus.molecule_type(), a_virus.is_single_stranded(), a_virus.is_positive_stranded())
    virus_id = cache_virus.get(virus_key)

    if not virus_id:
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
        virus_id = virus.virus_id
        cache_virus[virus_key] = virus_id
    return virus_id


def create_or_get_experiment(session, sample: VirusSample):
    global cache_experiment_type

    sequencing_technology = sample.sequencing_technology()
    assembly_method = sample.assembly_method()
    coverage = sample.coverage()

    exp_key = (sequencing_technology, assembly_method, coverage)
    experiment_id = cache_experiment_type.get(exp_key)

    if not experiment_id:
        experiment = session.query(ExperimentType).filter(
            ExperimentType.sequencing_technology == sequencing_technology,
            ExperimentType.assembly_method == assembly_method,
            ExperimentType.coverage == coverage).one_or_none()
        if not experiment:
            experiment = ExperimentType(
                sequencing_technology=sequencing_technology,
                assembly_method=assembly_method,
                coverage=coverage)
            session.add(experiment)
            session.flush()
        experiment_id = experiment.experiment_type_id
        cache_experiment_type[exp_key] = experiment_id
    return experiment_id


def create_or_get_sequencing_project(session, sample: VirusSample):
    global cache_sequencing_project

    submission_date = sample.submission_date()
    sequencing_lab = sample.sequencing_lab()
    bioproject_id = sample.bioproject_id()
    database_source = sample.database_source()

    project_key = (submission_date, sequencing_lab, bioproject_id, database_source)
    sequencing_project_id = cache_sequencing_project.get(project_key)

    if not sequencing_project_id:
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
        sequencing_project_id = sequencing_project.sequencing_project_id
        cache_sequencing_project[project_key] = sequencing_project_id
    return sequencing_project_id


def create_or_get_host_specie(session, sample: VirusSample) -> int:
    global cache_host_specie

    host_taxon_name = sample.host_taxon_name()
    if host_taxon_name:
        host_taxon_name = host_taxon_name.lower()
    host_taxon_id = sample.host_taxon_id()

    host_specie_key = host_taxon_name
    host_specie_id = cache_host_specie.get(host_specie_key)
    if not host_specie_id:
        host_specie = session.query(HostSpecie).filter(HostSpecie.host_taxon_name == host_taxon_name).one_or_none()
        if not host_specie:
            host_specie = HostSpecie(host_taxon_id=host_taxon_id,
                                     host_taxon_name=host_taxon_name)
            session.add(host_specie)
            session.flush()
        host_specie_id = host_specie.host_id
        cache_host_specie[host_specie_key] = host_specie_id
    return host_specie_id


def create_or_get_host_specie_alt(session, organism_name: str, organism_ncbi_id:int):
    global cache_host_specie

    host_specie_id = cache_host_specie.get(organism_name)
    if not host_specie_id:
        host_specie = session.query(HostSpecie).filter(HostSpecie.host_taxon_name == organism_name).one_or_none()
        if not host_specie:
            host_specie = HostSpecie(host_taxon_id=organism_ncbi_id,
                                     host_taxon_name=organism_name)
            session.add(host_specie)
            session.flush()
        host_specie_id = host_specie.host_id
        cache_host_specie[organism_name] = host_specie_id
    return host_specie_id


def create_or_get_host_sample(session, sample: VirusSample, host_specie_id: int) -> int:
    global cache_host_sample

    originating_lab = sample.originating_lab()
    collection_date = sample.collection_date()
    isolation_source = sample.isolation_source()

    gender = sample.gender()
    age = sample.age()

    country, region, geo_group = sample.country__region__geo_group()

    host_sample_key = (host_specie_id, gender, age, originating_lab, collection_date, isolation_source, country, region, geo_group)

    if host_sample_key in cache_host_sample:
        host_sample_id = cache_host_sample[host_sample_key]
    else:
        host_sample = session.query(HostSample).filter(HostSample.host_id == host_specie_id,
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
            host_sample = HostSample(host_id=host_specie_id,
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
        host_sample_id = host_sample.host_sample_id
        cache_host_sample[host_sample_key] = host_sample_id
    return host_sample_id


def get_sequence(session, virus_sample: VirusSample, virus_id) -> Sequence:
    # data from sample
    accession_id = virus_sample.primary_accession_number()
    alternative_accession_id = virus_sample.alternative_accession_number()

    sequence = session.query(Sequence).filter(Sequence.alternative_accession_id == alternative_accession_id,
                                              Sequence.accession_id == accession_id,
                                              Sequence.virus_id == virus_id
                                              ).one_or_none()
    return sequence


def create_and_get_sequence(session, virus_sample: VirusSample, virus_id, experiment_id, host_sample_id, sequencing_project_id):
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
                        experiment_type_id=experiment_id,
                        sequencing_project_id=sequencing_project_id,
                        virus_id=virus_id,
                        host_sample_id=host_sample_id)
    session.add(sequence)
    session.flush()
    return sequence


def create_annotation_and_aa_variants(session, sample: VirusSample, sequence: Sequence, reference_sample: VirusSample):
    for start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, aa_variants in sample.annotations_and_amino_acid_variants(reference_sample):
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
                aa_variant = AminoAcidVariant(annotation_id=annotation.annotation_id,
                                              sequence_aa_original=original,
                                              sequence_aa_alternative=alternative,
                                              start_aa_original=mutpos,
                                              variant_aa_length=mut_len,
                                              variant_aa_type=mut_type)
                session.add(aa_variant)
            session.flush()


def create_annotation_and_amino_acid_variants(session, sequence_id, *args):
    # print(args)
    gene_name, product, protein_id, feature_type, start, stop, nuc_seq, amino_acid_seq, aa_variants = args
    annotation = Annotation(start=start,
                            stop=stop,
                            gene_name=gene_name,
                            feature_type=feature_type,
                            product=product,
                            sequence_id=sequence_id,
                            external_reference=protein_id,
                            annotation_nucleotide_sequence=nuc_seq,
                            aminoacid_sequence=amino_acid_seq)
    session.add(annotation)
    session.flush()
    if aa_variants:
        for gen_name, protein, protein_id, start_pos, sequence_aa_original, sequence_aa_alternative, variant_aa_type in aa_variants:
            aa_variant = AminoAcidVariant(annotation_id=annotation.annotation_id,
                                          sequence_aa_original=sequence_aa_original,
                                          sequence_aa_alternative=sequence_aa_alternative,
                                          start_aa_original=int(start_pos) if start_pos else None,
                                          variant_aa_length=max(len(sequence_aa_original), len(sequence_aa_alternative)),
                                          variant_aa_type=variant_aa_type)
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


def create_nuc_variants_and_impacts(session, sequence_id, args):
    seq_original = args['sequence_original']
    seq_alternative = args['sequence_alternative']
    start_original = args['start_original']
    start_alternative = args['start_alternative']
    variant_length = args['variant_length']
    variant_type = args['variant_type']
    impacts = args['annotations']
    # insert nucleotide variant
    nuc_variant_db_row = NucleotideVariant(sequence_id=sequence_id,
                                           sequence_original=seq_original,
                                           sequence_alternative=seq_alternative,
                                           start_original=start_original,
                                           start_alternative=start_alternative,
                                           variant_length=variant_length,
                                           variant_type=variant_type)
    session.add(nuc_variant_db_row)
    session.flush()
    # insert related impact
    for effect, putative_impact, impact_gene_name in impacts:
        impact_db_row = VariantImpact(nucleotide_variant_id=nuc_variant_db_row.nucleotide_variant_id,
                                      effect=effect,
                                      putative_impact=putative_impact,
                                      impact_gene_name=impact_gene_name)
        session.add(impact_db_row)
    session.flush()


def create_epitopes(session, epitopes: List[Tuple], virus_id, _host_specie_id):
    global epitope_id_mappings
    epitope_id_mappings = dict()
    for elem in epitopes:
        pseudo_epi_id, _virus_taxon_id, _, protein_ncbi_id, _type, hla_restriction, response_frequency, epitope_sequence, \
        epi_annotation_start, epi_annotation_stop, is_imported, external_link, perdition_process, is_linear = elem
        epitope = Epitope(virus_id=virus_id,
                          host_id=_host_specie_id,
                          protein_ncbi_id=protein_ncbi_id,
                          epitope_type=_type,
                          hla_restriction=hla_restriction,
                          response_frequency=response_frequency,
                          epitope_sequence=epitope_sequence,
                          epi_annotation_start=epi_annotation_start,
                          epi_annotation_stop=epi_annotation_stop,
                          is_imported=is_imported,
                          external_link=external_link,
                          prediction_process=perdition_process,
                          is_linear=is_linear
                          )
        session.add(epitope)
        session.flush()
        epitope_id_mappings[pseudo_epi_id] = epitope.epitope_id


def create_epitope(session, epitope: Tuple):
    db_virus_id, host_specie_db_id, host_name, host_iri, protein_ncbi_id, cell_type, \
    mhc_class, mhc_allele, response_frequency_positive, assay_type, seq, start, stop, ext_links, \
    prediction_process, is_linear = epitope

    epitope = Epitope(virus_id=db_virus_id,
                      host_id=host_specie_db_id,
                      source_host_name=host_name,
                      source_host_iri=host_iri,
                      protein_ncbi_id=protein_ncbi_id,
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
    session.add(epitope)
    session.flush()
    return epitope.epitope_id


def create_epitope_fragment(session, epi_fragment: Tuple):
    epitope_id, seq, start, stop = epi_fragment
    fragment = EpitopeFragment(epitope_id=epitope_id,
                               epi_fragment_sequence=seq,
                               epi_frag_annotation_start=start,
                               epi_frag_annotation_stop=stop)
    session.add(fragment)


def create_epitopes_fragments(session, epi_fragments: List[Tuple]):
    global epitope_id_mappings
    if not epitope_id_mappings or len(epitope_id_mappings.keys()) == 0:
        raise ValueError('You should create epitopes rows before calling this method')
    for elem in epi_fragments:
        _, epitope_pseudo_id, epi_fragment_sequence, epi_fragment_annot_start, epi_fragment_annot_stop = elem
        try:
            real_id = epitope_id_mappings[epitope_pseudo_id]
        except KeyError as e:
            logger.error(f'the epitope fragment ID {epitope_pseudo_id} does not appear in the epitope IDs. This epitope fragment'
                         f' will be not inserted into the DB.')
            continue
        fragment = EpitopeFragment(epitope_id=real_id,
                                   epi_fragment_sequence=epi_fragment_sequence,
                                   epi_frag_annotation_start=epi_fragment_annot_start,
                                   epi_frag_annotation_stop=epi_fragment_annot_stop)
        session.add(fragment)


def get_virus(session, a_virus) -> Optional[Virus]:
    return session.query(Virus).filter(Virus.taxon_id == a_virus.taxon_id()).one_or_none()


def get_specie_id(session, organism_taxon_id:int):
    global cache_host_specie

    host_specie_id = cache_host_specie.get(organism_taxon_id)
    if not host_specie_id:
        host_specie = session.query(HostSpecie).filter(HostSpecie.host_taxon_id == organism_taxon_id).one_or_none()
        if host_specie:
            host_specie_id = host_specie.host_id
            cache_host_specie[organism_taxon_id] = host_specie_id
    return host_specie_id


def get_reference_sequence_of_virus(session, virus_id) -> Optional[Sequence]:
    return session.query(Sequence).filter(
        Sequence.virus_id == virus_id,
        # noqa              # == ignore warning on " == True" for this case
        Sequence.is_reference == True
    ).one_or_none()


def sequence_alternative_accession_ids(session, virus_id: int, sources: Optional[List[str]] = None):
    query = session.query(Sequence.alternative_accession_id).filter(
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


def remove_sequence_and_meta(session, primary_sequence_accession_id: Optional[str]=None, alternative_sequence_accession_id: Optional[str]=None):
    # get metadata ids
    query = session.query(Sequence.sequence_id, Sequence.experiment_type_id, Sequence.sequencing_project_id,
                          Sequence.host_sample_id)
    if primary_sequence_accession_id:
        query = query\
            .filter(Sequence.accession_id == primary_sequence_accession_id)
    elif alternative_sequence_accession_id:
        query = query\
            .filter(Sequence.alternative_accession_id == alternative_sequence_accession_id)
    else:
        raise ValueError('one between primary_sequence_accession_id and alternative_sequence_accession_id arguments must '
                         'be specified')
    sequence_id, experiment_id, sequence_project_id, host_sample_id = query.one()

    # delete aa variants and annotations
    annotation_ids = session.query(Annotation.annotation_id).filter(Annotation.sequence_id == sequence_id).all()
    annotation_ids = [_[0] for _ in annotation_ids]
    if annotation_ids:
        session.query(AminoAcidVariant).filter(AminoAcidVariant.annotation_id.in_(annotation_ids)).delete(synchronize_session=False)
    session.query(Annotation).filter(Annotation.sequence_id == sequence_id).delete(synchronize_session=False)

    # delete impacts and nuc variants
    nuc_variant_ids = session.query(NucleotideVariant.nucleotide_variant_id).filter(NucleotideVariant.sequence_id == sequence_id).all()
    nuc_variant_ids = [_[0] for _ in nuc_variant_ids]
    if nuc_variant_ids:
        session.query(VariantImpact).filter(VariantImpact.nucleotide_variant_id.in_(nuc_variant_ids)).delete(synchronize_session=False)
    session.query(NucleotideVariant).filter(NucleotideVariant.sequence_id == sequence_id).delete(synchronize_session=False)

    # delete sequence
    session.query(Sequence).filter(Sequence.sequence_id == sequence_id).delete(synchronize_session=False)

    # delete related meta
    # (host sample)
    sequences_with_same_host = session.query(Sequence.sequence_id)\
        .filter(Sequence.host_sample_id == host_sample_id,
                Sequence.sequence_id != sequence_id).first()
    if not sequences_with_same_host:
        session.query(HostSample).filter(HostSample.host_sample_id == host_sample_id).delete(synchronize_session=False)

    # (experiment)
    sequences_with_same_experiment = session.query(Sequence.sequence_id) \
        .filter(Sequence.experiment_type_id == experiment_id,
                Sequence.sequence_id != sequence_id).first()
    if not sequences_with_same_experiment:
        session.query(ExperimentType).filter(ExperimentType.experiment_type_id == experiment_id).delete(synchronize_session=False)

    # (seq project)
    sequences_with_same_sequencing_project = session.query(Sequence.sequence_id) \
        .filter(Sequence.sequencing_project_id == sequence_project_id,
                Sequence.sequence_id != sequence_id).first()
    if not sequences_with_same_sequencing_project:
        session.query(SequencingProject).filter(SequencingProject.sequencing_project_id == sequence_project_id).delete(synchronize_session=False)


def remove_sequence_and_meta_list(session, primary_sequence_accession_id: Optional[List[str]]=None, alternative_sequence_accession_id: Optional[List[str]]=None):
    # get sequence_id of all
    query = session.query(Sequence.sequence_id)
    if primary_sequence_accession_id:
        query = query\
            .filter(Sequence.accession_id.in_(primary_sequence_accession_id))
    elif alternative_sequence_accession_id:
        query = query\
            .filter(Sequence.alternative_accession_id.in_(alternative_sequence_accession_id))
    else:
        raise ValueError('one between primary_sequence_accession_id and alternative_sequence_accession_id arguments must '
                         'be specified')
    sequence_ids = query.all()
    sequence_ids = [_[0] for _ in sequence_ids]

    # delete aa variants and annotations
    annotation_ids = session.query(Annotation.annotation_id).filter(Annotation.sequence_id.in_(sequence_ids)).all()
    annotation_ids = [_[0] for _ in annotation_ids]
    if annotation_ids:
        session.query(AminoAcidVariant).filter(AminoAcidVariant.annotation_id.in_(annotation_ids)).delete(synchronize_session=False)
    session.query(Annotation).filter(Annotation.sequence_id.in_(sequence_ids)).delete(synchronize_session=False)

    # delete impacts and nuc variants
    nuc_variant_ids = session.query(NucleotideVariant.nucleotide_variant_id).filter(NucleotideVariant.sequence_id.in_(sequence_ids)).all()
    nuc_variant_ids = [_[0] for _ in nuc_variant_ids]
    if nuc_variant_ids:
        session.query(VariantImpact).filter(VariantImpact.nucleotide_variant_id.in_(nuc_variant_ids)).delete(synchronize_session=False)
    session.query(NucleotideVariant).filter(NucleotideVariant.sequence_id.in_(sequence_ids)).delete(synchronize_session=False)

    # delete sequence
    session.query(Sequence).filter(Sequence.sequence_id.in_(sequence_ids)).delete(synchronize_session=False)

    # delete unused meta
    # (host sample)
    session.query(HostSample)\
        .filter(HostSample.host_sample_id.notin_(session.query(Sequence.host_sample_id))).delete(synchronize_session=False)

    # (experiment)
    session.query(ExperimentType) \
        .filter(ExperimentType.experiment_type_id.notin_(session.query(Sequence.experiment_type_id))).delete(synchronize_session=False)

    # (seq project)
    session.query(SequencingProject) \
        .filter(SequencingProject.sequencing_project_id.notin_(session.query(Sequence.sequencing_project_id))).delete(synchronize_session=False)

    # (host specie)
    session.query(HostSpecie) \
        .filter(
        (HostSpecie.host_id.notin_(session.query(HostSample.host_id))) &
        (HostSpecie.host_id.notin_(session.query(Epitope.host_id)))).delete(synchronize_session=False)


def check_existence_epitopes(session, virus_id):
    one_epitope = session.query(Epitope).filter(Epitope.virus_id == virus_id).first()
    return one_epitope is not None


def update_db_metadata(session, virus_db_id: int, database_source: str):
    current_date = datetime.date(datetime.now())
    last_update = session.query(DBMeta).filter(DBMeta.virus_id == virus_db_id, DBMeta.source == database_source).one_or_none()
    if not last_update:
        meta = DBMeta(virus_id=virus_db_id, date_of_import=current_date, source=database_source)
        session.add(meta)
    else:
        last_update.date_of_import = current_date
        session.merge(last_update)
