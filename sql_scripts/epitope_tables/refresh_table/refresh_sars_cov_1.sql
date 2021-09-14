-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1ab polyprotein
TRUNCATE public.epitope_694009_orf1ab_polyprotein;
-- 694009 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf1ab_polyprotein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF1ab polyprotein'
            AND epi.protein_name = 'ORF1ab polyprotein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_694009_rna_dependent_rna_polymerase;
-- 694009 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_rna_dependent_rna_polymerase (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'RNA-dependent RNA polymerase'
            AND epi.protein_name = 'RNA-dependent RNA polymerase'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT helicase/NTPase
TRUNCATE public.epitope_694009_helicase_ntpase;
-- 694009 can be replaced with the virus taxon id, while helicase_ntpase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_helicase_ntpase (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'helicase/NTPase'
            AND epi.protein_name = 'helicase/NTPase'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3' to 5' exonuclease
TRUNCATE public.epitope_694009_3_to_5_exonuclease;
-- 694009 can be replaced with the virus taxon id, while 3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_3_to_5_exonuclease (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = '3'' to 5'' exonuclease'
            AND epi.protein_name = '3'' to 5'' exonuclease'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT endoribonuclease
TRUNCATE public.epitope_694009_endoribonuclease;
-- 694009 can be replaced with the virus taxon id, while endoribonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_endoribonuclease (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'endoribonuclease'
            AND epi.protein_name = 'endoribonuclease'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 2'-O-MTase
TRUNCATE public.epitope_694009_2_o_mtase;
-- 694009 can be replaced with the virus taxon id, while 2_o_mtase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_2_o_mtase (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = '2''-O-MTase'
            AND epi.protein_name = '2''-O-MTase'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1a polyprotein
TRUNCATE public.epitope_694009_orf1a_polyprotein;
-- 694009 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf1a_polyprotein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF1a polyprotein'
            AND epi.protein_name = 'ORF1a polyprotein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp1
TRUNCATE public.epitope_694009_nsp1;
-- 694009 can be replaced with the virus taxon id, while nsp1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp1 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp1'
            AND epi.protein_name = 'nsp1'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp2
TRUNCATE public.epitope_694009_nsp2;
-- 694009 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp2 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp2'
            AND epi.protein_name = 'nsp2'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp3
TRUNCATE public.epitope_694009_nsp3;
-- 694009 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp3 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp3'
            AND epi.protein_name = 'nsp3'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp4
TRUNCATE public.epitope_694009_nsp4;
-- 694009 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp4 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp4'
            AND epi.protein_name = 'nsp4'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3C-like protease
TRUNCATE public.epitope_694009_3c_like_protease;
-- 694009 can be replaced with the virus taxon id, while 3c_like_protease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_3c_like_protease (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = '3C-like protease'
            AND epi.protein_name = '3C-like protease'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp6
TRUNCATE public.epitope_694009_nsp6;
-- 694009 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp6 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp6'
            AND epi.protein_name = 'nsp6'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp7
TRUNCATE public.epitope_694009_nsp7;
-- 694009 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp7 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp7'
            AND epi.protein_name = 'nsp7'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp8
TRUNCATE public.epitope_694009_nsp8;
-- 694009 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp8 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp8'
            AND epi.protein_name = 'nsp8'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp9
TRUNCATE public.epitope_694009_nsp9;
-- 694009 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp9 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp9'
            AND epi.protein_name = 'nsp9'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp10
TRUNCATE public.epitope_694009_nsp10;
-- 694009 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nsp10 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nsp10'
            AND epi.protein_name = 'nsp10'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ndp11
TRUNCATE public.epitope_694009_ndp11;
-- 694009 can be replaced with the virus taxon id, while ndp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_ndp11 (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ndp11'
            AND epi.protein_name = 'ndp11'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT spike glycoprotein
TRUNCATE public.epitope_694009_spike_glycoprotein;
-- 694009 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_spike_glycoprotein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'spike glycoprotein'
            AND epi.protein_name = 'spike glycoprotein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3a protein
TRUNCATE public.epitope_694009_orf3a_protein;
-- 694009 can be replaced with the virus taxon id, while orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf3a_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF3a protein'
            AND epi.protein_name = 'ORF3a protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3b protein
TRUNCATE public.epitope_694009_orf3b_protein;
-- 694009 can be replaced with the virus taxon id, while orf3b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf3b_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF3b protein'
            AND epi.protein_name = 'ORF3b protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT small envelope protein
TRUNCATE public.epitope_694009_small_envelope_protein;
-- 694009 can be replaced with the virus taxon id, while small_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_small_envelope_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'small envelope protein'
            AND epi.protein_name = 'small envelope protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT membrane glycoprotein M
TRUNCATE public.epitope_694009_membrane_glycoprotein_m;
-- 694009 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_membrane_glycoprotein_m (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'membrane glycoprotein M'
            AND epi.protein_name = 'membrane glycoprotein M'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF6 protein
TRUNCATE public.epitope_694009_orf6_protein;
-- 694009 can be replaced with the virus taxon id, while orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf6_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF6 protein'
            AND epi.protein_name = 'ORF6 protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7a protein
TRUNCATE public.epitope_694009_orf7a_protein;
-- 694009 can be replaced with the virus taxon id, while orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf7a_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF7a protein'
            AND epi.protein_name = 'ORF7a protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7b protein
TRUNCATE public.epitope_694009_orf7b_protein;
-- 694009 can be replaced with the virus taxon id, while orf7b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf7b_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF7b protein'
            AND epi.protein_name = 'ORF7b protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8a protein
TRUNCATE public.epitope_694009_orf8a_protein;
-- 694009 can be replaced with the virus taxon id, while orf8a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf8a_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF8a protein'
            AND epi.protein_name = 'ORF8a protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8b protein
TRUNCATE public.epitope_694009_orf8b_protein;
-- 694009 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf8b_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF8b protein'
            AND epi.protein_name = 'ORF8b protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nucleocapsid protein
TRUNCATE public.epitope_694009_nucleocapsid_protein;
-- 694009 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_nucleocapsid_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'nucleocapsid protein'
            AND epi.protein_name = 'nucleocapsid protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9b protein
TRUNCATE public.epitope_694009_orf9b_protein;
-- 694009 can be replaced with the virus taxon id, while orf9b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf9b_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF9b protein'
            AND epi.protein_name = 'ORF9b protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9a protein
TRUNCATE public.epitope_694009_orf9a_protein;
-- 694009 can be replaced with the virus taxon id, while orf9a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf9a_protein (
    iedb_epitope_id,
    epitope_iri,
    cell_type,
    mhc_class,
    mhc_allele,
    response_frequency_pos,
    epi_annotation_start,
    epi_annotation_stop,
    is_linear,
    assay_type,
    epi_fragment_sequence,
    epi_frag_annotation_start,
    epi_frag_annotation_stop,
    taxon_id,
    taxon_name,
    host_taxon_id,
    host_taxon_name,
    sequence_id,
    product,
    aminoacid_variant_id,
    start_aa_original,
    sequence_aa_original,
    sequence_aa_alternative,
    variant_aa_length,
    variant_aa_type
) SELECT DISTINCT
    epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
FROM (((((((epitope epi
    JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
    JOIN virus vir ON ((epi.virus_id = vir.virus_id)))
    JOIN host_specie hspec ON ((epi.host_id = hspec.host_id)))
    JOIN host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
    JOIN sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
    JOIN annotation ann ON ((seq.sequence_id = ann.sequence_id)))
    JOIN aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
WHERE (epi.protein_name::text = ann.product::text
            AND amin.start_aa_original <= epif.epi_frag_annotation_stop
            AND amin.start_aa_original >= epif.epi_frag_annotation_start
            AND ann.product = 'ORF9a protein'
            AND epi.protein_name = 'ORF9a protein'
            AND vir.taxon_id = 694009)
ORDER BY epi.iedb_epitope_id;

