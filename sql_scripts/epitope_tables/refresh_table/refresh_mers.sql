-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 1AB polyprotein
TRUNCATE public.epitope_1335626_1ab_polyprotein;
-- 1335626 can be replaced with the virus taxon id, while 1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_1ab_polyprotein (
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
            AND ann.product = '1AB polyprotein'
            AND epi.protein_name = '1AB polyprotein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_1335626_rna_dependent_rna_polymerase;
-- 1335626 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_rna_dependent_rna_polymerase (
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
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT Hel
TRUNCATE public.epitope_1335626_hel;
-- 1335626 can be replaced with the virus taxon id, while hel can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_hel (
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
            AND ann.product = 'Hel'
            AND epi.protein_name = 'Hel'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT ExoN
TRUNCATE public.epitope_1335626_exon;
-- 1335626 can be replaced with the virus taxon id, while exon can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_exon (
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
            AND ann.product = 'ExoN'
            AND epi.protein_name = 'ExoN'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NendoU
TRUNCATE public.epitope_1335626_nendou;
-- 1335626 can be replaced with the virus taxon id, while nendou can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nendou (
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
            AND ann.product = 'NendoU'
            AND epi.protein_name = 'NendoU'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 2'-O-methyltransferase
TRUNCATE public.epitope_1335626_2_o_methyltransferase;
-- 1335626 can be replaced with the virus taxon id, while 2_o_methyltransferase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_2_o_methyltransferase (
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
            AND ann.product = '2''-O-methyltransferase'
            AND epi.protein_name = '2''-O-methyltransferase'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 1A polyprotein
TRUNCATE public.epitope_1335626_1a_polyprotein;
-- 1335626 can be replaced with the virus taxon id, while 1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_1a_polyprotein (
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
            AND ann.product = '1A polyprotein'
            AND epi.protein_name = '1A polyprotein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp1 protein
TRUNCATE public.epitope_1335626_nsp1_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp1_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp1_protein (
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
            AND ann.product = 'nsp1 protein'
            AND epi.protein_name = 'nsp1 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp2 protein
TRUNCATE public.epitope_1335626_nsp2_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp2_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp2_protein (
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
            AND ann.product = 'nsp2 protein'
            AND epi.protein_name = 'nsp2 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp3 protein
TRUNCATE public.epitope_1335626_nsp3_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp3_protein (
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
            AND ann.product = 'nsp3 protein'
            AND epi.protein_name = 'nsp3 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp4 protein
TRUNCATE public.epitope_1335626_nsp4_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp4_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp4_protein (
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
            AND ann.product = 'nsp4 protein'
            AND epi.protein_name = 'nsp4 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp5 protein
TRUNCATE public.epitope_1335626_nsp5_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp5_protein (
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
            AND ann.product = 'nsp5 protein'
            AND epi.protein_name = 'nsp5 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp6 protein
TRUNCATE public.epitope_1335626_nsp6_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp6_protein (
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
            AND ann.product = 'nsp6 protein'
            AND epi.protein_name = 'nsp6 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp7 protein
TRUNCATE public.epitope_1335626_nsp7_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp7_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp7_protein (
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
            AND ann.product = 'nsp7 protein'
            AND epi.protein_name = 'nsp7 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp8 protein
TRUNCATE public.epitope_1335626_nsp8_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp8_protein (
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
            AND ann.product = 'nsp8 protein'
            AND epi.protein_name = 'nsp8 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp9 protein
TRUNCATE public.epitope_1335626_nsp9_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp9_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp9_protein (
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
            AND ann.product = 'nsp9 protein'
            AND epi.protein_name = 'nsp9 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp10 protein
TRUNCATE public.epitope_1335626_nsp10_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp10_protein (
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
            AND ann.product = 'nsp10 protein'
            AND epi.protein_name = 'nsp10 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp11 protein
TRUNCATE public.epitope_1335626_nsp11_protein;
-- 1335626 can be replaced with the virus taxon id, while nsp11_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nsp11_protein (
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
            AND ann.product = 'nsp11 protein'
            AND epi.protein_name = 'nsp11 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT spike protein
TRUNCATE public.epitope_1335626_spike_protein;
-- 1335626 can be replaced with the virus taxon id, while spike_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_spike_protein (
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
            AND ann.product = 'spike protein'
            AND epi.protein_name = 'spike protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS3 protein
TRUNCATE public.epitope_1335626_ns3_protein;
-- 1335626 can be replaced with the virus taxon id, while ns3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_ns3_protein (
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
            AND ann.product = 'NS3 protein'
            AND epi.protein_name = 'NS3 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS4A protein
TRUNCATE public.epitope_1335626_ns4a_protein;
-- 1335626 can be replaced with the virus taxon id, while ns4a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_ns4a_protein (
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
            AND ann.product = 'NS4A protein'
            AND epi.protein_name = 'NS4A protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS4B protein
TRUNCATE public.epitope_1335626_ns4b_protein;
-- 1335626 can be replaced with the virus taxon id, while ns4b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_ns4b_protein (
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
            AND ann.product = 'NS4B protein'
            AND epi.protein_name = 'NS4B protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS5 protein
TRUNCATE public.epitope_1335626_ns5_protein;
-- 1335626 can be replaced with the virus taxon id, while ns5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_ns5_protein (
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
            AND ann.product = 'NS5 protein'
            AND epi.protein_name = 'NS5 protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT envelope protein
TRUNCATE public.epitope_1335626_envelope_protein;
-- 1335626 can be replaced with the virus taxon id, while envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_envelope_protein (
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
            AND ann.product = 'envelope protein'
            AND epi.protein_name = 'envelope protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT membrane protein
TRUNCATE public.epitope_1335626_membrane_protein;
-- 1335626 can be replaced with the virus taxon id, while membrane_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_membrane_protein (
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
            AND ann.product = 'membrane protein'
            AND epi.protein_name = 'membrane protein'
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nucleocapsid protein
TRUNCATE public.epitope_1335626_nucleocapsid_protein;
-- 1335626 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_nucleocapsid_protein (
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
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT ORF8b protein
TRUNCATE public.epitope_1335626_orf8b_protein;
-- 1335626 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_orf8b_protein (
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
            AND vir.taxon_id = 1335626)
ORDER BY epi.iedb_epitope_id;

