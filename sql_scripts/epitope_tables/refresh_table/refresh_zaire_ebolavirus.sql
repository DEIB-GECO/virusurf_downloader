-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_186538_nucleoprotein;
-- 186538 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_nucleoprotein (
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
            AND ann.product = 'nucleoprotein'
            AND epi.protein_name = 'nucleoprotein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_186538_polymerase_complex_protein;
-- 186538 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_polymerase_complex_protein (
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
            AND ann.product = 'polymerase complex protein'
            AND epi.protein_name = 'polymerase complex protein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_186538_matrix_protein;
-- 186538 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_matrix_protein (
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
            AND ann.product = 'matrix protein'
            AND epi.protein_name = 'matrix protein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_186538_spike_glycoprotein;
-- 186538 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_spike_glycoprotein (
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
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_186538_small_secreted_glycoprotein;
-- 186538 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_small_secreted_glycoprotein (
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
            AND ann.product = 'small secreted glycoprotein'
            AND epi.protein_name = 'small secreted glycoprotein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_186538_second_secreted_glycoprotein;
-- 186538 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_second_secreted_glycoprotein (
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
            AND ann.product = 'second secreted glycoprotein'
            AND epi.protein_name = 'second secreted glycoprotein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_186538_minor_nucleoprotein;
-- 186538 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_minor_nucleoprotein (
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
            AND ann.product = 'minor nucleoprotein'
            AND epi.protein_name = 'minor nucleoprotein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_186538_membrane_associated_protein;
-- 186538 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_membrane_associated_protein (
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
            AND ann.product = 'membrane-associated protein'
            AND epi.protein_name = 'membrane-associated protein'
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_186538_rna_dependent_rna_polymerase;
-- 186538 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_rna_dependent_rna_polymerase (
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
            AND vir.taxon_id = 186538)
ORDER BY epi.iedb_epitope_id;

