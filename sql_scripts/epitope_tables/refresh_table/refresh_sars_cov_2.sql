-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1ab polyprotein
TRUNCATE public.epitope_2697049_orf1ab_polyprotein;
-- 2697049 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_orf1ab_polyprotein (
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
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP12 (RNA-dependent RNA polymerase)
TRUNCATE public.epitope_2697049_nsp12_rna_dependent_rna_poly;
-- 2697049 can be replaced with the virus taxon id, while nsp12_rna_dependent_rna_poly can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp12_rna_dependent_rna_poly (
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
            AND ann.product = 'NSP12 (RNA-dependent RNA polymerase)'
            AND epi.protein_name = 'NSP12 (RNA-dependent RNA polymerase)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP13 (helicase)
TRUNCATE public.epitope_2697049_nsp13_helicase;
-- 2697049 can be replaced with the virus taxon id, while nsp13_helicase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp13_helicase (
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
            AND ann.product = 'NSP13 (helicase)'
            AND epi.protein_name = 'NSP13 (helicase)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP14 (3'-to-5' exonuclease)
TRUNCATE public.epitope_2697049_nsp14_3_to_5_exonuclease;
-- 2697049 can be replaced with the virus taxon id, while nsp14_3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp14_3_to_5_exonuclease (
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
            AND ann.product = 'NSP14 (3''-to-5'' exonuclease)'
            AND epi.protein_name = 'NSP14 (3''-to-5'' exonuclease)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP15 (endoRNAse)
TRUNCATE public.epitope_2697049_nsp15_endornase;
-- 2697049 can be replaced with the virus taxon id, while nsp15_endornase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp15_endornase (
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
            AND ann.product = 'NSP15 (endoRNAse)'
            AND epi.protein_name = 'NSP15 (endoRNAse)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP16 (2'-O-ribose methyltransferase)
TRUNCATE public.epitope_2697049_nsp16_2_o_ribose_methyltrans;
-- 2697049 can be replaced with the virus taxon id, while nsp16_2_o_ribose_methyltrans can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp16_2_o_ribose_methyltrans (
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
            AND ann.product = 'NSP16 (2''-O-ribose methyltransferase)'
            AND epi.protein_name = 'NSP16 (2''-O-ribose methyltransferase)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1a polyprotein
TRUNCATE public.epitope_2697049_orf1a_polyprotein;
-- 2697049 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_orf1a_polyprotein (
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
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP1 (leader protein)
TRUNCATE public.epitope_2697049_nsp1_leader_protein;
-- 2697049 can be replaced with the virus taxon id, while nsp1_leader_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp1_leader_protein (
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
            AND ann.product = 'NSP1 (leader protein)'
            AND epi.protein_name = 'NSP1 (leader protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP2
TRUNCATE public.epitope_2697049_nsp2;
-- 2697049 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp2 (
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
            AND ann.product = 'NSP2'
            AND epi.protein_name = 'NSP2'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP3
TRUNCATE public.epitope_2697049_nsp3;
-- 2697049 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp3 (
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
            AND ann.product = 'NSP3'
            AND epi.protein_name = 'NSP3'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP4
TRUNCATE public.epitope_2697049_nsp4;
-- 2697049 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp4 (
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
            AND ann.product = 'NSP4'
            AND epi.protein_name = 'NSP4'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP5 (3C-like proteinase)
TRUNCATE public.epitope_2697049_nsp5_3c_like_proteinase;
-- 2697049 can be replaced with the virus taxon id, while nsp5_3c_like_proteinase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp5_3c_like_proteinase (
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
            AND ann.product = 'NSP5 (3C-like proteinase)'
            AND epi.protein_name = 'NSP5 (3C-like proteinase)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP6
TRUNCATE public.epitope_2697049_nsp6;
-- 2697049 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp6 (
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
            AND ann.product = 'NSP6'
            AND epi.protein_name = 'NSP6'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP7
TRUNCATE public.epitope_2697049_nsp7;
-- 2697049 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp7 (
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
            AND ann.product = 'NSP7'
            AND epi.protein_name = 'NSP7'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP8
TRUNCATE public.epitope_2697049_nsp8;
-- 2697049 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp8 (
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
            AND ann.product = 'NSP8'
            AND epi.protein_name = 'NSP8'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP9
TRUNCATE public.epitope_2697049_nsp9;
-- 2697049 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp9 (
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
            AND ann.product = 'NSP9'
            AND epi.protein_name = 'NSP9'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP10
TRUNCATE public.epitope_2697049_nsp10;
-- 2697049 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp10 (
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
            AND ann.product = 'NSP10'
            AND epi.protein_name = 'NSP10'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP11
TRUNCATE public.epitope_2697049_nsp11;
-- 2697049 can be replaced with the virus taxon id, while nsp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_nsp11 (
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
            AND ann.product = 'NSP11'
            AND epi.protein_name = 'NSP11'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT Spike (surface glycoprotein)
TRUNCATE public.epitope_2697049_spike_surface_glycoprotein;
-- 2697049 can be replaced with the virus taxon id, while spike_surface_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_spike_surface_glycoprotein (
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
            AND ann.product = 'Spike (surface glycoprotein)'
            AND epi.protein_name = 'Spike (surface glycoprotein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS3 (ORF3a protein)
TRUNCATE public.epitope_2697049_ns3_orf3a_protein;
-- 2697049 can be replaced with the virus taxon id, while ns3_orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_ns3_orf3a_protein (
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
            AND ann.product = 'NS3 (ORF3a protein)'
            AND epi.protein_name = 'NS3 (ORF3a protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT E (envelope protein)
TRUNCATE public.epitope_2697049_e_envelope_protein;
-- 2697049 can be replaced with the virus taxon id, while e_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_e_envelope_protein (
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
            AND ann.product = 'E (envelope protein)'
            AND epi.protein_name = 'E (envelope protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT M (membrane glycoprotein)
TRUNCATE public.epitope_2697049_m_membrane_glycoprotein;
-- 2697049 can be replaced with the virus taxon id, while m_membrane_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_m_membrane_glycoprotein (
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
            AND ann.product = 'M (membrane glycoprotein)'
            AND epi.protein_name = 'M (membrane glycoprotein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS6 (ORF6 protein)
TRUNCATE public.epitope_2697049_ns6_orf6_protein;
-- 2697049 can be replaced with the virus taxon id, while ns6_orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_ns6_orf6_protein (
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
            AND ann.product = 'NS6 (ORF6 protein)'
            AND epi.protein_name = 'NS6 (ORF6 protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7a (ORF7a protein)
TRUNCATE public.epitope_2697049_ns7a_orf7a_protein;
-- 2697049 can be replaced with the virus taxon id, while ns7a_orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_ns7a_orf7a_protein (
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
            AND ann.product = 'NS7a (ORF7a protein)'
            AND epi.protein_name = 'NS7a (ORF7a protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7b (ORF7b)
TRUNCATE public.epitope_2697049_ns7b_orf7b;
-- 2697049 can be replaced with the virus taxon id, while ns7b_orf7b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_ns7b_orf7b (
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
            AND ann.product = 'NS7b (ORF7b)'
            AND epi.protein_name = 'NS7b (ORF7b)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS8 (ORF8 protein)
TRUNCATE public.epitope_2697049_ns8_orf8_protein;
-- 2697049 can be replaced with the virus taxon id, while ns8_orf8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_ns8_orf8_protein (
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
            AND ann.product = 'NS8 (ORF8 protein)'
            AND epi.protein_name = 'NS8 (ORF8 protein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT N (nucleocapsid phosphoprotein)
TRUNCATE public.epitope_2697049_n_nucleocapsid_phosphoprotei;
-- 2697049 can be replaced with the virus taxon id, while n_nucleocapsid_phosphoprotei can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_n_nucleocapsid_phosphoprotei (
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
            AND ann.product = 'N (nucleocapsid phosphoprotein)'
            AND epi.protein_name = 'N (nucleocapsid phosphoprotein)'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF10 protein
TRUNCATE public.epitope_2697049_orf10_protein;
-- 2697049 can be replaced with the virus taxon id, while orf10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_orf10_protein (
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
            AND ann.product = 'ORF10 protein'
            AND epi.protein_name = 'ORF10 protein'
            AND vir.taxon_id = 2697049)
ORDER BY epi.iedb_epitope_id;

