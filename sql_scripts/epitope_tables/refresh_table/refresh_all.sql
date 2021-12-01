-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1ab polyprotein
TRUNCATE public.epitope_2697049_orf1ab_polyprotein;
-- 2697049 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2697049_orf1ab_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp12_rna_dependent_rna_poly VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp13_helicase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp14_3_to_5_exonuclease VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp15_endornase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp16_2_o_ribose_methyltrans VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_orf1a_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp1_leader_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp2 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp3 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp4 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp5_3c_like_proteinase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp6 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp7 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp8 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp9 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp10 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_nsp11 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_spike_surface_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_ns3_orf3a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_e_envelope_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_m_membrane_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_ns6_orf6_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_ns7a_orf7a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_ns7b_orf7b VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_ns8_orf8_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_n_nucleocapsid_phosphoprotei VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_2697049_orf10_protein VALUES(
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
) 
SELECT DISTINCT 
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

-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1ab polyprotein
TRUNCATE public.epitope_694009_orf1ab_polyprotein;
-- 694009 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_694009_orf1ab_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_helicase_ntpase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_3_to_5_exonuclease VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_endoribonuclease VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_2_o_mtase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf1a_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp1 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp2 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp3 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp4 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_3c_like_protease VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp6 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp7 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp8 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp9 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nsp10 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_ndp11 VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf3a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf3b_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_small_envelope_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_membrane_glycoprotein_m VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf6_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf7a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf7b_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf8a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf8b_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_nucleocapsid_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf9b_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_694009_orf9a_protein VALUES(
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
) 
SELECT DISTINCT 
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

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT polyprotein
TRUNCATE public.epitope_11053_polyprotein;
-- 11053 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'polyprotein'
            AND epi.protein_name = 'polyprotein'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT anchored capsid protein ancC
TRUNCATE public.epitope_11053_anchored_capsid_protein_ancc;
-- 11053 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_anchored_capsid_protein_ancc VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'anchored capsid protein ancC'
            AND epi.protein_name = 'anchored capsid protein ancC'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT capsid protein C
TRUNCATE public.epitope_11053_capsid_protein_c;
-- 11053 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_capsid_protein_c VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'capsid protein C'
            AND epi.protein_name = 'capsid protein C'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT membrane glycoprotein precursor prM
TRUNCATE public.epitope_11053_membrane_glycoprotein_precur;
-- 11053 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_membrane_glycoprotein_precur VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'membrane glycoprotein precursor prM'
            AND epi.protein_name = 'membrane glycoprotein precursor prM'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT protein pr
TRUNCATE public.epitope_11053_protein_pr;
-- 11053 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_protein_pr VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein pr'
            AND epi.protein_name = 'protein pr'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT membrane glycoprotein M
TRUNCATE public.epitope_11053_membrane_glycoprotein_m;
-- 11053 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_membrane_glycoprotein_m VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT envelope protein E
TRUNCATE public.epitope_11053_envelope_protein_e;
-- 11053 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_envelope_protein_e VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'envelope protein E'
            AND epi.protein_name = 'envelope protein E'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS1
TRUNCATE public.epitope_11053_nonstructural_protein_ns1;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns1 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS1'
            AND epi.protein_name = 'nonstructural protein NS1'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS2A
TRUNCATE public.epitope_11053_nonstructural_protein_ns2a;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns2a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2A'
            AND epi.protein_name = 'nonstructural protein NS2A'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS2B
TRUNCATE public.epitope_11053_nonstructural_protein_ns2b;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns2b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2B'
            AND epi.protein_name = 'nonstructural protein NS2B'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS3
TRUNCATE public.epitope_11053_nonstructural_protein_ns3;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns3 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS3'
            AND epi.protein_name = 'nonstructural protein NS3'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS4A
TRUNCATE public.epitope_11053_nonstructural_protein_ns4a;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns4a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4A'
            AND epi.protein_name = 'nonstructural protein NS4A'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT protein 2K
TRUNCATE public.epitope_11053_protein_2k;
-- 11053 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_protein_2k VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein 2K'
            AND epi.protein_name = 'protein 2K'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS4B
TRUNCATE public.epitope_11053_nonstructural_protein_ns4b;
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_nonstructural_protein_ns4b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4B'
            AND epi.protein_name = 'nonstructural protein NS4B'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_1 and PROT RNA-dependent RNA polymerase NS5
TRUNCATE public.epitope_11053_rna_dependent_rna_polymerase;
-- 11053 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11053_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'RNA-dependent RNA polymerase NS5'
            AND epi.protein_name = 'RNA-dependent RNA polymerase NS5'
            AND vir.taxon_id = 11053)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT polyprotein
TRUNCATE public.epitope_11060_polyprotein;
-- 11060 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'polyprotein'
            AND epi.protein_name = 'polyprotein'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT anchored capsid protein ancC
TRUNCATE public.epitope_11060_anchored_capsid_protein_ancc;
-- 11060 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_anchored_capsid_protein_ancc VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'anchored capsid protein ancC'
            AND epi.protein_name = 'anchored capsid protein ancC'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT capsid protein C
TRUNCATE public.epitope_11060_capsid_protein_c;
-- 11060 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_capsid_protein_c VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'capsid protein C'
            AND epi.protein_name = 'capsid protein C'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT membrane glycoprotein precursor prM
TRUNCATE public.epitope_11060_membrane_glycoprotein_precur;
-- 11060 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_membrane_glycoprotein_precur VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'membrane glycoprotein precursor prM'
            AND epi.protein_name = 'membrane glycoprotein precursor prM'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT protein pr
TRUNCATE public.epitope_11060_protein_pr;
-- 11060 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_protein_pr VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein pr'
            AND epi.protein_name = 'protein pr'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT membrane glycoprotein M
TRUNCATE public.epitope_11060_membrane_glycoprotein_m;
-- 11060 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_membrane_glycoprotein_m VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT envelope protein E
TRUNCATE public.epitope_11060_envelope_protein_e;
-- 11060 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_envelope_protein_e VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'envelope protein E'
            AND epi.protein_name = 'envelope protein E'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS1
TRUNCATE public.epitope_11060_nonstructural_protein_ns1;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns1 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS1'
            AND epi.protein_name = 'nonstructural protein NS1'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS2A
TRUNCATE public.epitope_11060_nonstructural_protein_ns2a;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns2a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2A'
            AND epi.protein_name = 'nonstructural protein NS2A'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS2B
TRUNCATE public.epitope_11060_nonstructural_protein_ns2b;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns2b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2B'
            AND epi.protein_name = 'nonstructural protein NS2B'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS3
TRUNCATE public.epitope_11060_nonstructural_protein_ns3;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns3 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS3'
            AND epi.protein_name = 'nonstructural protein NS3'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS4A
TRUNCATE public.epitope_11060_nonstructural_protein_ns4a;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns4a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4A'
            AND epi.protein_name = 'nonstructural protein NS4A'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT protein 2K
TRUNCATE public.epitope_11060_protein_2k;
-- 11060 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_protein_2k VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein 2K'
            AND epi.protein_name = 'protein 2K'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS4B
TRUNCATE public.epitope_11060_nonstructural_protein_ns4b;
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_nonstructural_protein_ns4b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4B'
            AND epi.protein_name = 'nonstructural protein NS4B'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_2 and PROT RNA-dependent RNA polymerase NS5
TRUNCATE public.epitope_11060_rna_dependent_rna_polymerase;
-- 11060 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11060_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'RNA-dependent RNA polymerase NS5'
            AND epi.protein_name = 'RNA-dependent RNA polymerase NS5'
            AND vir.taxon_id = 11060)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT polyprotein
TRUNCATE public.epitope_11069_polyprotein;
-- 11069 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'polyprotein'
            AND epi.protein_name = 'polyprotein'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT anchored capsid protein ancC
TRUNCATE public.epitope_11069_anchored_capsid_protein_ancc;
-- 11069 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_anchored_capsid_protein_ancc VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'anchored capsid protein ancC'
            AND epi.protein_name = 'anchored capsid protein ancC'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT capsid protein C
TRUNCATE public.epitope_11069_capsid_protein_c;
-- 11069 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_capsid_protein_c VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'capsid protein C'
            AND epi.protein_name = 'capsid protein C'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT membrane glycoprotein precursor prM
TRUNCATE public.epitope_11069_membrane_glycoprotein_precur;
-- 11069 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_membrane_glycoprotein_precur VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'membrane glycoprotein precursor prM'
            AND epi.protein_name = 'membrane glycoprotein precursor prM'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT protein pr
TRUNCATE public.epitope_11069_protein_pr;
-- 11069 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_protein_pr VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein pr'
            AND epi.protein_name = 'protein pr'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT membrane glycoprotein M
TRUNCATE public.epitope_11069_membrane_glycoprotein_m;
-- 11069 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_membrane_glycoprotein_m VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT envelope protein E
TRUNCATE public.epitope_11069_envelope_protein_e;
-- 11069 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_envelope_protein_e VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'envelope protein E'
            AND epi.protein_name = 'envelope protein E'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS1
TRUNCATE public.epitope_11069_nonstructural_protein_ns1;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns1 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS1'
            AND epi.protein_name = 'nonstructural protein NS1'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS2A
TRUNCATE public.epitope_11069_nonstructural_protein_ns2a;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns2a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2A'
            AND epi.protein_name = 'nonstructural protein NS2A'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS2B
TRUNCATE public.epitope_11069_nonstructural_protein_ns2b;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns2b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2B'
            AND epi.protein_name = 'nonstructural protein NS2B'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS3
TRUNCATE public.epitope_11069_nonstructural_protein_ns3;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns3 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS3'
            AND epi.protein_name = 'nonstructural protein NS3'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS4A
TRUNCATE public.epitope_11069_nonstructural_protein_ns4a;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns4a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4A'
            AND epi.protein_name = 'nonstructural protein NS4A'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT protein 2K
TRUNCATE public.epitope_11069_protein_2k;
-- 11069 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_protein_2k VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein 2K'
            AND epi.protein_name = 'protein 2K'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS4B
TRUNCATE public.epitope_11069_nonstructural_protein_ns4b;
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_nonstructural_protein_ns4b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4B'
            AND epi.protein_name = 'nonstructural protein NS4B'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_3 and PROT RNA-dependent RNA polymerase NS5
TRUNCATE public.epitope_11069_rna_dependent_rna_polymerase;
-- 11069 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11069_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'RNA-dependent RNA polymerase NS5'
            AND epi.protein_name = 'RNA-dependent RNA polymerase NS5'
            AND vir.taxon_id = 11069)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT polyprotein
TRUNCATE public.epitope_11070_polyprotein;
-- 11070 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'polyprotein'
            AND epi.protein_name = 'polyprotein'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT anchored capsid protein ancC
TRUNCATE public.epitope_11070_anchored_capsid_protein_ancc;
-- 11070 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_anchored_capsid_protein_ancc VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'anchored capsid protein ancC'
            AND epi.protein_name = 'anchored capsid protein ancC'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT capsid protein C
TRUNCATE public.epitope_11070_capsid_protein_c;
-- 11070 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_capsid_protein_c VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'capsid protein C'
            AND epi.protein_name = 'capsid protein C'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein precursor prM
TRUNCATE public.epitope_11070_membrane_glycoprotein_precur;
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_membrane_glycoprotein_precur VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'membrane glycoprotein precursor prM'
            AND epi.protein_name = 'membrane glycoprotein precursor prM'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT protein pr
TRUNCATE public.epitope_11070_protein_pr;
-- 11070 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_protein_pr VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein pr'
            AND epi.protein_name = 'protein pr'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein M
TRUNCATE public.epitope_11070_membrane_glycoprotein_m;
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_membrane_glycoprotein_m VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT envelope protein E
TRUNCATE public.epitope_11070_envelope_protein_e;
-- 11070 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_envelope_protein_e VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'envelope protein E'
            AND epi.protein_name = 'envelope protein E'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS1
TRUNCATE public.epitope_11070_nonstructural_protein_ns1;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns1 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS1'
            AND epi.protein_name = 'nonstructural protein NS1'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2A
TRUNCATE public.epitope_11070_nonstructural_protein_ns2a;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns2a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2A'
            AND epi.protein_name = 'nonstructural protein NS2A'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2B
TRUNCATE public.epitope_11070_nonstructural_protein_ns2b;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns2b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS2B'
            AND epi.protein_name = 'nonstructural protein NS2B'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS3
TRUNCATE public.epitope_11070_nonstructural_protein_ns3;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns3 VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS3'
            AND epi.protein_name = 'nonstructural protein NS3'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4A
TRUNCATE public.epitope_11070_nonstructural_protein_ns4a;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns4a VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4A'
            AND epi.protein_name = 'nonstructural protein NS4A'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT protein 2K
TRUNCATE public.epitope_11070_protein_2k;
-- 11070 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_protein_2k VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'protein 2K'
            AND epi.protein_name = 'protein 2K'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4B
TRUNCATE public.epitope_11070_nonstructural_protein_ns4b;
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_nonstructural_protein_ns4b VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'nonstructural protein NS4B'
            AND epi.protein_name = 'nonstructural protein NS4B'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT RNA-dependent RNA polymerase NS5
TRUNCATE public.epitope_11070_rna_dependent_rna_polymerase;
-- 11070 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_11070_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND ann.product = 'RNA-dependent RNA polymerase NS5'
            AND epi.protein_name = 'RNA-dependent RNA polymerase NS5'
            AND vir.taxon_id = 11070)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 1AB polyprotein
TRUNCATE public.epitope_1335626_1ab_polyprotein;
-- 1335626 can be replaced with the virus taxon id, while 1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_1335626_1ab_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_hel VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_exon VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nendou VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_2_o_methyltransferase VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_1a_polyprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp1_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp2_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp3_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp4_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp5_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp6_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp7_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp8_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp9_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp10_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nsp11_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_spike_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_ns3_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_ns4a_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_ns4b_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_ns5_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_envelope_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_membrane_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_nucleocapsid_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_1335626_orf8b_protein VALUES(
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
) 
SELECT DISTINCT 
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

-- CREATE TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_186538_nucleoprotein;
-- 186538 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186538_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
INSERT INTO public.epitope_186538_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_186540_nucleoprotein;
-- 186540 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_186540_polymerase_complex_protein;
-- 186540 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_186540_matrix_protein;
-- 186540 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_186540_spike_glycoprotein;
-- 186540 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_186540_small_secreted_glycoprotein;
-- 186540 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_186540_second_secreted_glycoprotein;
-- 186540 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_186540_minor_nucleoprotein;
-- 186540 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_186540_membrane_associated_protein;
-- 186540 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_186540_rna_dependent_rna_polymerase;
-- 186540 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186540_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186540)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_186539_nucleoprotein;
-- 186539 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_186539_polymerase_complex_protein;
-- 186539 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_186539_matrix_protein;
-- 186539 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_186539_spike_glycoprotein;
-- 186539 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_186539_small_secreted_glycoprotein;
-- 186539 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_186539_second_secreted_glycoprotein;
-- 186539 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_186539_minor_nucleoprotein;
-- 186539 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_186539_membrane_associated_protein;
-- 186539 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_186539_rna_dependent_rna_polymerase;
-- 186539 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186539_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186539)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_565995_nucleoprotein;
-- 565995 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_565995_polymerase_complex_protein;
-- 565995 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_565995_matrix_protein;
-- 565995 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_565995_spike_glycoprotein;
-- 565995 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_565995_small_secreted_glycoprotein;
-- 565995 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_565995_second_secreted_glycoprotein;
-- 565995 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_565995_minor_nucleoprotein;
-- 565995 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_565995_membrane_associated_protein;
-- 565995 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_565995_rna_dependent_rna_polymerase;
-- 565995 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_565995_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 565995)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_2010960_nucleoprotein;
-- 2010960 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_2010960_polymerase_complex_protein;
-- 2010960 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_2010960_matrix_protein;
-- 2010960 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_2010960_spike_glycoprotein;
-- 2010960 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_2010960_small_secreted_glycoprotein;
-- 2010960 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_2010960_second_secreted_glycoprotein;
-- 2010960 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_2010960_minor_nucleoprotein;
-- 2010960 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_2010960_membrane_associated_protein;
-- 2010960 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_2010960_rna_dependent_rna_polymerase;
-- 2010960 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_2010960_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 2010960)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT nucleoprotein
TRUNCATE public.epitope_186541_nucleoprotein;
-- 186541 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT polymerase complex protein
TRUNCATE public.epitope_186541_polymerase_complex_protein;
-- 186541 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_polymerase_complex_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT matrix protein
TRUNCATE public.epitope_186541_matrix_protein;
-- 186541 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_matrix_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT spike glycoprotein
TRUNCATE public.epitope_186541_spike_glycoprotein;
-- 186541 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_spike_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT small secreted glycoprotein
TRUNCATE public.epitope_186541_small_secreted_glycoprotein;
-- 186541 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_small_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT second secreted glycoprotein
TRUNCATE public.epitope_186541_second_secreted_glycoprotein;
-- 186541 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_second_secreted_glycoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT minor nucleoprotein
TRUNCATE public.epitope_186541_minor_nucleoprotein;
-- 186541 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_minor_nucleoprotein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT membrane-associated protein
TRUNCATE public.epitope_186541_membrane_associated_protein;
-- 186541 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_membrane_associated_protein VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

-- CREATE TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT RNA-dependent RNA polymerase
TRUNCATE public.epitope_186541_rna_dependent_rna_polymerase;
-- 186541 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
INSERT INTO public.epitope_186541_rna_dependent_rna_polymerase VALUES(
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
) 
SELECT DISTINCT 
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
            AND vir.taxon_id = 186541)
ORDER BY epi.iedb_epitope_id;

