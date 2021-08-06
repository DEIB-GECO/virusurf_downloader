-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN ORF1ab polyprotein
CREATE MATERIALIZED VIEW public.epitope_2697049_orf1ab_poly
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_orf1ab_poly
    OWNER TO geco;

CREATE INDEX epi_2697049_orf1ab_poly__cell_type__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__epi_annotation_start__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__epi_annotation_stop__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__epi_frag_annotation_start__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__epi_frag_annotation_stop__id
    ON public.epitope_2697049_orf1ab_poly USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__host_taxon_id__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__host_taxon_name_lower__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__iedb_epitope_id__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__is_linear__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__mhc_allele__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__mhc_class_lower__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__product_lower__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__response_frequency_pos__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__sequence_aa_alternative__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__sequence_aa_original__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__start_aa_original__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__taxon_id__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__taxon_name_lower__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__variant_aa_length__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__variant_aa_type__idx
    ON public.epitope_2697049_orf1ab_poly USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_cell_type__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_epi_start__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_epi_stop__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_is_linear__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_mhc_allele__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_product__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_poly__virus_host_resp_freq__i
    ON public.epitope_2697049_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP12 (RNA-dependent RNA polymerase)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp12_rna_d
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp12_rna_d
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp12_rna_d__cell_type__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__epi_annotation_start__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__host_taxon_id__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__is_linear__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__mhc_allele__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__mhc_class_lower__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__product_lower__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__response_frequency_pos__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__sequence_aa_original__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__start_aa_original__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__taxon_id__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__taxon_name_lower__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__variant_aa_length__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__variant_aa_type__idx
    ON public.epitope_2697049_nsp12_rna_d USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_cell_type__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_epi_start__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_is_linear__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_product__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp12_rna_d USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP13 (helicase)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp13_helic
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp13_helic
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp13_helic__cell_type__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__epi_annotation_start__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp13_helic USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__host_taxon_id__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__is_linear__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__mhc_allele__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__mhc_class_lower__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__product_lower__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__response_frequency_pos__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__sequence_aa_original__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__start_aa_original__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__taxon_id__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__taxon_name_lower__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__variant_aa_length__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__variant_aa_type__idx
    ON public.epitope_2697049_nsp13_helic USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_cell_type__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_epi_start__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_is_linear__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_product__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helic__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp13_helic USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP14 (3'-to-5' exonuclease)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp14_3_to_
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp14_3_to_
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp14_3_to___cell_type__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___epi_annotation_start__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___epi_annotation_stop__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___host_taxon_id__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___iedb_epitope_id__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___is_linear__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___mhc_allele__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___mhc_class_lower__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___product_lower__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___response_frequency_pos__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___sequence_aa_original__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___start_aa_original__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___taxon_id__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___taxon_name_lower__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___variant_aa_length__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___variant_aa_type__idx
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_cell_type__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_epi_start__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_epi_stop__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_is_linear__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_product__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to___virus_host_resp_freq__i
    ON public.epitope_2697049_nsp14_3_to_ USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP15 (endoRNAse)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp15_endor
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp15_endor
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp15_endor__cell_type__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__epi_annotation_start__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp15_endor USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__host_taxon_id__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__is_linear__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__mhc_allele__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__mhc_class_lower__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__product_lower__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__response_frequency_pos__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__sequence_aa_original__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__start_aa_original__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__taxon_id__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__taxon_name_lower__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__variant_aa_length__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__variant_aa_type__idx
    ON public.epitope_2697049_nsp15_endor USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_cell_type__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_epi_start__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_is_linear__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_product__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endor__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp15_endor USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP16 (2'-O-ribose methyltransferase)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp16_2_o_r
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp16_2_o_r
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp16_2_o_r__cell_type__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__epi_annotation_start__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__host_taxon_id__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__is_linear__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__mhc_allele__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__mhc_class_lower__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__product_lower__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__response_frequency_pos__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__sequence_aa_original__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__start_aa_original__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__taxon_id__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__taxon_name_lower__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__variant_aa_length__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__variant_aa_type__idx
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_cell_type__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_epi_start__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_is_linear__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_product__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp16_2_o_r USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN ORF1a polyprotein
CREATE MATERIALIZED VIEW public.epitope_2697049_orf1a_polyp
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_orf1a_polyp
    OWNER TO geco;

CREATE INDEX epi_2697049_orf1a_polyp__cell_type__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__epi_annotation_start__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__epi_annotation_stop__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__epi_frag_annotation_start__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__epi_frag_annotation_stop__id
    ON public.epitope_2697049_orf1a_polyp USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__host_taxon_id__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__host_taxon_name_lower__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__iedb_epitope_id__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__is_linear__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__mhc_allele__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__mhc_class_lower__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__product_lower__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__response_frequency_pos__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__sequence_aa_alternative__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__sequence_aa_original__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__start_aa_original__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__taxon_id__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__taxon_name_lower__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__variant_aa_length__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__variant_aa_type__idx
    ON public.epitope_2697049_orf1a_polyp USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_cell_type__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_epi_start__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_epi_stop__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_is_linear__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_mhc_allele__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_product__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyp__virus_host_resp_freq__i
    ON public.epitope_2697049_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP1 (leader protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp1_leader
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp1_leader
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp1_leader__cell_type__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__epi_annotation_start__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp1_leader USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__host_taxon_id__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__is_linear__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__mhc_allele__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__mhc_class_lower__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__product_lower__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__response_frequency_pos__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__sequence_aa_original__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__start_aa_original__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__taxon_id__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__taxon_name_lower__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__variant_aa_length__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__variant_aa_type__idx
    ON public.epitope_2697049_nsp1_leader USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_cell_type__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_epi_start__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_is_linear__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_product__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp1_leader USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP2
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp2
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp2
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp2__cell_type__idx
    ON public.epitope_2697049_nsp2 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_annotation_start__idx
    ON public.epitope_2697049_nsp2 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp2 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp2 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp2 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__host_taxon_id__idx
    ON public.epitope_2697049_nsp2 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp2 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp2 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__is_linear__idx
    ON public.epitope_2697049_nsp2 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__mhc_allele__idx
    ON public.epitope_2697049_nsp2 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__mhc_class_lower__idx
    ON public.epitope_2697049_nsp2 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__product_lower__idx
    ON public.epitope_2697049_nsp2 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__response_frequency_pos__idx
    ON public.epitope_2697049_nsp2 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp2 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__sequence_aa_original__idx
    ON public.epitope_2697049_nsp2 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__start_aa_original__idx
    ON public.epitope_2697049_nsp2 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__taxon_id__idx
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__taxon_name_lower__idx
    ON public.epitope_2697049_nsp2 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__variant_aa_length__idx
    ON public.epitope_2697049_nsp2 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__variant_aa_type__idx
    ON public.epitope_2697049_nsp2 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_cell_type__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_epi_start__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_is_linear__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_product__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP3
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp3
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp3
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp3__cell_type__idx
    ON public.epitope_2697049_nsp3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_annotation_start__idx
    ON public.epitope_2697049_nsp3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__host_taxon_id__idx
    ON public.epitope_2697049_nsp3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__is_linear__idx
    ON public.epitope_2697049_nsp3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__mhc_allele__idx
    ON public.epitope_2697049_nsp3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__mhc_class_lower__idx
    ON public.epitope_2697049_nsp3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__product_lower__idx
    ON public.epitope_2697049_nsp3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__response_frequency_pos__idx
    ON public.epitope_2697049_nsp3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__sequence_aa_original__idx
    ON public.epitope_2697049_nsp3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__start_aa_original__idx
    ON public.epitope_2697049_nsp3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__taxon_id__idx
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__taxon_name_lower__idx
    ON public.epitope_2697049_nsp3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__variant_aa_length__idx
    ON public.epitope_2697049_nsp3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__variant_aa_type__idx
    ON public.epitope_2697049_nsp3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_cell_type__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_epi_start__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_is_linear__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_product__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP4
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp4
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp4
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp4__cell_type__idx
    ON public.epitope_2697049_nsp4 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_annotation_start__idx
    ON public.epitope_2697049_nsp4 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp4 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp4 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp4 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__host_taxon_id__idx
    ON public.epitope_2697049_nsp4 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp4 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp4 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__is_linear__idx
    ON public.epitope_2697049_nsp4 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__mhc_allele__idx
    ON public.epitope_2697049_nsp4 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__mhc_class_lower__idx
    ON public.epitope_2697049_nsp4 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__product_lower__idx
    ON public.epitope_2697049_nsp4 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__response_frequency_pos__idx
    ON public.epitope_2697049_nsp4 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp4 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__sequence_aa_original__idx
    ON public.epitope_2697049_nsp4 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__start_aa_original__idx
    ON public.epitope_2697049_nsp4 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__taxon_id__idx
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__taxon_name_lower__idx
    ON public.epitope_2697049_nsp4 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__variant_aa_length__idx
    ON public.epitope_2697049_nsp4 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__variant_aa_type__idx
    ON public.epitope_2697049_nsp4 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_cell_type__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_epi_start__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_is_linear__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_product__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP5 (3C-like proteinase)
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp5_3c_lik
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp5_3c_lik
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp5_3c_lik__cell_type__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__epi_annotation_start__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__host_taxon_id__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__is_linear__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__mhc_allele__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__mhc_class_lower__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__product_lower__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__response_frequency_pos__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__sequence_aa_original__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__start_aa_original__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__taxon_id__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__taxon_name_lower__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__variant_aa_length__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__variant_aa_type__idx
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_cell_type__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_epi_start__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_is_linear__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_product__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp5_3c_lik USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP6
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp6
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp6
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp6__cell_type__idx
    ON public.epitope_2697049_nsp6 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_annotation_start__idx
    ON public.epitope_2697049_nsp6 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp6 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp6 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp6 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__host_taxon_id__idx
    ON public.epitope_2697049_nsp6 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp6 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp6 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__is_linear__idx
    ON public.epitope_2697049_nsp6 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__mhc_allele__idx
    ON public.epitope_2697049_nsp6 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__mhc_class_lower__idx
    ON public.epitope_2697049_nsp6 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__product_lower__idx
    ON public.epitope_2697049_nsp6 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__response_frequency_pos__idx
    ON public.epitope_2697049_nsp6 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp6 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__sequence_aa_original__idx
    ON public.epitope_2697049_nsp6 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__start_aa_original__idx
    ON public.epitope_2697049_nsp6 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__taxon_id__idx
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__taxon_name_lower__idx
    ON public.epitope_2697049_nsp6 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__variant_aa_length__idx
    ON public.epitope_2697049_nsp6 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__variant_aa_type__idx
    ON public.epitope_2697049_nsp6 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_cell_type__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_epi_start__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_is_linear__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_product__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP7
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp7
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp7
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp7__cell_type__idx
    ON public.epitope_2697049_nsp7 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_annotation_start__idx
    ON public.epitope_2697049_nsp7 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp7 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp7 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp7 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__host_taxon_id__idx
    ON public.epitope_2697049_nsp7 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp7 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp7 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__is_linear__idx
    ON public.epitope_2697049_nsp7 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__mhc_allele__idx
    ON public.epitope_2697049_nsp7 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__mhc_class_lower__idx
    ON public.epitope_2697049_nsp7 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__product_lower__idx
    ON public.epitope_2697049_nsp7 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__response_frequency_pos__idx
    ON public.epitope_2697049_nsp7 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp7 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__sequence_aa_original__idx
    ON public.epitope_2697049_nsp7 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__start_aa_original__idx
    ON public.epitope_2697049_nsp7 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__taxon_id__idx
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__taxon_name_lower__idx
    ON public.epitope_2697049_nsp7 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__variant_aa_length__idx
    ON public.epitope_2697049_nsp7 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__variant_aa_type__idx
    ON public.epitope_2697049_nsp7 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_cell_type__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_epi_start__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_is_linear__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_product__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP8
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp8
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp8
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp8__cell_type__idx
    ON public.epitope_2697049_nsp8 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_annotation_start__idx
    ON public.epitope_2697049_nsp8 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp8 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp8 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp8 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__host_taxon_id__idx
    ON public.epitope_2697049_nsp8 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp8 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp8 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__is_linear__idx
    ON public.epitope_2697049_nsp8 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__mhc_allele__idx
    ON public.epitope_2697049_nsp8 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__mhc_class_lower__idx
    ON public.epitope_2697049_nsp8 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__product_lower__idx
    ON public.epitope_2697049_nsp8 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__response_frequency_pos__idx
    ON public.epitope_2697049_nsp8 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp8 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__sequence_aa_original__idx
    ON public.epitope_2697049_nsp8 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__start_aa_original__idx
    ON public.epitope_2697049_nsp8 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__taxon_id__idx
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__taxon_name_lower__idx
    ON public.epitope_2697049_nsp8 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__variant_aa_length__idx
    ON public.epitope_2697049_nsp8 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__variant_aa_type__idx
    ON public.epitope_2697049_nsp8 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_cell_type__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_epi_start__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_is_linear__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_product__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP9
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp9
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp9
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp9__cell_type__idx
    ON public.epitope_2697049_nsp9 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_annotation_start__idx
    ON public.epitope_2697049_nsp9 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp9 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp9 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp9 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__host_taxon_id__idx
    ON public.epitope_2697049_nsp9 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp9 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp9 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__is_linear__idx
    ON public.epitope_2697049_nsp9 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__mhc_allele__idx
    ON public.epitope_2697049_nsp9 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__mhc_class_lower__idx
    ON public.epitope_2697049_nsp9 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__product_lower__idx
    ON public.epitope_2697049_nsp9 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__response_frequency_pos__idx
    ON public.epitope_2697049_nsp9 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp9 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__sequence_aa_original__idx
    ON public.epitope_2697049_nsp9 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__start_aa_original__idx
    ON public.epitope_2697049_nsp9 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__taxon_id__idx
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__taxon_name_lower__idx
    ON public.epitope_2697049_nsp9 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__variant_aa_length__idx
    ON public.epitope_2697049_nsp9 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__variant_aa_type__idx
    ON public.epitope_2697049_nsp9 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_cell_type__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_epi_start__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_is_linear__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_product__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP10
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp10
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp10
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp10__cell_type__idx
    ON public.epitope_2697049_nsp10 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_annotation_start__idx
    ON public.epitope_2697049_nsp10 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp10 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp10 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp10 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__host_taxon_id__idx
    ON public.epitope_2697049_nsp10 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp10 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp10 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__is_linear__idx
    ON public.epitope_2697049_nsp10 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__mhc_allele__idx
    ON public.epitope_2697049_nsp10 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__mhc_class_lower__idx
    ON public.epitope_2697049_nsp10 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__product_lower__idx
    ON public.epitope_2697049_nsp10 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__response_frequency_pos__idx
    ON public.epitope_2697049_nsp10 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp10 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__sequence_aa_original__idx
    ON public.epitope_2697049_nsp10 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__start_aa_original__idx
    ON public.epitope_2697049_nsp10 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__taxon_id__idx
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__taxon_name_lower__idx
    ON public.epitope_2697049_nsp10 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__variant_aa_length__idx
    ON public.epitope_2697049_nsp10 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__variant_aa_type__idx
    ON public.epitope_2697049_nsp10 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_cell_type__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_epi_start__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_is_linear__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_product__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NSP11
CREATE MATERIALIZED VIEW public.epitope_2697049_nsp11
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_nsp11
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp11__cell_type__idx
    ON public.epitope_2697049_nsp11 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_annotation_start__idx
    ON public.epitope_2697049_nsp11 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_annotation_stop__idx
    ON public.epitope_2697049_nsp11 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_frag_annotation_start__i
    ON public.epitope_2697049_nsp11 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_frag_annotation_stop__id
    ON public.epitope_2697049_nsp11 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__host_taxon_id__idx
    ON public.epitope_2697049_nsp11 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__host_taxon_name_lower__idx
    ON public.epitope_2697049_nsp11 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__iedb_epitope_id__idx
    ON public.epitope_2697049_nsp11 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__is_linear__idx
    ON public.epitope_2697049_nsp11 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__mhc_allele__idx
    ON public.epitope_2697049_nsp11 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__mhc_class_lower__idx
    ON public.epitope_2697049_nsp11 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__product_lower__idx
    ON public.epitope_2697049_nsp11 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__response_frequency_pos__idx
    ON public.epitope_2697049_nsp11 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__sequence_aa_alternative__idx
    ON public.epitope_2697049_nsp11 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__sequence_aa_original__idx
    ON public.epitope_2697049_nsp11 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__start_aa_original__idx
    ON public.epitope_2697049_nsp11 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__taxon_id__idx
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__taxon_name_lower__idx
    ON public.epitope_2697049_nsp11 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__variant_aa_length__idx
    ON public.epitope_2697049_nsp11 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__variant_aa_type__idx
    ON public.epitope_2697049_nsp11 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_cell_type__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_epi_start__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_epi_stop__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_is_linear__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_mhc_allele__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_product__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_resp_freq__i
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN Spike (surface glycoprotein)
CREATE MATERIALIZED VIEW public.epitope_2697049_spike_surfa
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_spike_surfa
    OWNER TO geco;

CREATE INDEX epi_2697049_spike_surfa__cell_type__idx
    ON public.epitope_2697049_spike_surfa USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__epi_annotation_start__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__epi_annotation_stop__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__epi_frag_annotation_start__i
    ON public.epitope_2697049_spike_surfa USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__epi_frag_annotation_stop__id
    ON public.epitope_2697049_spike_surfa USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__host_taxon_id__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__host_taxon_name_lower__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__iedb_epitope_id__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__is_linear__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__mhc_allele__idx
    ON public.epitope_2697049_spike_surfa USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__mhc_class_lower__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__product_lower__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__response_frequency_pos__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__sequence_aa_alternative__idx
    ON public.epitope_2697049_spike_surfa USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__sequence_aa_original__idx
    ON public.epitope_2697049_spike_surfa USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__start_aa_original__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__taxon_id__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__taxon_name_lower__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__variant_aa_length__idx
    ON public.epitope_2697049_spike_surfa USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__variant_aa_type__idx
    ON public.epitope_2697049_spike_surfa USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_cell_type__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_epi_start__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_epi_stop__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_is_linear__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_mhc_allele__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_product__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surfa__virus_host_resp_freq__i
    ON public.epitope_2697049_spike_surfa USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NS3 (ORF3a protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_ns3_orf3a_p
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_ns3_orf3a_p
    OWNER TO geco;

CREATE INDEX epi_2697049_ns3_orf3a_p__cell_type__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__epi_annotation_start__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__epi_annotation_stop__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__epi_frag_annotation_start__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__epi_frag_annotation_stop__id
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__host_taxon_id__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__host_taxon_name_lower__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__iedb_epitope_id__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__is_linear__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__mhc_allele__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__mhc_class_lower__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__product_lower__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__response_frequency_pos__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__sequence_aa_alternative__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__sequence_aa_original__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__start_aa_original__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__taxon_id__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__taxon_name_lower__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__variant_aa_length__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__variant_aa_type__idx
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_cell_type__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_epi_start__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_epi_stop__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_is_linear__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_mhc_allele__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_product__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_resp_freq__i
    ON public.epitope_2697049_ns3_orf3a_p USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN E (envelope protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_e_envelope_
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_e_envelope_
    OWNER TO geco;

CREATE INDEX epi_2697049_e_envelope___cell_type__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___epi_annotation_start__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___epi_annotation_stop__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___epi_frag_annotation_start__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___epi_frag_annotation_stop__id
    ON public.epitope_2697049_e_envelope_ USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___host_taxon_id__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___host_taxon_name_lower__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___iedb_epitope_id__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___is_linear__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___mhc_allele__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___mhc_class_lower__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___product_lower__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___response_frequency_pos__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___sequence_aa_alternative__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___sequence_aa_original__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___start_aa_original__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___taxon_id__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___taxon_name_lower__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___variant_aa_length__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___variant_aa_type__idx
    ON public.epitope_2697049_e_envelope_ USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_cell_type__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_epi_start__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_epi_stop__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_is_linear__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_mhc_allele__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_product__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope___virus_host_resp_freq__i
    ON public.epitope_2697049_e_envelope_ USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN M (membrane glycoprotein)
CREATE MATERIALIZED VIEW public.epitope_2697049_m_membrane_
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_m_membrane_
    OWNER TO geco;

CREATE INDEX epi_2697049_m_membrane___cell_type__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___epi_annotation_start__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___epi_annotation_stop__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___epi_frag_annotation_start__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___epi_frag_annotation_stop__id
    ON public.epitope_2697049_m_membrane_ USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___host_taxon_id__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___host_taxon_name_lower__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___iedb_epitope_id__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___is_linear__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___mhc_allele__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___mhc_class_lower__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___product_lower__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___response_frequency_pos__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___sequence_aa_alternative__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___sequence_aa_original__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___start_aa_original__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___taxon_id__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___taxon_name_lower__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___variant_aa_length__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___variant_aa_type__idx
    ON public.epitope_2697049_m_membrane_ USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_cell_type__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_epi_start__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_epi_stop__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_is_linear__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_mhc_allele__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_product__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane___virus_host_resp_freq__i
    ON public.epitope_2697049_m_membrane_ USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NS6 (ORF6 protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_ns6_orf6_pr
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_ns6_orf6_pr
    OWNER TO geco;

CREATE INDEX epi_2697049_ns6_orf6_pr__cell_type__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__epi_annotation_start__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__epi_annotation_stop__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__epi_frag_annotation_start__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__epi_frag_annotation_stop__id
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__host_taxon_id__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__host_taxon_name_lower__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__iedb_epitope_id__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__is_linear__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__mhc_allele__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__mhc_class_lower__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__product_lower__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__response_frequency_pos__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__sequence_aa_alternative__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__sequence_aa_original__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__start_aa_original__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__taxon_id__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__taxon_name_lower__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__variant_aa_length__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__variant_aa_type__idx
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_cell_type__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_epi_start__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_epi_stop__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_is_linear__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_mhc_allele__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_product__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_resp_freq__i
    ON public.epitope_2697049_ns6_orf6_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NS7a (ORF7a protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_ns7a_orf7a_
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_ns7a_orf7a_
    OWNER TO geco;

CREATE INDEX epi_2697049_ns7a_orf7a___cell_type__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___epi_annotation_start__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___epi_annotation_stop__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___epi_frag_annotation_start__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___epi_frag_annotation_stop__id
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___host_taxon_id__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___host_taxon_name_lower__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___iedb_epitope_id__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___is_linear__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___mhc_allele__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___mhc_class_lower__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___product_lower__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___response_frequency_pos__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___sequence_aa_alternative__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___sequence_aa_original__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___start_aa_original__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___taxon_id__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___taxon_name_lower__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___variant_aa_length__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___variant_aa_type__idx
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_cell_type__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_epi_start__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_epi_stop__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_is_linear__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_mhc_allele__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_product__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_resp_freq__i
    ON public.epitope_2697049_ns7a_orf7a_ USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NS7b (ORF7b)
CREATE MATERIALIZED VIEW public.epitope_2697049_ns7b_orf7b
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_ns7b_orf7b
    OWNER TO geco;

CREATE INDEX epi_2697049_ns7b_orf7b__cell_type__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_annotation_start__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_annotation_stop__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_annotation_start__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_annotation_stop__id
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__host_taxon_id__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__host_taxon_name_lower__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__iedb_epitope_id__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__is_linear__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__mhc_allele__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__mhc_class_lower__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__product_lower__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__response_frequency_pos__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__sequence_aa_alternative__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__sequence_aa_original__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__start_aa_original__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__taxon_id__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__taxon_name_lower__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_length__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_type__idx
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_cell_type__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_start__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_stop__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_is_linear__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_mhc_allele__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_product__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_resp_freq__i
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN NS8 (ORF8 protein)
CREATE MATERIALIZED VIEW public.epitope_2697049_ns8_orf8_pr
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_ns8_orf8_pr
    OWNER TO geco;

CREATE INDEX epi_2697049_ns8_orf8_pr__cell_type__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__epi_annotation_start__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__epi_annotation_stop__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__epi_frag_annotation_start__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__epi_frag_annotation_stop__id
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__host_taxon_id__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__host_taxon_name_lower__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__iedb_epitope_id__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__is_linear__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__mhc_allele__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__mhc_class_lower__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__product_lower__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__response_frequency_pos__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__sequence_aa_alternative__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__sequence_aa_original__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__start_aa_original__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__taxon_id__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__taxon_name_lower__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__variant_aa_length__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__variant_aa_type__idx
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_cell_type__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_epi_start__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_epi_stop__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_is_linear__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_mhc_allele__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_product__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_resp_freq__i
    ON public.epitope_2697049_ns8_orf8_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN N (nucleocapsid phosphoprotein)
CREATE MATERIALIZED VIEW public.epitope_2697049_n_nucleocap
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_n_nucleocap
    OWNER TO geco;

CREATE INDEX epi_2697049_n_nucleocap__cell_type__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__epi_annotation_start__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__epi_annotation_stop__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__epi_frag_annotation_start__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__epi_frag_annotation_stop__id
    ON public.epitope_2697049_n_nucleocap USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__host_taxon_id__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__host_taxon_name_lower__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__iedb_epitope_id__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__is_linear__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__mhc_allele__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__mhc_class_lower__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__product_lower__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__response_frequency_pos__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__sequence_aa_alternative__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__sequence_aa_original__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__start_aa_original__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__taxon_id__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__taxon_name_lower__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__variant_aa_length__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__variant_aa_type__idx
    ON public.epitope_2697049_n_nucleocap USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_cell_type__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_epi_start__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_epi_stop__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_is_linear__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_mhc_allele__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_product__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocap__virus_host_resp_freq__i
    ON public.epitope_2697049_n_nucleocap USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 2697049 AND PROTEIN ORF10 protein
CREATE MATERIALIZED VIEW public.epitope_2697049_orf10_prote
    TABLESPACE default_ts
AS
 SELECT DISTINCT epi.iedb_epitope_id,
    epi.epitope_iri,
    epi.cell_type,
    epi.mhc_class,
    epi.mhc_allele,
    epi.response_frequency_pos,
    epi.epi_annotation_start,
    epi.epi_annotation_stop,
    epi.is_linear,
    epi.assay_type,
    epif.epi_fragment_sequence,
    epif.epi_frag_annotation_start,
    epif.epi_frag_annotation_stop,
    vir.taxon_id,
    vir.taxon_name,
    hspec.host_taxon_id,
    hspec.host_taxon_name,
    seq.sequence_id,
    ann.product,
    amin.aminoacid_variant_id,
    amin.start_aa_original,
    amin.sequence_aa_original,
    amin.sequence_aa_alternative,
    amin.variant_aa_length,
    amin.variant_aa_type
   FROM (((((((epitope epi
     JOIN epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
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
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_2697049_orf10_prote
    OWNER TO geco;

CREATE INDEX epi_2697049_orf10_prote__cell_type__idx
    ON public.epitope_2697049_orf10_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__epi_annotation_start__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__epi_annotation_stop__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__epi_frag_annotation_start__i
    ON public.epitope_2697049_orf10_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__epi_frag_annotation_stop__id
    ON public.epitope_2697049_orf10_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__host_taxon_id__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__host_taxon_name_lower__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__iedb_epitope_id__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__is_linear__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__mhc_allele__idx
    ON public.epitope_2697049_orf10_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__mhc_class_lower__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__product_lower__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__response_frequency_pos__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__sequence_aa_alternative__idx
    ON public.epitope_2697049_orf10_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__sequence_aa_original__idx
    ON public.epitope_2697049_orf10_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__start_aa_original__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__taxon_id__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__taxon_name_lower__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__variant_aa_length__idx
    ON public.epitope_2697049_orf10_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__variant_aa_type__idx
    ON public.epitope_2697049_orf10_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_cell_type__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_epi_start__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_epi_stop__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_is_linear__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_mhc_allele__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_product__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_prote__virus_host_resp_freq__i
    ON public.epitope_2697049_orf10_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


