-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN 1AB polyprotein
CREATE MATERIALIZED VIEW public.epitope_1335626_1ab_polypro
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
             AND ann.product = '1AB polyprotein'
             AND epi.protein_name = '1AB polyprotein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_1ab_polypro
    OWNER TO geco;

CREATE INDEX epi_1335626_1ab_polypro__cell_type__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__epi_annotation_start__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__epi_annotation_stop__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__epi_frag_annotation_start__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__epi_frag_annotation_stop__id
    ON public.epitope_1335626_1ab_polypro USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__host_taxon_id__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__host_taxon_name_lower__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__iedb_epitope_id__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__is_linear__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__mhc_allele__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__mhc_class_lower__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__product_lower__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__response_frequency_pos__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__sequence_aa_alternative__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__sequence_aa_original__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__start_aa_original__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__taxon_id__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__taxon_name_lower__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__variant_aa_length__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__variant_aa_type__idx
    ON public.epitope_1335626_1ab_polypro USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_cell_type__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_epi_start__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_epi_stop__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_is_linear__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_mhc_allele__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_product__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polypro__virus_host_resp_freq__i
    ON public.epitope_1335626_1ab_polypro USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN RNA-dependent RNA polymerase
CREATE MATERIALIZED VIEW public.epitope_1335626_rna_depende
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
             AND ann.product = 'RNA-dependent RNA polymerase'
             AND epi.protein_name = 'RNA-dependent RNA polymerase'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_rna_depende
    OWNER TO geco;

CREATE INDEX epi_1335626_rna_depende__cell_type__idx
    ON public.epitope_1335626_rna_depende USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__epi_annotation_start__idx
    ON public.epitope_1335626_rna_depende USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__epi_annotation_stop__idx
    ON public.epitope_1335626_rna_depende USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__epi_frag_annotation_start__i
    ON public.epitope_1335626_rna_depende USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__epi_frag_annotation_stop__id
    ON public.epitope_1335626_rna_depende USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__host_taxon_id__idx
    ON public.epitope_1335626_rna_depende USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__host_taxon_name_lower__idx
    ON public.epitope_1335626_rna_depende USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__iedb_epitope_id__idx
    ON public.epitope_1335626_rna_depende USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__is_linear__idx
    ON public.epitope_1335626_rna_depende USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__mhc_allele__idx
    ON public.epitope_1335626_rna_depende USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__mhc_class_lower__idx
    ON public.epitope_1335626_rna_depende USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__product_lower__idx
    ON public.epitope_1335626_rna_depende USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__response_frequency_pos__idx
    ON public.epitope_1335626_rna_depende USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__sequence_aa_alternative__idx
    ON public.epitope_1335626_rna_depende USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__sequence_aa_original__idx
    ON public.epitope_1335626_rna_depende USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__start_aa_original__idx
    ON public.epitope_1335626_rna_depende USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__taxon_id__idx
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__taxon_name_lower__idx
    ON public.epitope_1335626_rna_depende USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__variant_aa_length__idx
    ON public.epitope_1335626_rna_depende USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__variant_aa_type__idx
    ON public.epitope_1335626_rna_depende USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_cell_type__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_epi_start__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_epi_stop__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_is_linear__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_mhc_allele__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_product__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_depende__virus_host_resp_freq__i
    ON public.epitope_1335626_rna_depende USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN Hel
CREATE MATERIALIZED VIEW public.epitope_1335626_hel
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
             AND ann.product = 'Hel'
             AND epi.protein_name = 'Hel'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_hel
    OWNER TO geco;

CREATE INDEX epi_1335626_hel__cell_type__idx
    ON public.epitope_1335626_hel USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_annotation_start__idx
    ON public.epitope_1335626_hel USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_annotation_stop__idx
    ON public.epitope_1335626_hel USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_frag_annotation_start__i
    ON public.epitope_1335626_hel USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_frag_annotation_stop__id
    ON public.epitope_1335626_hel USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__host_taxon_id__idx
    ON public.epitope_1335626_hel USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__host_taxon_name_lower__idx
    ON public.epitope_1335626_hel USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__iedb_epitope_id__idx
    ON public.epitope_1335626_hel USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__is_linear__idx
    ON public.epitope_1335626_hel USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__mhc_allele__idx
    ON public.epitope_1335626_hel USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__mhc_class_lower__idx
    ON public.epitope_1335626_hel USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__product_lower__idx
    ON public.epitope_1335626_hel USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__response_frequency_pos__idx
    ON public.epitope_1335626_hel USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__sequence_aa_alternative__idx
    ON public.epitope_1335626_hel USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__sequence_aa_original__idx
    ON public.epitope_1335626_hel USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__start_aa_original__idx
    ON public.epitope_1335626_hel USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__taxon_id__idx
    ON public.epitope_1335626_hel USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__taxon_name_lower__idx
    ON public.epitope_1335626_hel USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__variant_aa_length__idx
    ON public.epitope_1335626_hel USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__variant_aa_type__idx
    ON public.epitope_1335626_hel USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_cell_type__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_epi_start__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_epi_stop__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_is_linear__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_mhc_allele__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_product__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_resp_freq__i
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN ExoN
CREATE MATERIALIZED VIEW public.epitope_1335626_exon
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
             AND ann.product = 'ExoN'
             AND epi.protein_name = 'ExoN'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_exon
    OWNER TO geco;

CREATE INDEX epi_1335626_exon__cell_type__idx
    ON public.epitope_1335626_exon USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_annotation_start__idx
    ON public.epitope_1335626_exon USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_annotation_stop__idx
    ON public.epitope_1335626_exon USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_frag_annotation_start__i
    ON public.epitope_1335626_exon USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_frag_annotation_stop__id
    ON public.epitope_1335626_exon USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__host_taxon_id__idx
    ON public.epitope_1335626_exon USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__host_taxon_name_lower__idx
    ON public.epitope_1335626_exon USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__iedb_epitope_id__idx
    ON public.epitope_1335626_exon USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__is_linear__idx
    ON public.epitope_1335626_exon USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__mhc_allele__idx
    ON public.epitope_1335626_exon USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__mhc_class_lower__idx
    ON public.epitope_1335626_exon USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__product_lower__idx
    ON public.epitope_1335626_exon USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__response_frequency_pos__idx
    ON public.epitope_1335626_exon USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__sequence_aa_alternative__idx
    ON public.epitope_1335626_exon USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__sequence_aa_original__idx
    ON public.epitope_1335626_exon USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__start_aa_original__idx
    ON public.epitope_1335626_exon USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__taxon_id__idx
    ON public.epitope_1335626_exon USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__taxon_name_lower__idx
    ON public.epitope_1335626_exon USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__variant_aa_length__idx
    ON public.epitope_1335626_exon USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__variant_aa_type__idx
    ON public.epitope_1335626_exon USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_cell_type__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_epi_start__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_epi_stop__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_is_linear__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_mhc_allele__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_product__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_resp_freq__i
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN NendoU
CREATE MATERIALIZED VIEW public.epitope_1335626_nendou
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
             AND ann.product = 'NendoU'
             AND epi.protein_name = 'NendoU'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nendou
    OWNER TO geco;

CREATE INDEX epi_1335626_nendou__cell_type__idx
    ON public.epitope_1335626_nendou USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_annotation_start__idx
    ON public.epitope_1335626_nendou USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_annotation_stop__idx
    ON public.epitope_1335626_nendou USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_frag_annotation_start__i
    ON public.epitope_1335626_nendou USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nendou USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__host_taxon_id__idx
    ON public.epitope_1335626_nendou USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__host_taxon_name_lower__idx
    ON public.epitope_1335626_nendou USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__iedb_epitope_id__idx
    ON public.epitope_1335626_nendou USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__is_linear__idx
    ON public.epitope_1335626_nendou USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__mhc_allele__idx
    ON public.epitope_1335626_nendou USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__mhc_class_lower__idx
    ON public.epitope_1335626_nendou USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__product_lower__idx
    ON public.epitope_1335626_nendou USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__response_frequency_pos__idx
    ON public.epitope_1335626_nendou USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__sequence_aa_alternative__idx
    ON public.epitope_1335626_nendou USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__sequence_aa_original__idx
    ON public.epitope_1335626_nendou USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__start_aa_original__idx
    ON public.epitope_1335626_nendou USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__taxon_id__idx
    ON public.epitope_1335626_nendou USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__taxon_name_lower__idx
    ON public.epitope_1335626_nendou USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__variant_aa_length__idx
    ON public.epitope_1335626_nendou USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__variant_aa_type__idx
    ON public.epitope_1335626_nendou USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_cell_type__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_epi_start__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_epi_stop__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_is_linear__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_mhc_allele__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_product__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_resp_freq__i
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN 2'-O-methyltransferase
CREATE MATERIALIZED VIEW public.epitope_1335626_2_o_methylt
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
             AND ann.product = '2''-O-methyltransferase'
             AND epi.protein_name = '2''-O-methyltransferase'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_2_o_methylt
    OWNER TO geco;

CREATE INDEX epi_1335626_2_o_methylt__cell_type__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__epi_annotation_start__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__epi_annotation_stop__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__epi_frag_annotation_start__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__epi_frag_annotation_stop__id
    ON public.epitope_1335626_2_o_methylt USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__host_taxon_id__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__host_taxon_name_lower__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__iedb_epitope_id__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__is_linear__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__mhc_allele__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__mhc_class_lower__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__product_lower__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__response_frequency_pos__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__sequence_aa_alternative__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__sequence_aa_original__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__start_aa_original__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__taxon_id__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__taxon_name_lower__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__variant_aa_length__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__variant_aa_type__idx
    ON public.epitope_1335626_2_o_methylt USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_cell_type__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_epi_start__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_epi_stop__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_is_linear__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_mhc_allele__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_product__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methylt__virus_host_resp_freq__i
    ON public.epitope_1335626_2_o_methylt USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN 1A polyprotein
CREATE MATERIALIZED VIEW public.epitope_1335626_1a_polyprot
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
             AND ann.product = '1A polyprotein'
             AND epi.protein_name = '1A polyprotein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_1a_polyprot
    OWNER TO geco;

CREATE INDEX epi_1335626_1a_polyprot__cell_type__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__epi_annotation_start__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__epi_annotation_stop__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__epi_frag_annotation_start__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__epi_frag_annotation_stop__id
    ON public.epitope_1335626_1a_polyprot USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__host_taxon_id__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__host_taxon_name_lower__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__iedb_epitope_id__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__is_linear__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__mhc_allele__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__mhc_class_lower__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__product_lower__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__response_frequency_pos__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__sequence_aa_alternative__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__sequence_aa_original__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__start_aa_original__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__taxon_id__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__taxon_name_lower__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__variant_aa_length__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__variant_aa_type__idx
    ON public.epitope_1335626_1a_polyprot USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_cell_type__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_epi_start__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_epi_stop__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_is_linear__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_mhc_allele__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_product__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprot__virus_host_resp_freq__i
    ON public.epitope_1335626_1a_polyprot USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp1 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp1_protei
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
             AND ann.product = 'nsp1 protein'
             AND epi.protein_name = 'nsp1 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp1_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp1_protei__cell_type__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp1_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__is_linear__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__product_lower__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__taxon_id__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp1_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_product__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp1_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp2 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp2_protei
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
             AND ann.product = 'nsp2 protein'
             AND epi.protein_name = 'nsp2 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp2_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp2_protei__cell_type__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp2_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__is_linear__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__product_lower__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__taxon_id__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp2_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_product__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp2_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp3 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp3_protei
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
             AND ann.product = 'nsp3 protein'
             AND epi.protein_name = 'nsp3 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp3_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp3_protei__cell_type__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp3_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__is_linear__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__product_lower__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__taxon_id__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp3_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_product__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp3_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp4 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp4_protei
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
             AND ann.product = 'nsp4 protein'
             AND epi.protein_name = 'nsp4 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp4_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp4_protei__cell_type__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp4_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__is_linear__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__product_lower__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__taxon_id__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp4_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_product__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp4_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp5 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp5_protei
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
             AND ann.product = 'nsp5 protein'
             AND epi.protein_name = 'nsp5 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp5_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp5_protei__cell_type__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp5_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__is_linear__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__product_lower__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__taxon_id__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp5_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_product__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp5_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp6 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp6_protei
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
             AND ann.product = 'nsp6 protein'
             AND epi.protein_name = 'nsp6 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp6_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp6_protei__cell_type__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp6_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__is_linear__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__product_lower__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__taxon_id__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp6_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_product__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp6_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp7 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp7_protei
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
             AND ann.product = 'nsp7 protein'
             AND epi.protein_name = 'nsp7 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp7_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp7_protei__cell_type__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp7_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__is_linear__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__product_lower__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__taxon_id__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp7_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_product__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp7_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp8 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp8_protei
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
             AND ann.product = 'nsp8 protein'
             AND epi.protein_name = 'nsp8 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp8_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp8_protei__cell_type__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp8_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__is_linear__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__product_lower__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__taxon_id__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp8_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_product__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp8_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp9 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp9_protei
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
             AND ann.product = 'nsp9 protein'
             AND epi.protein_name = 'nsp9 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp9_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp9_protei__cell_type__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__epi_annotation_start__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp9_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__host_taxon_id__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__is_linear__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__mhc_allele__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__mhc_class_lower__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__product_lower__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__response_frequency_pos__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__sequence_aa_original__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__start_aa_original__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__taxon_id__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__taxon_name_lower__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__variant_aa_length__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__variant_aa_type__idx
    ON public.epitope_1335626_nsp9_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_cell_type__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_epi_start__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_is_linear__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_product__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp9_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp10 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp10_prote
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
             AND ann.product = 'nsp10 protein'
             AND epi.protein_name = 'nsp10 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp10_prote
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp10_prote__cell_type__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__epi_annotation_start__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp10_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__host_taxon_id__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__is_linear__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__mhc_allele__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__mhc_class_lower__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__product_lower__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__response_frequency_pos__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__sequence_aa_original__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__start_aa_original__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__taxon_id__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__taxon_name_lower__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__variant_aa_length__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__variant_aa_type__idx
    ON public.epitope_1335626_nsp10_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_cell_type__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_epi_start__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_is_linear__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_product__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_prote__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp10_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nsp11 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nsp11_prote
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
             AND ann.product = 'nsp11 protein'
             AND epi.protein_name = 'nsp11 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nsp11_prote
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp11_prote__cell_type__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__epi_annotation_start__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__epi_annotation_stop__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__epi_frag_annotation_start__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nsp11_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__host_taxon_id__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__host_taxon_name_lower__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__iedb_epitope_id__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__is_linear__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__mhc_allele__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__mhc_class_lower__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__product_lower__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__response_frequency_pos__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__sequence_aa_alternative__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__sequence_aa_original__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__start_aa_original__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__taxon_id__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__taxon_name_lower__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__variant_aa_length__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__variant_aa_type__idx
    ON public.epitope_1335626_nsp11_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_cell_type__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_epi_start__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_epi_stop__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_is_linear__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_mhc_allele__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_product__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_prote__virus_host_resp_freq__i
    ON public.epitope_1335626_nsp11_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN spike protein
CREATE MATERIALIZED VIEW public.epitope_1335626_spike_prote
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
             AND ann.product = 'spike protein'
             AND epi.protein_name = 'spike protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_spike_prote
    OWNER TO geco;

CREATE INDEX epi_1335626_spike_prote__cell_type__idx
    ON public.epitope_1335626_spike_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__epi_annotation_start__idx
    ON public.epitope_1335626_spike_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__epi_annotation_stop__idx
    ON public.epitope_1335626_spike_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__epi_frag_annotation_start__i
    ON public.epitope_1335626_spike_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__epi_frag_annotation_stop__id
    ON public.epitope_1335626_spike_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__host_taxon_id__idx
    ON public.epitope_1335626_spike_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__host_taxon_name_lower__idx
    ON public.epitope_1335626_spike_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__iedb_epitope_id__idx
    ON public.epitope_1335626_spike_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__is_linear__idx
    ON public.epitope_1335626_spike_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__mhc_allele__idx
    ON public.epitope_1335626_spike_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__mhc_class_lower__idx
    ON public.epitope_1335626_spike_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__product_lower__idx
    ON public.epitope_1335626_spike_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__response_frequency_pos__idx
    ON public.epitope_1335626_spike_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__sequence_aa_alternative__idx
    ON public.epitope_1335626_spike_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__sequence_aa_original__idx
    ON public.epitope_1335626_spike_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__start_aa_original__idx
    ON public.epitope_1335626_spike_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__taxon_id__idx
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__taxon_name_lower__idx
    ON public.epitope_1335626_spike_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__variant_aa_length__idx
    ON public.epitope_1335626_spike_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__variant_aa_type__idx
    ON public.epitope_1335626_spike_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_cell_type__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_epi_start__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_epi_stop__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_is_linear__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_mhc_allele__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_product__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_prote__virus_host_resp_freq__i
    ON public.epitope_1335626_spike_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN NS3 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_ns3_protein
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
             AND ann.product = 'NS3 protein'
             AND epi.protein_name = 'NS3 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_ns3_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns3_protein__cell_type__idx
    ON public.epitope_1335626_ns3_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_annotation_start__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_annotation_stop__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_frag_annotation_start__i
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_frag_annotation_stop__id
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__host_taxon_id__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__host_taxon_name_lower__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__iedb_epitope_id__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__is_linear__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__mhc_allele__idx
    ON public.epitope_1335626_ns3_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__mhc_class_lower__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__product_lower__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__response_frequency_pos__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__sequence_aa_alternative__idx
    ON public.epitope_1335626_ns3_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__sequence_aa_original__idx
    ON public.epitope_1335626_ns3_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__start_aa_original__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__taxon_id__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__taxon_name_lower__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__variant_aa_length__idx
    ON public.epitope_1335626_ns3_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__variant_aa_type__idx
    ON public.epitope_1335626_ns3_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_cell_type__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_epi_start__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_epi_stop__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_is_linear__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_mhc_allele__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_product__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_resp_freq__i
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN NS4A protein
CREATE MATERIALIZED VIEW public.epitope_1335626_ns4a_protei
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
             AND ann.product = 'NS4A protein'
             AND epi.protein_name = 'NS4A protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_ns4a_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_ns4a_protei__cell_type__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__epi_annotation_start__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_ns4a_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__host_taxon_id__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__is_linear__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__mhc_allele__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__mhc_class_lower__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__product_lower__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__response_frequency_pos__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__sequence_aa_original__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__start_aa_original__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__taxon_id__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__taxon_name_lower__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__variant_aa_length__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__variant_aa_type__idx
    ON public.epitope_1335626_ns4a_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_cell_type__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_epi_start__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_is_linear__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_product__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_ns4a_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN NS4B protein
CREATE MATERIALIZED VIEW public.epitope_1335626_ns4b_protei
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
             AND ann.product = 'NS4B protein'
             AND epi.protein_name = 'NS4B protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_ns4b_protei
    OWNER TO geco;

CREATE INDEX epi_1335626_ns4b_protei__cell_type__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__epi_annotation_start__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__epi_annotation_stop__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__epi_frag_annotation_start__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__epi_frag_annotation_stop__id
    ON public.epitope_1335626_ns4b_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__host_taxon_id__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__host_taxon_name_lower__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__iedb_epitope_id__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__is_linear__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__mhc_allele__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__mhc_class_lower__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__product_lower__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__response_frequency_pos__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__sequence_aa_alternative__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__sequence_aa_original__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__start_aa_original__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__taxon_id__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__taxon_name_lower__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__variant_aa_length__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__variant_aa_type__idx
    ON public.epitope_1335626_ns4b_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_cell_type__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_epi_start__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_epi_stop__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_is_linear__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_mhc_allele__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_product__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protei__virus_host_resp_freq__i
    ON public.epitope_1335626_ns4b_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN NS5 protein
CREATE MATERIALIZED VIEW public.epitope_1335626_ns5_protein
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
             AND ann.product = 'NS5 protein'
             AND epi.protein_name = 'NS5 protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_ns5_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns5_protein__cell_type__idx
    ON public.epitope_1335626_ns5_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_annotation_start__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_annotation_stop__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_frag_annotation_start__i
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_frag_annotation_stop__id
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__host_taxon_id__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__host_taxon_name_lower__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__iedb_epitope_id__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__is_linear__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__mhc_allele__idx
    ON public.epitope_1335626_ns5_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__mhc_class_lower__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__product_lower__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__response_frequency_pos__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__sequence_aa_alternative__idx
    ON public.epitope_1335626_ns5_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__sequence_aa_original__idx
    ON public.epitope_1335626_ns5_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__start_aa_original__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__taxon_id__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__taxon_name_lower__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__variant_aa_length__idx
    ON public.epitope_1335626_ns5_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__variant_aa_type__idx
    ON public.epitope_1335626_ns5_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_cell_type__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_epi_start__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_epi_stop__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_is_linear__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_mhc_allele__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_product__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_resp_freq__i
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN envelope protein
CREATE MATERIALIZED VIEW public.epitope_1335626_envelope_pr
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
             AND ann.product = 'envelope protein'
             AND epi.protein_name = 'envelope protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_envelope_pr
    OWNER TO geco;

CREATE INDEX epi_1335626_envelope_pr__cell_type__idx
    ON public.epitope_1335626_envelope_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__epi_annotation_start__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__epi_annotation_stop__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__epi_frag_annotation_start__i
    ON public.epitope_1335626_envelope_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__epi_frag_annotation_stop__id
    ON public.epitope_1335626_envelope_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__host_taxon_id__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__host_taxon_name_lower__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__iedb_epitope_id__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__is_linear__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__mhc_allele__idx
    ON public.epitope_1335626_envelope_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__mhc_class_lower__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__product_lower__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__response_frequency_pos__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__sequence_aa_alternative__idx
    ON public.epitope_1335626_envelope_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__sequence_aa_original__idx
    ON public.epitope_1335626_envelope_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__start_aa_original__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__taxon_id__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__taxon_name_lower__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__variant_aa_length__idx
    ON public.epitope_1335626_envelope_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__variant_aa_type__idx
    ON public.epitope_1335626_envelope_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_cell_type__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_epi_start__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_epi_stop__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_is_linear__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_mhc_allele__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_product__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_pr__virus_host_resp_freq__i
    ON public.epitope_1335626_envelope_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN membrane protein
CREATE MATERIALIZED VIEW public.epitope_1335626_membrane_pr
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
             AND ann.product = 'membrane protein'
             AND epi.protein_name = 'membrane protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_membrane_pr
    OWNER TO geco;

CREATE INDEX epi_1335626_membrane_pr__cell_type__idx
    ON public.epitope_1335626_membrane_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__epi_annotation_start__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__epi_annotation_stop__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__epi_frag_annotation_start__i
    ON public.epitope_1335626_membrane_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__epi_frag_annotation_stop__id
    ON public.epitope_1335626_membrane_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__host_taxon_id__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__host_taxon_name_lower__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__iedb_epitope_id__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__is_linear__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__mhc_allele__idx
    ON public.epitope_1335626_membrane_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__mhc_class_lower__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__product_lower__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__response_frequency_pos__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__sequence_aa_alternative__idx
    ON public.epitope_1335626_membrane_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__sequence_aa_original__idx
    ON public.epitope_1335626_membrane_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__start_aa_original__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__taxon_id__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__taxon_name_lower__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__variant_aa_length__idx
    ON public.epitope_1335626_membrane_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__variant_aa_type__idx
    ON public.epitope_1335626_membrane_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_cell_type__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_epi_start__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_epi_stop__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_is_linear__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_mhc_allele__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_product__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_pr__virus_host_resp_freq__i
    ON public.epitope_1335626_membrane_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN nucleocapsid protein
CREATE MATERIALIZED VIEW public.epitope_1335626_nucleocapsi
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
             AND ann.product = 'nucleocapsid protein'
             AND epi.protein_name = 'nucleocapsid protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_nucleocapsi
    OWNER TO geco;

CREATE INDEX epi_1335626_nucleocapsi__cell_type__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__epi_annotation_start__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__epi_annotation_stop__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__epi_frag_annotation_start__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__epi_frag_annotation_stop__id
    ON public.epitope_1335626_nucleocapsi USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__host_taxon_id__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__host_taxon_name_lower__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__iedb_epitope_id__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__is_linear__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__mhc_allele__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__mhc_class_lower__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__product_lower__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__response_frequency_pos__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__sequence_aa_alternative__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__sequence_aa_original__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__start_aa_original__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__taxon_id__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__taxon_name_lower__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__variant_aa_length__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__variant_aa_type__idx
    ON public.epitope_1335626_nucleocapsi USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_cell_type__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_epi_start__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_epi_stop__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_is_linear__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_mhc_allele__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_product__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsi__virus_host_resp_freq__i
    ON public.epitope_1335626_nucleocapsi USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 1335626 AND PROTEIN ORF8b protein
CREATE MATERIALIZED VIEW public.epitope_1335626_orf8b_prote
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
             AND ann.product = 'ORF8b protein'
             AND epi.protein_name = 'ORF8b protein'
             AND vir.taxon_id = 1335626)
  ORDER BY epi.iedb_epitope_id
WITH NO DATA;

ALTER TABLE public.epitope_1335626_orf8b_prote
    OWNER TO geco;

CREATE INDEX epi_1335626_orf8b_prote__cell_type__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__epi_annotation_start__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__epi_annotation_stop__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__epi_frag_annotation_start__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__epi_frag_annotation_stop__id
    ON public.epitope_1335626_orf8b_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__host_taxon_id__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__host_taxon_name_lower__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__iedb_epitope_id__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__is_linear__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__mhc_allele__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__mhc_class_lower__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__product_lower__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__response_frequency_pos__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__sequence_aa_alternative__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__sequence_aa_original__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__start_aa_original__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__taxon_id__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__taxon_name_lower__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__variant_aa_length__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__variant_aa_type__idx
    ON public.epitope_1335626_orf8b_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_cell_type__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_epi_start__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_epi_stop__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_is_linear__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_mhc_allele__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_product__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_prote__virus_host_resp_freq__i
    ON public.epitope_1335626_orf8b_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


