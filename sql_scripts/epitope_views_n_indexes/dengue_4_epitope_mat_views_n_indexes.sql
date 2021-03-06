-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN polyprotein
CREATE MATERIALIZED VIEW public.epitope_11070_polyprotein
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
             AND ann.product = 'polyprotein'
             AND epi.protein_name = 'polyprotein'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_polyprotein
    OWNER TO geco;

CREATE INDEX epi_11070_polyprotein__cell_type__idx
    ON public.epitope_11070_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_annotation_start__idx
    ON public.epitope_11070_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_annotation_stop__idx
    ON public.epitope_11070_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_frag_annotation_start__i
    ON public.epitope_11070_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_frag_annotation_stop__id
    ON public.epitope_11070_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__host_taxon_id__idx
    ON public.epitope_11070_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__host_taxon_name_lower__idx
    ON public.epitope_11070_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__iedb_epitope_id__idx
    ON public.epitope_11070_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__is_linear__idx
    ON public.epitope_11070_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__mhc_allele__idx
    ON public.epitope_11070_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__mhc_class_lower__idx
    ON public.epitope_11070_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__product_lower__idx
    ON public.epitope_11070_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__response_frequency_pos__idx
    ON public.epitope_11070_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__sequence_aa_alternative__idx
    ON public.epitope_11070_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__sequence_aa_original__idx
    ON public.epitope_11070_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__start_aa_original__idx
    ON public.epitope_11070_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__taxon_id__idx
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__taxon_name_lower__idx
    ON public.epitope_11070_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__variant_aa_length__idx
    ON public.epitope_11070_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__variant_aa_type__idx
    ON public.epitope_11070_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_cell_type__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_epi_start__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_epi_stop__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_is_linear__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_mhc_allele__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_product__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_resp_freq__i
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN anchored capsid protein ancC
CREATE MATERIALIZED VIEW public.epitope_11070_anchored_ca
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
             AND ann.product = 'anchored capsid protein ancC'
             AND epi.protein_name = 'anchored capsid protein ancC'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_anchored_ca
    OWNER TO geco;

CREATE INDEX epi_11070_anchored_ca__cell_type__idx
    ON public.epitope_11070_anchored_ca USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__epi_annotation_start__idx
    ON public.epitope_11070_anchored_ca USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__epi_annotation_stop__idx
    ON public.epitope_11070_anchored_ca USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__epi_frag_annotation_start__i
    ON public.epitope_11070_anchored_ca USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__epi_frag_annotation_stop__id
    ON public.epitope_11070_anchored_ca USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__host_taxon_id__idx
    ON public.epitope_11070_anchored_ca USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__host_taxon_name_lower__idx
    ON public.epitope_11070_anchored_ca USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__iedb_epitope_id__idx
    ON public.epitope_11070_anchored_ca USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__is_linear__idx
    ON public.epitope_11070_anchored_ca USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__mhc_allele__idx
    ON public.epitope_11070_anchored_ca USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__mhc_class_lower__idx
    ON public.epitope_11070_anchored_ca USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__product_lower__idx
    ON public.epitope_11070_anchored_ca USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__response_frequency_pos__idx
    ON public.epitope_11070_anchored_ca USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__sequence_aa_alternative__idx
    ON public.epitope_11070_anchored_ca USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__sequence_aa_original__idx
    ON public.epitope_11070_anchored_ca USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__start_aa_original__idx
    ON public.epitope_11070_anchored_ca USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__taxon_id__idx
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__taxon_name_lower__idx
    ON public.epitope_11070_anchored_ca USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__variant_aa_length__idx
    ON public.epitope_11070_anchored_ca USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__variant_aa_type__idx
    ON public.epitope_11070_anchored_ca USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_cell_type__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_epi_start__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_epi_stop__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_is_linear__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_mhc_allele__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_product__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_ca__virus_host_resp_freq__i
    ON public.epitope_11070_anchored_ca USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN capsid protein C
CREATE MATERIALIZED VIEW public.epitope_11070_capsid_prot
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
             AND ann.product = 'capsid protein C'
             AND epi.protein_name = 'capsid protein C'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_capsid_prot
    OWNER TO geco;

CREATE INDEX epi_11070_capsid_prot__cell_type__idx
    ON public.epitope_11070_capsid_prot USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__epi_annotation_start__idx
    ON public.epitope_11070_capsid_prot USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__epi_annotation_stop__idx
    ON public.epitope_11070_capsid_prot USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__epi_frag_annotation_start__i
    ON public.epitope_11070_capsid_prot USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__epi_frag_annotation_stop__id
    ON public.epitope_11070_capsid_prot USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__host_taxon_id__idx
    ON public.epitope_11070_capsid_prot USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__host_taxon_name_lower__idx
    ON public.epitope_11070_capsid_prot USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__iedb_epitope_id__idx
    ON public.epitope_11070_capsid_prot USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__is_linear__idx
    ON public.epitope_11070_capsid_prot USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__mhc_allele__idx
    ON public.epitope_11070_capsid_prot USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__mhc_class_lower__idx
    ON public.epitope_11070_capsid_prot USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__product_lower__idx
    ON public.epitope_11070_capsid_prot USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__response_frequency_pos__idx
    ON public.epitope_11070_capsid_prot USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__sequence_aa_alternative__idx
    ON public.epitope_11070_capsid_prot USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__sequence_aa_original__idx
    ON public.epitope_11070_capsid_prot USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__start_aa_original__idx
    ON public.epitope_11070_capsid_prot USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__taxon_id__idx
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__taxon_name_lower__idx
    ON public.epitope_11070_capsid_prot USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__variant_aa_length__idx
    ON public.epitope_11070_capsid_prot USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__variant_aa_type__idx
    ON public.epitope_11070_capsid_prot USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_cell_type__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_epi_start__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_epi_stop__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_is_linear__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_mhc_allele__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_product__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_prot__virus_host_resp_freq__i
    ON public.epitope_11070_capsid_prot USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN membrane glycoprotein precursor prM
CREATE MATERIALIZED VIEW public.epitope_11070_membrane_gl
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
             AND ann.product = 'membrane glycoprotein precursor prM'
             AND epi.protein_name = 'membrane glycoprotein precursor prM'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_membrane_gl
    OWNER TO geco;

CREATE INDEX epi_11070_membrane_gl__cell_type__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_annotation_start__idx
    ON public.epitope_11070_membrane_gl USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_annotation_stop__idx
    ON public.epitope_11070_membrane_gl USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_frag_annotation_start__i
    ON public.epitope_11070_membrane_gl USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_frag_annotation_stop__id
    ON public.epitope_11070_membrane_gl USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__host_taxon_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__host_taxon_name_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__iedb_epitope_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__is_linear__idx
    ON public.epitope_11070_membrane_gl USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__mhc_allele__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__mhc_class_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__product_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__response_frequency_pos__idx
    ON public.epitope_11070_membrane_gl USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__sequence_aa_alternative__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__sequence_aa_original__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__start_aa_original__idx
    ON public.epitope_11070_membrane_gl USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__taxon_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__taxon_name_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__variant_aa_length__idx
    ON public.epitope_11070_membrane_gl USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__variant_aa_type__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_cell_type__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_epi_start__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_epi_stop__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_is_linear__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_mhc_allele__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_product__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_resp_freq__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN protein pr
CREATE MATERIALIZED VIEW public.epitope_11070_protein_pr
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
             AND ann.product = 'protein pr'
             AND epi.protein_name = 'protein pr'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_protein_pr
    OWNER TO geco;

CREATE INDEX epi_11070_protein_pr__cell_type__idx
    ON public.epitope_11070_protein_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_annotation_start__idx
    ON public.epitope_11070_protein_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_annotation_stop__idx
    ON public.epitope_11070_protein_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_frag_annotation_start__i
    ON public.epitope_11070_protein_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_frag_annotation_stop__id
    ON public.epitope_11070_protein_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__host_taxon_id__idx
    ON public.epitope_11070_protein_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__host_taxon_name_lower__idx
    ON public.epitope_11070_protein_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__iedb_epitope_id__idx
    ON public.epitope_11070_protein_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__is_linear__idx
    ON public.epitope_11070_protein_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__mhc_allele__idx
    ON public.epitope_11070_protein_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__mhc_class_lower__idx
    ON public.epitope_11070_protein_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__product_lower__idx
    ON public.epitope_11070_protein_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__response_frequency_pos__idx
    ON public.epitope_11070_protein_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__sequence_aa_alternative__idx
    ON public.epitope_11070_protein_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__sequence_aa_original__idx
    ON public.epitope_11070_protein_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__start_aa_original__idx
    ON public.epitope_11070_protein_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__taxon_id__idx
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__taxon_name_lower__idx
    ON public.epitope_11070_protein_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__variant_aa_length__idx
    ON public.epitope_11070_protein_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__variant_aa_type__idx
    ON public.epitope_11070_protein_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_cell_type__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_epi_start__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_epi_stop__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_is_linear__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_mhc_allele__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_product__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_resp_freq__i
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN membrane glycoprotein M
CREATE MATERIALIZED VIEW public.epitope_11070_membrane_gl
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
             AND ann.product = 'membrane glycoprotein M'
             AND epi.protein_name = 'membrane glycoprotein M'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_membrane_gl
    OWNER TO geco;

CREATE INDEX epi_11070_membrane_gl__cell_type__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_annotation_start__idx
    ON public.epitope_11070_membrane_gl USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_annotation_stop__idx
    ON public.epitope_11070_membrane_gl USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_frag_annotation_start__i
    ON public.epitope_11070_membrane_gl USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__epi_frag_annotation_stop__id
    ON public.epitope_11070_membrane_gl USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__host_taxon_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__host_taxon_name_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__iedb_epitope_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__is_linear__idx
    ON public.epitope_11070_membrane_gl USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__mhc_allele__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__mhc_class_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__product_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__response_frequency_pos__idx
    ON public.epitope_11070_membrane_gl USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__sequence_aa_alternative__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__sequence_aa_original__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__start_aa_original__idx
    ON public.epitope_11070_membrane_gl USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__taxon_id__idx
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__taxon_name_lower__idx
    ON public.epitope_11070_membrane_gl USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__variant_aa_length__idx
    ON public.epitope_11070_membrane_gl USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__variant_aa_type__idx
    ON public.epitope_11070_membrane_gl USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_cell_type__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_epi_start__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_epi_stop__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_is_linear__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_mhc_allele__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_product__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_gl__virus_host_resp_freq__i
    ON public.epitope_11070_membrane_gl USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN envelope protein E
CREATE MATERIALIZED VIEW public.epitope_11070_envelope_pr
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
             AND ann.product = 'envelope protein E'
             AND epi.protein_name = 'envelope protein E'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_envelope_pr
    OWNER TO geco;

CREATE INDEX epi_11070_envelope_pr__cell_type__idx
    ON public.epitope_11070_envelope_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__epi_annotation_start__idx
    ON public.epitope_11070_envelope_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__epi_annotation_stop__idx
    ON public.epitope_11070_envelope_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__epi_frag_annotation_start__i
    ON public.epitope_11070_envelope_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__epi_frag_annotation_stop__id
    ON public.epitope_11070_envelope_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__host_taxon_id__idx
    ON public.epitope_11070_envelope_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__host_taxon_name_lower__idx
    ON public.epitope_11070_envelope_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__iedb_epitope_id__idx
    ON public.epitope_11070_envelope_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__is_linear__idx
    ON public.epitope_11070_envelope_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__mhc_allele__idx
    ON public.epitope_11070_envelope_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__mhc_class_lower__idx
    ON public.epitope_11070_envelope_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__product_lower__idx
    ON public.epitope_11070_envelope_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__response_frequency_pos__idx
    ON public.epitope_11070_envelope_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__sequence_aa_alternative__idx
    ON public.epitope_11070_envelope_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__sequence_aa_original__idx
    ON public.epitope_11070_envelope_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__start_aa_original__idx
    ON public.epitope_11070_envelope_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__taxon_id__idx
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__taxon_name_lower__idx
    ON public.epitope_11070_envelope_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__variant_aa_length__idx
    ON public.epitope_11070_envelope_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__variant_aa_type__idx
    ON public.epitope_11070_envelope_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_cell_type__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_epi_start__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_epi_stop__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_is_linear__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_mhc_allele__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_product__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_pr__virus_host_resp_freq__i
    ON public.epitope_11070_envelope_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS1
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS1'
             AND epi.protein_name = 'nonstructural protein NS1'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS2A
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS2A'
             AND epi.protein_name = 'nonstructural protein NS2A'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS2B
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS2B'
             AND epi.protein_name = 'nonstructural protein NS2B'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS3
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS3'
             AND epi.protein_name = 'nonstructural protein NS3'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS4A
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS4A'
             AND epi.protein_name = 'nonstructural protein NS4A'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN protein 2K
CREATE MATERIALIZED VIEW public.epitope_11070_protein_2k
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
             AND ann.product = 'protein 2K'
             AND epi.protein_name = 'protein 2K'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_protein_2k
    OWNER TO geco;

CREATE INDEX epi_11070_protein_2k__cell_type__idx
    ON public.epitope_11070_protein_2k USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_annotation_start__idx
    ON public.epitope_11070_protein_2k USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_annotation_stop__idx
    ON public.epitope_11070_protein_2k USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_frag_annotation_start__i
    ON public.epitope_11070_protein_2k USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_frag_annotation_stop__id
    ON public.epitope_11070_protein_2k USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__host_taxon_id__idx
    ON public.epitope_11070_protein_2k USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__host_taxon_name_lower__idx
    ON public.epitope_11070_protein_2k USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__iedb_epitope_id__idx
    ON public.epitope_11070_protein_2k USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__is_linear__idx
    ON public.epitope_11070_protein_2k USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__mhc_allele__idx
    ON public.epitope_11070_protein_2k USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__mhc_class_lower__idx
    ON public.epitope_11070_protein_2k USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__product_lower__idx
    ON public.epitope_11070_protein_2k USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__response_frequency_pos__idx
    ON public.epitope_11070_protein_2k USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__sequence_aa_alternative__idx
    ON public.epitope_11070_protein_2k USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__sequence_aa_original__idx
    ON public.epitope_11070_protein_2k USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__start_aa_original__idx
    ON public.epitope_11070_protein_2k USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__taxon_id__idx
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__taxon_name_lower__idx
    ON public.epitope_11070_protein_2k USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__variant_aa_length__idx
    ON public.epitope_11070_protein_2k USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__variant_aa_type__idx
    ON public.epitope_11070_protein_2k USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_cell_type__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_epi_start__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_epi_stop__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_is_linear__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_mhc_allele__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_product__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_resp_freq__i
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN nonstructural protein NS4B
CREATE MATERIALIZED VIEW public.epitope_11070_nonstructur
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
             AND ann.product = 'nonstructural protein NS4B'
             AND epi.protein_name = 'nonstructural protein NS4B'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_nonstructur
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructur__cell_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_start__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_annotation_stop__idx
    ON public.epitope_11070_nonstructur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_start__i
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__epi_frag_annotation_stop__id
    ON public.epitope_11070_nonstructur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__host_taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__iedb_epitope_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__is_linear__idx
    ON public.epitope_11070_nonstructur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_allele__idx
    ON public.epitope_11070_nonstructur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__mhc_class_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__product_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__response_frequency_pos__idx
    ON public.epitope_11070_nonstructur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_alternative__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__sequence_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__start_aa_original__idx
    ON public.epitope_11070_nonstructur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_id__idx
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__taxon_name_lower__idx
    ON public.epitope_11070_nonstructur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_length__idx
    ON public.epitope_11070_nonstructur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__variant_aa_type__idx
    ON public.epitope_11070_nonstructur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_cell_type__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_start__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_epi_stop__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_is_linear__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_mhc_allele__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_product__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructur__virus_host_resp_freq__i
    ON public.epitope_11070_nonstructur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN RNA-dependent RNA polymerase NS5
CREATE MATERIALIZED VIEW public.epitope_11070_rna_depende
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
             AND ann.product = 'RNA-dependent RNA polymerase NS5'
             AND epi.protein_name = 'RNA-dependent RNA polymerase NS5'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_rna_depende
    OWNER TO geco;

CREATE INDEX epi_11070_rna_depende__cell_type__idx
    ON public.epitope_11070_rna_depende USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__epi_annotation_start__idx
    ON public.epitope_11070_rna_depende USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__epi_annotation_stop__idx
    ON public.epitope_11070_rna_depende USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__epi_frag_annotation_start__i
    ON public.epitope_11070_rna_depende USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__epi_frag_annotation_stop__id
    ON public.epitope_11070_rna_depende USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__host_taxon_id__idx
    ON public.epitope_11070_rna_depende USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__host_taxon_name_lower__idx
    ON public.epitope_11070_rna_depende USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__iedb_epitope_id__idx
    ON public.epitope_11070_rna_depende USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__is_linear__idx
    ON public.epitope_11070_rna_depende USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__mhc_allele__idx
    ON public.epitope_11070_rna_depende USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__mhc_class_lower__idx
    ON public.epitope_11070_rna_depende USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__product_lower__idx
    ON public.epitope_11070_rna_depende USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__response_frequency_pos__idx
    ON public.epitope_11070_rna_depende USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__sequence_aa_alternative__idx
    ON public.epitope_11070_rna_depende USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__sequence_aa_original__idx
    ON public.epitope_11070_rna_depende USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__start_aa_original__idx
    ON public.epitope_11070_rna_depende USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__taxon_id__idx
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__taxon_name_lower__idx
    ON public.epitope_11070_rna_depende USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__variant_aa_length__idx
    ON public.epitope_11070_rna_depende USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__variant_aa_type__idx
    ON public.epitope_11070_rna_depende USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_cell_type__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_epi_start__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_epi_stop__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_is_linear__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_mhc_allele__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_product__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_depende__virus_host_resp_freq__i
    ON public.epitope_11070_rna_depende USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN sfRNA2
CREATE MATERIALIZED VIEW public.epitope_11070_sfrna2
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
             AND ann.product = 'sfRNA2'
             AND epi.protein_name = 'sfRNA2'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_sfrna2
    OWNER TO geco;

CREATE INDEX epi_11070_sfrna2__cell_type__idx
    ON public.epitope_11070_sfrna2 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__epi_annotation_start__idx
    ON public.epitope_11070_sfrna2 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__epi_annotation_stop__idx
    ON public.epitope_11070_sfrna2 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__epi_frag_annotation_start__i
    ON public.epitope_11070_sfrna2 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__epi_frag_annotation_stop__id
    ON public.epitope_11070_sfrna2 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__host_taxon_id__idx
    ON public.epitope_11070_sfrna2 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__host_taxon_name_lower__idx
    ON public.epitope_11070_sfrna2 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__iedb_epitope_id__idx
    ON public.epitope_11070_sfrna2 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__is_linear__idx
    ON public.epitope_11070_sfrna2 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__mhc_allele__idx
    ON public.epitope_11070_sfrna2 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__mhc_class_lower__idx
    ON public.epitope_11070_sfrna2 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__product_lower__idx
    ON public.epitope_11070_sfrna2 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__response_frequency_pos__idx
    ON public.epitope_11070_sfrna2 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__sequence_aa_alternative__idx
    ON public.epitope_11070_sfrna2 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__sequence_aa_original__idx
    ON public.epitope_11070_sfrna2 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__start_aa_original__idx
    ON public.epitope_11070_sfrna2 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__taxon_id__idx
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__taxon_name_lower__idx
    ON public.epitope_11070_sfrna2 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__variant_aa_length__idx
    ON public.epitope_11070_sfrna2 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__variant_aa_type__idx
    ON public.epitope_11070_sfrna2 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_cell_type__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_epi_start__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_epi_stop__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_is_linear__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_mhc_allele__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_product__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna2__virus_host_resp_freq__i
    ON public.epitope_11070_sfrna2 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN sfRNA3
CREATE MATERIALIZED VIEW public.epitope_11070_sfrna3
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
             AND ann.product = 'sfRNA3'
             AND epi.protein_name = 'sfRNA3'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_sfrna3
    OWNER TO geco;

CREATE INDEX epi_11070_sfrna3__cell_type__idx
    ON public.epitope_11070_sfrna3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__epi_annotation_start__idx
    ON public.epitope_11070_sfrna3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__epi_annotation_stop__idx
    ON public.epitope_11070_sfrna3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__epi_frag_annotation_start__i
    ON public.epitope_11070_sfrna3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__epi_frag_annotation_stop__id
    ON public.epitope_11070_sfrna3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__host_taxon_id__idx
    ON public.epitope_11070_sfrna3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__host_taxon_name_lower__idx
    ON public.epitope_11070_sfrna3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__iedb_epitope_id__idx
    ON public.epitope_11070_sfrna3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__is_linear__idx
    ON public.epitope_11070_sfrna3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__mhc_allele__idx
    ON public.epitope_11070_sfrna3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__mhc_class_lower__idx
    ON public.epitope_11070_sfrna3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__product_lower__idx
    ON public.epitope_11070_sfrna3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__response_frequency_pos__idx
    ON public.epitope_11070_sfrna3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__sequence_aa_alternative__idx
    ON public.epitope_11070_sfrna3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__sequence_aa_original__idx
    ON public.epitope_11070_sfrna3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__start_aa_original__idx
    ON public.epitope_11070_sfrna3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__taxon_id__idx
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__taxon_name_lower__idx
    ON public.epitope_11070_sfrna3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__variant_aa_length__idx
    ON public.epitope_11070_sfrna3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__variant_aa_type__idx
    ON public.epitope_11070_sfrna3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_cell_type__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_epi_start__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_epi_stop__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_is_linear__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_mhc_allele__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_product__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna3__virus_host_resp_freq__i
    ON public.epitope_11070_sfrna3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 11070 AND PROTEIN sfRNA4
CREATE MATERIALIZED VIEW public.epitope_11070_sfrna4
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
             AND ann.product = 'sfRNA4'
             AND epi.protein_name = 'sfRNA4'
             AND vir.taxon_id = 11070)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_11070_sfrna4
    OWNER TO geco;

CREATE INDEX epi_11070_sfrna4__cell_type__idx
    ON public.epitope_11070_sfrna4 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__epi_annotation_start__idx
    ON public.epitope_11070_sfrna4 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__epi_annotation_stop__idx
    ON public.epitope_11070_sfrna4 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__epi_frag_annotation_start__i
    ON public.epitope_11070_sfrna4 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__epi_frag_annotation_stop__id
    ON public.epitope_11070_sfrna4 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__host_taxon_id__idx
    ON public.epitope_11070_sfrna4 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__host_taxon_name_lower__idx
    ON public.epitope_11070_sfrna4 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__iedb_epitope_id__idx
    ON public.epitope_11070_sfrna4 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__is_linear__idx
    ON public.epitope_11070_sfrna4 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__mhc_allele__idx
    ON public.epitope_11070_sfrna4 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__mhc_class_lower__idx
    ON public.epitope_11070_sfrna4 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__product_lower__idx
    ON public.epitope_11070_sfrna4 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__response_frequency_pos__idx
    ON public.epitope_11070_sfrna4 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__sequence_aa_alternative__idx
    ON public.epitope_11070_sfrna4 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__sequence_aa_original__idx
    ON public.epitope_11070_sfrna4 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__start_aa_original__idx
    ON public.epitope_11070_sfrna4 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__taxon_id__idx
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__taxon_name_lower__idx
    ON public.epitope_11070_sfrna4 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__variant_aa_length__idx
    ON public.epitope_11070_sfrna4 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__variant_aa_type__idx
    ON public.epitope_11070_sfrna4 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_taxon_and_host_taxon_id__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_cell_type__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_epi_start__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_epi_stop__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_is_linear__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_mhc_allele__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_product__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_sfrna4__virus_host_resp_freq__i
    ON public.epitope_11070_sfrna4 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


