-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF1ab polyprotein
CREATE MATERIALIZED VIEW public.epitope_694009_orf1ab_poly
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf1ab_poly
    OWNER TO geco;

CREATE INDEX epi_694009_orf1ab_poly__cell_type__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__epi_annotation_start__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__epi_annotation_stop__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__epi_frag_annotation_start__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf1ab_poly USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__host_taxon_id__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__host_taxon_name_lower__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__iedb_epitope_id__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__is_linear__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__mhc_allele__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__mhc_class_lower__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__product_lower__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__response_frequency_pos__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__sequence_aa_alternative__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__sequence_aa_original__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__start_aa_original__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__taxon_id__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__taxon_name_lower__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__variant_aa_length__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__variant_aa_type__idx
    ON public.epitope_694009_orf1ab_poly USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_cell_type__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_epi_start__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_epi_stop__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_is_linear__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_mhc_allele__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_product__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_poly__virus_host_resp_freq__i
    ON public.epitope_694009_orf1ab_poly USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN RNA-dependent RNA polymerase
CREATE MATERIALIZED VIEW public.epitope_694009_rna_depende
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_rna_depende
    OWNER TO geco;

CREATE INDEX epi_694009_rna_depende__cell_type__idx
    ON public.epitope_694009_rna_depende USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__epi_annotation_start__idx
    ON public.epitope_694009_rna_depende USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__epi_annotation_stop__idx
    ON public.epitope_694009_rna_depende USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__epi_frag_annotation_start__i
    ON public.epitope_694009_rna_depende USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__epi_frag_annotation_stop__id
    ON public.epitope_694009_rna_depende USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__host_taxon_id__idx
    ON public.epitope_694009_rna_depende USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__host_taxon_name_lower__idx
    ON public.epitope_694009_rna_depende USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__iedb_epitope_id__idx
    ON public.epitope_694009_rna_depende USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__is_linear__idx
    ON public.epitope_694009_rna_depende USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__mhc_allele__idx
    ON public.epitope_694009_rna_depende USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__mhc_class_lower__idx
    ON public.epitope_694009_rna_depende USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__product_lower__idx
    ON public.epitope_694009_rna_depende USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__response_frequency_pos__idx
    ON public.epitope_694009_rna_depende USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__sequence_aa_alternative__idx
    ON public.epitope_694009_rna_depende USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__sequence_aa_original__idx
    ON public.epitope_694009_rna_depende USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__start_aa_original__idx
    ON public.epitope_694009_rna_depende USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__taxon_id__idx
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__taxon_name_lower__idx
    ON public.epitope_694009_rna_depende USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__variant_aa_length__idx
    ON public.epitope_694009_rna_depende USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__variant_aa_type__idx
    ON public.epitope_694009_rna_depende USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_cell_type__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_epi_start__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_epi_stop__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_is_linear__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_mhc_allele__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_product__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_depende__virus_host_resp_freq__i
    ON public.epitope_694009_rna_depende USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN helicase/NTPase
CREATE MATERIALIZED VIEW public.epitope_694009_helicase_nt
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
             AND ann.product = 'helicase/NTPase'
             AND epi.protein_name = 'helicase/NTPase'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_helicase_nt
    OWNER TO geco;

CREATE INDEX epi_694009_helicase_nt__cell_type__idx
    ON public.epitope_694009_helicase_nt USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__epi_annotation_start__idx
    ON public.epitope_694009_helicase_nt USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__epi_annotation_stop__idx
    ON public.epitope_694009_helicase_nt USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__epi_frag_annotation_start__i
    ON public.epitope_694009_helicase_nt USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__epi_frag_annotation_stop__id
    ON public.epitope_694009_helicase_nt USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__host_taxon_id__idx
    ON public.epitope_694009_helicase_nt USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__host_taxon_name_lower__idx
    ON public.epitope_694009_helicase_nt USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__iedb_epitope_id__idx
    ON public.epitope_694009_helicase_nt USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__is_linear__idx
    ON public.epitope_694009_helicase_nt USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__mhc_allele__idx
    ON public.epitope_694009_helicase_nt USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__mhc_class_lower__idx
    ON public.epitope_694009_helicase_nt USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__product_lower__idx
    ON public.epitope_694009_helicase_nt USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__response_frequency_pos__idx
    ON public.epitope_694009_helicase_nt USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__sequence_aa_alternative__idx
    ON public.epitope_694009_helicase_nt USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__sequence_aa_original__idx
    ON public.epitope_694009_helicase_nt USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__start_aa_original__idx
    ON public.epitope_694009_helicase_nt USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__taxon_id__idx
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__taxon_name_lower__idx
    ON public.epitope_694009_helicase_nt USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__variant_aa_length__idx
    ON public.epitope_694009_helicase_nt USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__variant_aa_type__idx
    ON public.epitope_694009_helicase_nt USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_cell_type__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_epi_start__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_epi_stop__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_is_linear__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_mhc_allele__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_product__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_nt__virus_host_resp_freq__i
    ON public.epitope_694009_helicase_nt USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN 3' to 5' exonuclease
CREATE MATERIALIZED VIEW public.epitope_694009_3_to_5_exon
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
             AND ann.product = '3' to 5' exonuclease'
             AND epi.protein_name = '3' to 5' exonuclease'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_3_to_5_exon
    OWNER TO geco;

CREATE INDEX epi_694009_3_to_5_exon__cell_type__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__epi_annotation_start__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__epi_annotation_stop__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__epi_frag_annotation_start__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__epi_frag_annotation_stop__id
    ON public.epitope_694009_3_to_5_exon USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__host_taxon_id__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__host_taxon_name_lower__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__iedb_epitope_id__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__is_linear__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__mhc_allele__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__mhc_class_lower__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__product_lower__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__response_frequency_pos__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__sequence_aa_alternative__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__sequence_aa_original__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__start_aa_original__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__taxon_id__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__taxon_name_lower__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__variant_aa_length__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__variant_aa_type__idx
    ON public.epitope_694009_3_to_5_exon USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_cell_type__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_epi_start__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_epi_stop__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_is_linear__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_mhc_allele__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_product__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exon__virus_host_resp_freq__i
    ON public.epitope_694009_3_to_5_exon USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN endoribonuclease
CREATE MATERIALIZED VIEW public.epitope_694009_endoribonuc
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
             AND ann.product = 'endoribonuclease'
             AND epi.protein_name = 'endoribonuclease'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_endoribonuc
    OWNER TO geco;

CREATE INDEX epi_694009_endoribonuc__cell_type__idx
    ON public.epitope_694009_endoribonuc USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__epi_annotation_start__idx
    ON public.epitope_694009_endoribonuc USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__epi_annotation_stop__idx
    ON public.epitope_694009_endoribonuc USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__epi_frag_annotation_start__i
    ON public.epitope_694009_endoribonuc USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__epi_frag_annotation_stop__id
    ON public.epitope_694009_endoribonuc USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__host_taxon_id__idx
    ON public.epitope_694009_endoribonuc USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__host_taxon_name_lower__idx
    ON public.epitope_694009_endoribonuc USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__iedb_epitope_id__idx
    ON public.epitope_694009_endoribonuc USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__is_linear__idx
    ON public.epitope_694009_endoribonuc USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__mhc_allele__idx
    ON public.epitope_694009_endoribonuc USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__mhc_class_lower__idx
    ON public.epitope_694009_endoribonuc USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__product_lower__idx
    ON public.epitope_694009_endoribonuc USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__response_frequency_pos__idx
    ON public.epitope_694009_endoribonuc USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__sequence_aa_alternative__idx
    ON public.epitope_694009_endoribonuc USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__sequence_aa_original__idx
    ON public.epitope_694009_endoribonuc USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__start_aa_original__idx
    ON public.epitope_694009_endoribonuc USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__taxon_id__idx
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__taxon_name_lower__idx
    ON public.epitope_694009_endoribonuc USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__variant_aa_length__idx
    ON public.epitope_694009_endoribonuc USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__variant_aa_type__idx
    ON public.epitope_694009_endoribonuc USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_cell_type__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_epi_start__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_epi_stop__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_is_linear__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_mhc_allele__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_product__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuc__virus_host_resp_freq__i
    ON public.epitope_694009_endoribonuc USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN 2'-O-MTase
CREATE MATERIALIZED VIEW public.epitope_694009_2_o_mtase
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
             AND ann.product = '2'-O-MTase'
             AND epi.protein_name = '2'-O-MTase'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_2_o_mtase
    OWNER TO geco;

CREATE INDEX epi_694009_2_o_mtase__cell_type__idx
    ON public.epitope_694009_2_o_mtase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_annotation_start__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_annotation_stop__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_frag_annotation_start__i
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_frag_annotation_stop__id
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__host_taxon_id__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__host_taxon_name_lower__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__iedb_epitope_id__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__is_linear__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__mhc_allele__idx
    ON public.epitope_694009_2_o_mtase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__mhc_class_lower__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__product_lower__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__response_frequency_pos__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__sequence_aa_alternative__idx
    ON public.epitope_694009_2_o_mtase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__sequence_aa_original__idx
    ON public.epitope_694009_2_o_mtase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__start_aa_original__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__taxon_id__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__taxon_name_lower__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__variant_aa_length__idx
    ON public.epitope_694009_2_o_mtase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__variant_aa_type__idx
    ON public.epitope_694009_2_o_mtase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_cell_type__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_epi_start__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_epi_stop__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_is_linear__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_mhc_allele__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_product__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_resp_freq__i
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF1a polyprotein
CREATE MATERIALIZED VIEW public.epitope_694009_orf1a_polyp
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf1a_polyp
    OWNER TO geco;

CREATE INDEX epi_694009_orf1a_polyp__cell_type__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__epi_annotation_start__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__epi_annotation_stop__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__epi_frag_annotation_start__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf1a_polyp USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__host_taxon_id__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__host_taxon_name_lower__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__iedb_epitope_id__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__is_linear__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__mhc_allele__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__mhc_class_lower__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__product_lower__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__response_frequency_pos__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__sequence_aa_alternative__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__sequence_aa_original__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__start_aa_original__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__taxon_id__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__taxon_name_lower__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__variant_aa_length__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__variant_aa_type__idx
    ON public.epitope_694009_orf1a_polyp USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_cell_type__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_epi_start__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_epi_stop__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_is_linear__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_mhc_allele__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_product__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyp__virus_host_resp_freq__i
    ON public.epitope_694009_orf1a_polyp USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp1
CREATE MATERIALIZED VIEW public.epitope_694009_nsp1
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
             AND ann.product = 'nsp1'
             AND epi.protein_name = 'nsp1'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp1
    OWNER TO geco;

CREATE INDEX epi_694009_nsp1__cell_type__idx
    ON public.epitope_694009_nsp1 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_annotation_start__idx
    ON public.epitope_694009_nsp1 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_annotation_stop__idx
    ON public.epitope_694009_nsp1 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp1 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp1 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__host_taxon_id__idx
    ON public.epitope_694009_nsp1 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp1 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__iedb_epitope_id__idx
    ON public.epitope_694009_nsp1 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__is_linear__idx
    ON public.epitope_694009_nsp1 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__mhc_allele__idx
    ON public.epitope_694009_nsp1 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__mhc_class_lower__idx
    ON public.epitope_694009_nsp1 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__product_lower__idx
    ON public.epitope_694009_nsp1 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__response_frequency_pos__idx
    ON public.epitope_694009_nsp1 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp1 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__sequence_aa_original__idx
    ON public.epitope_694009_nsp1 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__start_aa_original__idx
    ON public.epitope_694009_nsp1 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__taxon_id__idx
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__taxon_name_lower__idx
    ON public.epitope_694009_nsp1 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__variant_aa_length__idx
    ON public.epitope_694009_nsp1 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__variant_aa_type__idx
    ON public.epitope_694009_nsp1 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_cell_type__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_epi_start__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_epi_stop__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_is_linear__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_product__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_resp_freq__i
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp2
CREATE MATERIALIZED VIEW public.epitope_694009_nsp2
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
             AND ann.product = 'nsp2'
             AND epi.protein_name = 'nsp2'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp2
    OWNER TO geco;

CREATE INDEX epi_694009_nsp2__cell_type__idx
    ON public.epitope_694009_nsp2 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_annotation_start__idx
    ON public.epitope_694009_nsp2 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_annotation_stop__idx
    ON public.epitope_694009_nsp2 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp2 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp2 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__host_taxon_id__idx
    ON public.epitope_694009_nsp2 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp2 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__iedb_epitope_id__idx
    ON public.epitope_694009_nsp2 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__is_linear__idx
    ON public.epitope_694009_nsp2 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__mhc_allele__idx
    ON public.epitope_694009_nsp2 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__mhc_class_lower__idx
    ON public.epitope_694009_nsp2 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__product_lower__idx
    ON public.epitope_694009_nsp2 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__response_frequency_pos__idx
    ON public.epitope_694009_nsp2 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp2 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__sequence_aa_original__idx
    ON public.epitope_694009_nsp2 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__start_aa_original__idx
    ON public.epitope_694009_nsp2 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__taxon_id__idx
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__taxon_name_lower__idx
    ON public.epitope_694009_nsp2 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__variant_aa_length__idx
    ON public.epitope_694009_nsp2 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__variant_aa_type__idx
    ON public.epitope_694009_nsp2 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_cell_type__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_epi_start__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_epi_stop__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_is_linear__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_product__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_resp_freq__i
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp3
CREATE MATERIALIZED VIEW public.epitope_694009_nsp3
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
             AND ann.product = 'nsp3'
             AND epi.protein_name = 'nsp3'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp3
    OWNER TO geco;

CREATE INDEX epi_694009_nsp3__cell_type__idx
    ON public.epitope_694009_nsp3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_annotation_start__idx
    ON public.epitope_694009_nsp3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_annotation_stop__idx
    ON public.epitope_694009_nsp3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__host_taxon_id__idx
    ON public.epitope_694009_nsp3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__iedb_epitope_id__idx
    ON public.epitope_694009_nsp3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__is_linear__idx
    ON public.epitope_694009_nsp3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__mhc_allele__idx
    ON public.epitope_694009_nsp3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__mhc_class_lower__idx
    ON public.epitope_694009_nsp3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__product_lower__idx
    ON public.epitope_694009_nsp3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__response_frequency_pos__idx
    ON public.epitope_694009_nsp3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__sequence_aa_original__idx
    ON public.epitope_694009_nsp3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__start_aa_original__idx
    ON public.epitope_694009_nsp3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__taxon_id__idx
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__taxon_name_lower__idx
    ON public.epitope_694009_nsp3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__variant_aa_length__idx
    ON public.epitope_694009_nsp3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__variant_aa_type__idx
    ON public.epitope_694009_nsp3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_cell_type__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_epi_start__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_epi_stop__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_is_linear__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_product__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_resp_freq__i
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp4
CREATE MATERIALIZED VIEW public.epitope_694009_nsp4
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
             AND ann.product = 'nsp4'
             AND epi.protein_name = 'nsp4'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp4
    OWNER TO geco;

CREATE INDEX epi_694009_nsp4__cell_type__idx
    ON public.epitope_694009_nsp4 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_annotation_start__idx
    ON public.epitope_694009_nsp4 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_annotation_stop__idx
    ON public.epitope_694009_nsp4 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp4 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp4 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__host_taxon_id__idx
    ON public.epitope_694009_nsp4 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp4 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__iedb_epitope_id__idx
    ON public.epitope_694009_nsp4 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__is_linear__idx
    ON public.epitope_694009_nsp4 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__mhc_allele__idx
    ON public.epitope_694009_nsp4 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__mhc_class_lower__idx
    ON public.epitope_694009_nsp4 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__product_lower__idx
    ON public.epitope_694009_nsp4 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__response_frequency_pos__idx
    ON public.epitope_694009_nsp4 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp4 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__sequence_aa_original__idx
    ON public.epitope_694009_nsp4 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__start_aa_original__idx
    ON public.epitope_694009_nsp4 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__taxon_id__idx
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__taxon_name_lower__idx
    ON public.epitope_694009_nsp4 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__variant_aa_length__idx
    ON public.epitope_694009_nsp4 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__variant_aa_type__idx
    ON public.epitope_694009_nsp4 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_cell_type__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_epi_start__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_epi_stop__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_is_linear__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_product__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_resp_freq__i
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN 3C-like protease
CREATE MATERIALIZED VIEW public.epitope_694009_3c_like_pro
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
             AND ann.product = '3C-like protease'
             AND epi.protein_name = '3C-like protease'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_3c_like_pro
    OWNER TO geco;

CREATE INDEX epi_694009_3c_like_pro__cell_type__idx
    ON public.epitope_694009_3c_like_pro USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__epi_annotation_start__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__epi_annotation_stop__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__epi_frag_annotation_start__i
    ON public.epitope_694009_3c_like_pro USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__epi_frag_annotation_stop__id
    ON public.epitope_694009_3c_like_pro USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__host_taxon_id__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__host_taxon_name_lower__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__iedb_epitope_id__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__is_linear__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__mhc_allele__idx
    ON public.epitope_694009_3c_like_pro USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__mhc_class_lower__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__product_lower__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__response_frequency_pos__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__sequence_aa_alternative__idx
    ON public.epitope_694009_3c_like_pro USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__sequence_aa_original__idx
    ON public.epitope_694009_3c_like_pro USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__start_aa_original__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__taxon_id__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__taxon_name_lower__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__variant_aa_length__idx
    ON public.epitope_694009_3c_like_pro USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__variant_aa_type__idx
    ON public.epitope_694009_3c_like_pro USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_cell_type__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_epi_start__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_epi_stop__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_is_linear__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_mhc_allele__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_product__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_pro__virus_host_resp_freq__i
    ON public.epitope_694009_3c_like_pro USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp6
CREATE MATERIALIZED VIEW public.epitope_694009_nsp6
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
             AND ann.product = 'nsp6'
             AND epi.protein_name = 'nsp6'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp6
    OWNER TO geco;

CREATE INDEX epi_694009_nsp6__cell_type__idx
    ON public.epitope_694009_nsp6 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_annotation_start__idx
    ON public.epitope_694009_nsp6 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_annotation_stop__idx
    ON public.epitope_694009_nsp6 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp6 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp6 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__host_taxon_id__idx
    ON public.epitope_694009_nsp6 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp6 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__iedb_epitope_id__idx
    ON public.epitope_694009_nsp6 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__is_linear__idx
    ON public.epitope_694009_nsp6 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__mhc_allele__idx
    ON public.epitope_694009_nsp6 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__mhc_class_lower__idx
    ON public.epitope_694009_nsp6 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__product_lower__idx
    ON public.epitope_694009_nsp6 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__response_frequency_pos__idx
    ON public.epitope_694009_nsp6 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp6 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__sequence_aa_original__idx
    ON public.epitope_694009_nsp6 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__start_aa_original__idx
    ON public.epitope_694009_nsp6 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__taxon_id__idx
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__taxon_name_lower__idx
    ON public.epitope_694009_nsp6 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__variant_aa_length__idx
    ON public.epitope_694009_nsp6 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__variant_aa_type__idx
    ON public.epitope_694009_nsp6 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_cell_type__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_epi_start__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_epi_stop__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_is_linear__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_product__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_resp_freq__i
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp7
CREATE MATERIALIZED VIEW public.epitope_694009_nsp7
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
             AND ann.product = 'nsp7'
             AND epi.protein_name = 'nsp7'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp7
    OWNER TO geco;

CREATE INDEX epi_694009_nsp7__cell_type__idx
    ON public.epitope_694009_nsp7 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_annotation_start__idx
    ON public.epitope_694009_nsp7 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_annotation_stop__idx
    ON public.epitope_694009_nsp7 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp7 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp7 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__host_taxon_id__idx
    ON public.epitope_694009_nsp7 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp7 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__iedb_epitope_id__idx
    ON public.epitope_694009_nsp7 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__is_linear__idx
    ON public.epitope_694009_nsp7 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__mhc_allele__idx
    ON public.epitope_694009_nsp7 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__mhc_class_lower__idx
    ON public.epitope_694009_nsp7 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__product_lower__idx
    ON public.epitope_694009_nsp7 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__response_frequency_pos__idx
    ON public.epitope_694009_nsp7 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp7 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__sequence_aa_original__idx
    ON public.epitope_694009_nsp7 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__start_aa_original__idx
    ON public.epitope_694009_nsp7 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__taxon_id__idx
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__taxon_name_lower__idx
    ON public.epitope_694009_nsp7 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__variant_aa_length__idx
    ON public.epitope_694009_nsp7 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__variant_aa_type__idx
    ON public.epitope_694009_nsp7 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_cell_type__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_epi_start__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_epi_stop__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_is_linear__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_product__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_resp_freq__i
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp8
CREATE MATERIALIZED VIEW public.epitope_694009_nsp8
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
             AND ann.product = 'nsp8'
             AND epi.protein_name = 'nsp8'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp8
    OWNER TO geco;

CREATE INDEX epi_694009_nsp8__cell_type__idx
    ON public.epitope_694009_nsp8 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_annotation_start__idx
    ON public.epitope_694009_nsp8 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_annotation_stop__idx
    ON public.epitope_694009_nsp8 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp8 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp8 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__host_taxon_id__idx
    ON public.epitope_694009_nsp8 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp8 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__iedb_epitope_id__idx
    ON public.epitope_694009_nsp8 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__is_linear__idx
    ON public.epitope_694009_nsp8 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__mhc_allele__idx
    ON public.epitope_694009_nsp8 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__mhc_class_lower__idx
    ON public.epitope_694009_nsp8 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__product_lower__idx
    ON public.epitope_694009_nsp8 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__response_frequency_pos__idx
    ON public.epitope_694009_nsp8 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp8 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__sequence_aa_original__idx
    ON public.epitope_694009_nsp8 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__start_aa_original__idx
    ON public.epitope_694009_nsp8 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__taxon_id__idx
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__taxon_name_lower__idx
    ON public.epitope_694009_nsp8 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__variant_aa_length__idx
    ON public.epitope_694009_nsp8 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__variant_aa_type__idx
    ON public.epitope_694009_nsp8 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_cell_type__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_epi_start__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_epi_stop__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_is_linear__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_product__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_resp_freq__i
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp9
CREATE MATERIALIZED VIEW public.epitope_694009_nsp9
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
             AND ann.product = 'nsp9'
             AND epi.protein_name = 'nsp9'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp9
    OWNER TO geco;

CREATE INDEX epi_694009_nsp9__cell_type__idx
    ON public.epitope_694009_nsp9 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_annotation_start__idx
    ON public.epitope_694009_nsp9 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_annotation_stop__idx
    ON public.epitope_694009_nsp9 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp9 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp9 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__host_taxon_id__idx
    ON public.epitope_694009_nsp9 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp9 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__iedb_epitope_id__idx
    ON public.epitope_694009_nsp9 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__is_linear__idx
    ON public.epitope_694009_nsp9 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__mhc_allele__idx
    ON public.epitope_694009_nsp9 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__mhc_class_lower__idx
    ON public.epitope_694009_nsp9 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__product_lower__idx
    ON public.epitope_694009_nsp9 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__response_frequency_pos__idx
    ON public.epitope_694009_nsp9 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp9 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__sequence_aa_original__idx
    ON public.epitope_694009_nsp9 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__start_aa_original__idx
    ON public.epitope_694009_nsp9 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__taxon_id__idx
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__taxon_name_lower__idx
    ON public.epitope_694009_nsp9 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__variant_aa_length__idx
    ON public.epitope_694009_nsp9 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__variant_aa_type__idx
    ON public.epitope_694009_nsp9 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_cell_type__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_epi_start__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_epi_stop__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_is_linear__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_product__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_resp_freq__i
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nsp10
CREATE MATERIALIZED VIEW public.epitope_694009_nsp10
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
             AND ann.product = 'nsp10'
             AND epi.protein_name = 'nsp10'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nsp10
    OWNER TO geco;

CREATE INDEX epi_694009_nsp10__cell_type__idx
    ON public.epitope_694009_nsp10 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_annotation_start__idx
    ON public.epitope_694009_nsp10 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_annotation_stop__idx
    ON public.epitope_694009_nsp10 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_frag_annotation_start__i
    ON public.epitope_694009_nsp10 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_frag_annotation_stop__id
    ON public.epitope_694009_nsp10 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__host_taxon_id__idx
    ON public.epitope_694009_nsp10 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__host_taxon_name_lower__idx
    ON public.epitope_694009_nsp10 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__iedb_epitope_id__idx
    ON public.epitope_694009_nsp10 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__is_linear__idx
    ON public.epitope_694009_nsp10 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__mhc_allele__idx
    ON public.epitope_694009_nsp10 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__mhc_class_lower__idx
    ON public.epitope_694009_nsp10 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__product_lower__idx
    ON public.epitope_694009_nsp10 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__response_frequency_pos__idx
    ON public.epitope_694009_nsp10 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__sequence_aa_alternative__idx
    ON public.epitope_694009_nsp10 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__sequence_aa_original__idx
    ON public.epitope_694009_nsp10 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__start_aa_original__idx
    ON public.epitope_694009_nsp10 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__taxon_id__idx
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__taxon_name_lower__idx
    ON public.epitope_694009_nsp10 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__variant_aa_length__idx
    ON public.epitope_694009_nsp10 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__variant_aa_type__idx
    ON public.epitope_694009_nsp10 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_cell_type__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_epi_start__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_epi_stop__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_is_linear__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_mhc_allele__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_product__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_resp_freq__i
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ndp11
CREATE MATERIALIZED VIEW public.epitope_694009_ndp11
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
             AND ann.product = 'ndp11'
             AND epi.protein_name = 'ndp11'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_ndp11
    OWNER TO geco;

CREATE INDEX epi_694009_ndp11__cell_type__idx
    ON public.epitope_694009_ndp11 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_annotation_start__idx
    ON public.epitope_694009_ndp11 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_annotation_stop__idx
    ON public.epitope_694009_ndp11 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_frag_annotation_start__i
    ON public.epitope_694009_ndp11 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_frag_annotation_stop__id
    ON public.epitope_694009_ndp11 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__host_taxon_id__idx
    ON public.epitope_694009_ndp11 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__host_taxon_name_lower__idx
    ON public.epitope_694009_ndp11 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__iedb_epitope_id__idx
    ON public.epitope_694009_ndp11 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__is_linear__idx
    ON public.epitope_694009_ndp11 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__mhc_allele__idx
    ON public.epitope_694009_ndp11 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__mhc_class_lower__idx
    ON public.epitope_694009_ndp11 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__product_lower__idx
    ON public.epitope_694009_ndp11 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__response_frequency_pos__idx
    ON public.epitope_694009_ndp11 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__sequence_aa_alternative__idx
    ON public.epitope_694009_ndp11 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__sequence_aa_original__idx
    ON public.epitope_694009_ndp11 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__start_aa_original__idx
    ON public.epitope_694009_ndp11 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__taxon_id__idx
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__taxon_name_lower__idx
    ON public.epitope_694009_ndp11 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__variant_aa_length__idx
    ON public.epitope_694009_ndp11 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__variant_aa_type__idx
    ON public.epitope_694009_ndp11 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_cell_type__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_epi_start__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_epi_stop__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_is_linear__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_mhc_allele__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_product__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_resp_freq__i
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN spike glycoprotein
CREATE MATERIALIZED VIEW public.epitope_694009_spike_glyco
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
             AND ann.product = 'spike glycoprotein'
             AND epi.protein_name = 'spike glycoprotein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_spike_glyco
    OWNER TO geco;

CREATE INDEX epi_694009_spike_glyco__cell_type__idx
    ON public.epitope_694009_spike_glyco USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__epi_annotation_start__idx
    ON public.epitope_694009_spike_glyco USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__epi_annotation_stop__idx
    ON public.epitope_694009_spike_glyco USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__epi_frag_annotation_start__i
    ON public.epitope_694009_spike_glyco USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__epi_frag_annotation_stop__id
    ON public.epitope_694009_spike_glyco USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__host_taxon_id__idx
    ON public.epitope_694009_spike_glyco USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__host_taxon_name_lower__idx
    ON public.epitope_694009_spike_glyco USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__iedb_epitope_id__idx
    ON public.epitope_694009_spike_glyco USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__is_linear__idx
    ON public.epitope_694009_spike_glyco USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__mhc_allele__idx
    ON public.epitope_694009_spike_glyco USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__mhc_class_lower__idx
    ON public.epitope_694009_spike_glyco USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__product_lower__idx
    ON public.epitope_694009_spike_glyco USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__response_frequency_pos__idx
    ON public.epitope_694009_spike_glyco USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__sequence_aa_alternative__idx
    ON public.epitope_694009_spike_glyco USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__sequence_aa_original__idx
    ON public.epitope_694009_spike_glyco USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__start_aa_original__idx
    ON public.epitope_694009_spike_glyco USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__taxon_id__idx
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__taxon_name_lower__idx
    ON public.epitope_694009_spike_glyco USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__variant_aa_length__idx
    ON public.epitope_694009_spike_glyco USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__variant_aa_type__idx
    ON public.epitope_694009_spike_glyco USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_cell_type__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_epi_start__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_epi_stop__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_is_linear__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_mhc_allele__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_product__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glyco__virus_host_resp_freq__i
    ON public.epitope_694009_spike_glyco USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF3a protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf3a_prote
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
             AND ann.product = 'ORF3a protein'
             AND epi.protein_name = 'ORF3a protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf3a_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf3a_prote__cell_type__idx
    ON public.epitope_694009_orf3a_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf3a_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf3a_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__host_taxon_id__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__is_linear__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__mhc_allele__idx
    ON public.epitope_694009_orf3a_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__product_lower__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf3a_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf3a_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__start_aa_original__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__taxon_id__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__variant_aa_length__idx
    ON public.epitope_694009_orf3a_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__variant_aa_type__idx
    ON public.epitope_694009_orf3a_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_product__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf3a_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF3b protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf3b_prote
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
             AND ann.product = 'ORF3b protein'
             AND epi.protein_name = 'ORF3b protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf3b_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf3b_prote__cell_type__idx
    ON public.epitope_694009_orf3b_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf3b_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf3b_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__host_taxon_id__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__is_linear__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__mhc_allele__idx
    ON public.epitope_694009_orf3b_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__product_lower__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf3b_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf3b_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__start_aa_original__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__taxon_id__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__variant_aa_length__idx
    ON public.epitope_694009_orf3b_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__variant_aa_type__idx
    ON public.epitope_694009_orf3b_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_product__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf3b_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN small envelope protein
CREATE MATERIALIZED VIEW public.epitope_694009_small_envel
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
             AND ann.product = 'small envelope protein'
             AND epi.protein_name = 'small envelope protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_small_envel
    OWNER TO geco;

CREATE INDEX epi_694009_small_envel__cell_type__idx
    ON public.epitope_694009_small_envel USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__epi_annotation_start__idx
    ON public.epitope_694009_small_envel USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__epi_annotation_stop__idx
    ON public.epitope_694009_small_envel USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__epi_frag_annotation_start__i
    ON public.epitope_694009_small_envel USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__epi_frag_annotation_stop__id
    ON public.epitope_694009_small_envel USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__host_taxon_id__idx
    ON public.epitope_694009_small_envel USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__host_taxon_name_lower__idx
    ON public.epitope_694009_small_envel USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__iedb_epitope_id__idx
    ON public.epitope_694009_small_envel USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__is_linear__idx
    ON public.epitope_694009_small_envel USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__mhc_allele__idx
    ON public.epitope_694009_small_envel USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__mhc_class_lower__idx
    ON public.epitope_694009_small_envel USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__product_lower__idx
    ON public.epitope_694009_small_envel USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__response_frequency_pos__idx
    ON public.epitope_694009_small_envel USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__sequence_aa_alternative__idx
    ON public.epitope_694009_small_envel USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__sequence_aa_original__idx
    ON public.epitope_694009_small_envel USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__start_aa_original__idx
    ON public.epitope_694009_small_envel USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__taxon_id__idx
    ON public.epitope_694009_small_envel USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__taxon_name_lower__idx
    ON public.epitope_694009_small_envel USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__variant_aa_length__idx
    ON public.epitope_694009_small_envel USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__variant_aa_type__idx
    ON public.epitope_694009_small_envel USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_cell_type__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_epi_start__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_epi_stop__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_is_linear__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_mhc_allele__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_product__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envel__virus_host_resp_freq__i
    ON public.epitope_694009_small_envel USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN membrane glycoprotein M
CREATE MATERIALIZED VIEW public.epitope_694009_membrane_gl
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_membrane_gl
    OWNER TO geco;

CREATE INDEX epi_694009_membrane_gl__cell_type__idx
    ON public.epitope_694009_membrane_gl USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__epi_annotation_start__idx
    ON public.epitope_694009_membrane_gl USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__epi_annotation_stop__idx
    ON public.epitope_694009_membrane_gl USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__epi_frag_annotation_start__i
    ON public.epitope_694009_membrane_gl USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__epi_frag_annotation_stop__id
    ON public.epitope_694009_membrane_gl USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__host_taxon_id__idx
    ON public.epitope_694009_membrane_gl USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__host_taxon_name_lower__idx
    ON public.epitope_694009_membrane_gl USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__iedb_epitope_id__idx
    ON public.epitope_694009_membrane_gl USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__is_linear__idx
    ON public.epitope_694009_membrane_gl USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__mhc_allele__idx
    ON public.epitope_694009_membrane_gl USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__mhc_class_lower__idx
    ON public.epitope_694009_membrane_gl USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__product_lower__idx
    ON public.epitope_694009_membrane_gl USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__response_frequency_pos__idx
    ON public.epitope_694009_membrane_gl USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__sequence_aa_alternative__idx
    ON public.epitope_694009_membrane_gl USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__sequence_aa_original__idx
    ON public.epitope_694009_membrane_gl USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__start_aa_original__idx
    ON public.epitope_694009_membrane_gl USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__taxon_id__idx
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__taxon_name_lower__idx
    ON public.epitope_694009_membrane_gl USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__variant_aa_length__idx
    ON public.epitope_694009_membrane_gl USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__variant_aa_type__idx
    ON public.epitope_694009_membrane_gl USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_cell_type__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_epi_start__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_epi_stop__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_is_linear__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_mhc_allele__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_product__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_gl__virus_host_resp_freq__i
    ON public.epitope_694009_membrane_gl USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF6 protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf6_protei
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
             AND ann.product = 'ORF6 protein'
             AND epi.protein_name = 'ORF6 protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf6_protei
    OWNER TO geco;

CREATE INDEX epi_694009_orf6_protei__cell_type__idx
    ON public.epitope_694009_orf6_protei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__epi_annotation_start__idx
    ON public.epitope_694009_orf6_protei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__epi_annotation_stop__idx
    ON public.epitope_694009_orf6_protei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__epi_frag_annotation_start__i
    ON public.epitope_694009_orf6_protei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf6_protei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__host_taxon_id__idx
    ON public.epitope_694009_orf6_protei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__host_taxon_name_lower__idx
    ON public.epitope_694009_orf6_protei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__iedb_epitope_id__idx
    ON public.epitope_694009_orf6_protei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__is_linear__idx
    ON public.epitope_694009_orf6_protei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__mhc_allele__idx
    ON public.epitope_694009_orf6_protei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__mhc_class_lower__idx
    ON public.epitope_694009_orf6_protei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__product_lower__idx
    ON public.epitope_694009_orf6_protei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__response_frequency_pos__idx
    ON public.epitope_694009_orf6_protei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__sequence_aa_alternative__idx
    ON public.epitope_694009_orf6_protei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__sequence_aa_original__idx
    ON public.epitope_694009_orf6_protei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__start_aa_original__idx
    ON public.epitope_694009_orf6_protei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__taxon_id__idx
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__taxon_name_lower__idx
    ON public.epitope_694009_orf6_protei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__variant_aa_length__idx
    ON public.epitope_694009_orf6_protei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__variant_aa_type__idx
    ON public.epitope_694009_orf6_protei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_cell_type__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_epi_start__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_epi_stop__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_is_linear__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_mhc_allele__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_product__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protei__virus_host_resp_freq__i
    ON public.epitope_694009_orf6_protei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF7a protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf7a_prote
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
             AND ann.product = 'ORF7a protein'
             AND epi.protein_name = 'ORF7a protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf7a_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf7a_prote__cell_type__idx
    ON public.epitope_694009_orf7a_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf7a_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf7a_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__host_taxon_id__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__is_linear__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__mhc_allele__idx
    ON public.epitope_694009_orf7a_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__product_lower__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf7a_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf7a_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__start_aa_original__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__taxon_id__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__variant_aa_length__idx
    ON public.epitope_694009_orf7a_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__variant_aa_type__idx
    ON public.epitope_694009_orf7a_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_product__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf7a_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF7b protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf7b_prote
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
             AND ann.product = 'ORF7b protein'
             AND epi.protein_name = 'ORF7b protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf7b_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf7b_prote__cell_type__idx
    ON public.epitope_694009_orf7b_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf7b_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf7b_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__host_taxon_id__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__is_linear__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__mhc_allele__idx
    ON public.epitope_694009_orf7b_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__product_lower__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf7b_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf7b_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__start_aa_original__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__taxon_id__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__variant_aa_length__idx
    ON public.epitope_694009_orf7b_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__variant_aa_type__idx
    ON public.epitope_694009_orf7b_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_product__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf7b_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF8a protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf8a_prote
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
             AND ann.product = 'ORF8a protein'
             AND epi.protein_name = 'ORF8a protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf8a_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf8a_prote__cell_type__idx
    ON public.epitope_694009_orf8a_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf8a_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf8a_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__host_taxon_id__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__is_linear__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__mhc_allele__idx
    ON public.epitope_694009_orf8a_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__product_lower__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf8a_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf8a_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__start_aa_original__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__taxon_id__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__variant_aa_length__idx
    ON public.epitope_694009_orf8a_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__variant_aa_type__idx
    ON public.epitope_694009_orf8a_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_product__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf8a_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF8b protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf8b_prote
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf8b_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf8b_prote__cell_type__idx
    ON public.epitope_694009_orf8b_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf8b_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf8b_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__host_taxon_id__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__is_linear__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__mhc_allele__idx
    ON public.epitope_694009_orf8b_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__product_lower__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf8b_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf8b_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__start_aa_original__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__taxon_id__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__variant_aa_length__idx
    ON public.epitope_694009_orf8b_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__variant_aa_type__idx
    ON public.epitope_694009_orf8b_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_product__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf8b_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN nucleocapsid protein
CREATE MATERIALIZED VIEW public.epitope_694009_nucleocapsi
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
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_nucleocapsi
    OWNER TO geco;

CREATE INDEX epi_694009_nucleocapsi__cell_type__idx
    ON public.epitope_694009_nucleocapsi USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__epi_annotation_start__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__epi_annotation_stop__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__epi_frag_annotation_start__i
    ON public.epitope_694009_nucleocapsi USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__epi_frag_annotation_stop__id
    ON public.epitope_694009_nucleocapsi USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__host_taxon_id__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__host_taxon_name_lower__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__iedb_epitope_id__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__is_linear__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__mhc_allele__idx
    ON public.epitope_694009_nucleocapsi USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__mhc_class_lower__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__product_lower__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__response_frequency_pos__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__sequence_aa_alternative__idx
    ON public.epitope_694009_nucleocapsi USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__sequence_aa_original__idx
    ON public.epitope_694009_nucleocapsi USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__start_aa_original__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__taxon_id__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__taxon_name_lower__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__variant_aa_length__idx
    ON public.epitope_694009_nucleocapsi USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__variant_aa_type__idx
    ON public.epitope_694009_nucleocapsi USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_cell_type__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_epi_start__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_epi_stop__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_is_linear__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_mhc_allele__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_product__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsi__virus_host_resp_freq__i
    ON public.epitope_694009_nucleocapsi USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF9b protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf9b_prote
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
             AND ann.product = 'ORF9b protein'
             AND epi.protein_name = 'ORF9b protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf9b_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf9b_prote__cell_type__idx
    ON public.epitope_694009_orf9b_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf9b_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf9b_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__host_taxon_id__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__is_linear__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__mhc_allele__idx
    ON public.epitope_694009_orf9b_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__product_lower__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf9b_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf9b_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__start_aa_original__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__taxon_id__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__variant_aa_length__idx
    ON public.epitope_694009_orf9b_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__variant_aa_type__idx
    ON public.epitope_694009_orf9b_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_product__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf9b_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 694009 AND PROTEIN ORF9a protein
CREATE MATERIALIZED VIEW public.epitope_694009_orf9a_prote
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
             AND ann.product = 'ORF9a protein'
             AND epi.protein_name = 'ORF9a protein'
             AND vir.taxon_id = 694009)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_694009_orf9a_prote
    OWNER TO geco;

CREATE INDEX epi_694009_orf9a_prote__cell_type__idx
    ON public.epitope_694009_orf9a_prote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__epi_annotation_start__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__epi_annotation_stop__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__epi_frag_annotation_start__i
    ON public.epitope_694009_orf9a_prote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__epi_frag_annotation_stop__id
    ON public.epitope_694009_orf9a_prote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__host_taxon_id__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__host_taxon_name_lower__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__iedb_epitope_id__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__is_linear__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__mhc_allele__idx
    ON public.epitope_694009_orf9a_prote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__mhc_class_lower__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__product_lower__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__response_frequency_pos__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__sequence_aa_alternative__idx
    ON public.epitope_694009_orf9a_prote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__sequence_aa_original__idx
    ON public.epitope_694009_orf9a_prote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__start_aa_original__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__taxon_id__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__taxon_name_lower__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__variant_aa_length__idx
    ON public.epitope_694009_orf9a_prote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__variant_aa_type__idx
    ON public.epitope_694009_orf9a_prote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_cell_type__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_epi_start__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_epi_stop__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_is_linear__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_mhc_allele__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_product__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_prote__virus_host_resp_freq__i
    ON public.epitope_694009_orf9a_prote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


