CREATE MATERIALIZED VIEW public.epitope_$virus_id_$short_prot_name
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
             AND ann.product = '$prot_name'
             AND epi.protein_name = '$prot_name'
             AND vir.taxon_id = $virus_id)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_$virus_id_$short_prot_name
    OWNER TO geco;

CREATE INDEX epi_$virus_id_$short_prot_name__cell_type__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_annotation_start__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_annotation_stop__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_frag_annotation_start__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_frag_annotation_stop__id
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__host_taxon_id__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__host_taxon_name_lower__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__iedb_epitope_id__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__is_linear__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__mhc_allele__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__mhc_class_lower__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__product_lower__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__response_frequency_pos__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__sequence_aa_alternative__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__sequence_aa_original__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__start_aa_original__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__taxon_id__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__taxon_name_lower__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__variant_aa_length__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__variant_aa_type__idx
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_taxon_and_host_taxon_id__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_cell_type__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_epi_start__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_epi_stop__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_is_linear__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_mhc_allele__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_product__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_resp_freq__i
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
