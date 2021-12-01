-- $virus_id can be replaced with the virus taxon id, while $short_prot_name can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_$virus_id_$short_prot_name (
    iedb_epitope_id int4 NULL,
    epitope_iri varchar NULL,
    cell_type varchar NULL,
    mhc_class varchar NULL,
    mhc_allele varchar NULL,
    response_frequency_pos float4 NULL,
    epi_annotation_start int4 NULL,
    epi_annotation_stop int4 NULL,
    is_linear bool NULL,
    assay_type varchar NULL,
    epi_fragment_sequence varchar NULL,
    epi_frag_annotation_start int4 NULL,
    epi_frag_annotation_stop int4 NULL,
    taxon_id int4 NULL,
    taxon_name varchar NULL,
    host_taxon_id int4 NULL,
    host_taxon_name varchar NULL,
    sequence_id int4 NULL,
    product varchar NULL,
    aminoacid_variant_id int4 NULL,
    start_aa_original int4 NULL,
    sequence_aa_original varchar NOT NULL,
    sequence_aa_alternative varchar NOT NULL,
    variant_aa_length int4 NOT NULL,
    variant_aa_type varchar NOT NULL
);

ALTER TABLE public.epitope_$virus_id_$short_prot_name
    OWNER TO geco;

CREATE INDEX epi_$virus_id_$short_prot_name__cell_type
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_an_start
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_an_nstop
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_frag_an_start
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__epi_frag_an_stop
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__host_tax_id
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__host_tax_name
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__iedb_id
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__is_linear
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__mhc_allele
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__mhc_class_lower
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__product_lower
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__response_freq_pos
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__seq_aa_alt
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__seq_aa_orig
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__start_aa_orig
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__taxon_id
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__taxon_name_lower
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__variant_aa_length
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__variant_aa_type
    ON public.epitope_$virus_id_$short_prot_name USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_n_host_tax_id
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_host_cell_type
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_host_epi_start
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_epi_stop
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__virus_host_is_linear
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_host_mhc_allele
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_host_product
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_$virus_id_$short_prot_name__vir_host_resp_freq
    ON public.epitope_$virus_id_$short_prot_name USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
