-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 1AB polyprotein
-- 1335626 can be replaced with the virus taxon id, while 1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_1ab_polyprotein (
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

ALTER TABLE public.epitope_1335626_1ab_polyprotein
    OWNER TO geco;

CREATE INDEX epi_1335626_1ab_polyprotein__cell_type
    ON public.epitope_1335626_1ab_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__epi_an_start
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__epi_an_nstop
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__epi_frag_an_start
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__epi_frag_an_stop
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__host_tax_id
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__host_tax_name
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__iedb_id
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__is_linear
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__mhc_allele
    ON public.epitope_1335626_1ab_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__mhc_class_lower
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__product_lower
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__response_freq_pos
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__seq_aa_alt
    ON public.epitope_1335626_1ab_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__seq_aa_orig
    ON public.epitope_1335626_1ab_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__start_aa_orig
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__taxon_id
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__taxon_name_lower
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__variant_aa_length
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__variant_aa_type
    ON public.epitope_1335626_1ab_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_n_host_tax_id
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_host_cell_type
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_host_epi_start
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__virus_host_epi_stop
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__virus_host_is_linear
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_host_mhc_allele
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_host_product
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1ab_polyprotein__vir_host_resp_freq
    ON public.epitope_1335626_1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT RNA-dependent RNA polymerase
-- 1335626 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_rna_dependent_rna_polymerase (
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

ALTER TABLE public.epitope_1335626_rna_dependent_rna_polymerase
    OWNER TO geco;

CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__cell_type
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__epi_an_start
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__epi_an_nstop
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__epi_frag_an_start
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__epi_frag_an_stop
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__host_tax_id
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__host_tax_name
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__iedb_id
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__is_linear
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__mhc_allele
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__mhc_class_lower
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__product_lower
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__response_freq_pos
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__seq_aa_alt
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__seq_aa_orig
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__start_aa_orig
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__taxon_id
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__taxon_name_lower
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__variant_aa_length
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__variant_aa_type
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_n_host_tax_id
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_host_cell_type
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_host_epi_start
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__virus_host_epi_stop
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__virus_host_is_linear
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_host_mhc_allele
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_host_product
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_rna_dependent_rna_polymerase__vir_host_resp_freq
    ON public.epitope_1335626_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT Hel
-- 1335626 can be replaced with the virus taxon id, while hel can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_hel (
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

ALTER TABLE public.epitope_1335626_hel
    OWNER TO geco;

CREATE INDEX epi_1335626_hel__cell_type
    ON public.epitope_1335626_hel USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_an_start
    ON public.epitope_1335626_hel USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_an_nstop
    ON public.epitope_1335626_hel USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_frag_an_start
    ON public.epitope_1335626_hel USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__epi_frag_an_stop
    ON public.epitope_1335626_hel USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__host_tax_id
    ON public.epitope_1335626_hel USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__host_tax_name
    ON public.epitope_1335626_hel USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__iedb_id
    ON public.epitope_1335626_hel USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__is_linear
    ON public.epitope_1335626_hel USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__mhc_allele
    ON public.epitope_1335626_hel USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__mhc_class_lower
    ON public.epitope_1335626_hel USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__product_lower
    ON public.epitope_1335626_hel USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__response_freq_pos
    ON public.epitope_1335626_hel USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__seq_aa_alt
    ON public.epitope_1335626_hel USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__seq_aa_orig
    ON public.epitope_1335626_hel USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__start_aa_orig
    ON public.epitope_1335626_hel USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__taxon_id
    ON public.epitope_1335626_hel USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__taxon_name_lower
    ON public.epitope_1335626_hel USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__variant_aa_length
    ON public.epitope_1335626_hel USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__variant_aa_type
    ON public.epitope_1335626_hel USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_n_host_tax_id
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_host_cell_type
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_host_epi_start
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_epi_stop
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__virus_host_is_linear
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_host_mhc_allele
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_host_product
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_hel__vir_host_resp_freq
    ON public.epitope_1335626_hel USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT ExoN
-- 1335626 can be replaced with the virus taxon id, while exon can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_exon (
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

ALTER TABLE public.epitope_1335626_exon
    OWNER TO geco;

CREATE INDEX epi_1335626_exon__cell_type
    ON public.epitope_1335626_exon USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_an_start
    ON public.epitope_1335626_exon USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_an_nstop
    ON public.epitope_1335626_exon USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_frag_an_start
    ON public.epitope_1335626_exon USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__epi_frag_an_stop
    ON public.epitope_1335626_exon USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__host_tax_id
    ON public.epitope_1335626_exon USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__host_tax_name
    ON public.epitope_1335626_exon USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__iedb_id
    ON public.epitope_1335626_exon USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__is_linear
    ON public.epitope_1335626_exon USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__mhc_allele
    ON public.epitope_1335626_exon USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__mhc_class_lower
    ON public.epitope_1335626_exon USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__product_lower
    ON public.epitope_1335626_exon USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__response_freq_pos
    ON public.epitope_1335626_exon USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__seq_aa_alt
    ON public.epitope_1335626_exon USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__seq_aa_orig
    ON public.epitope_1335626_exon USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__start_aa_orig
    ON public.epitope_1335626_exon USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__taxon_id
    ON public.epitope_1335626_exon USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__taxon_name_lower
    ON public.epitope_1335626_exon USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__variant_aa_length
    ON public.epitope_1335626_exon USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__variant_aa_type
    ON public.epitope_1335626_exon USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_n_host_tax_id
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_host_cell_type
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_host_epi_start
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_epi_stop
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__virus_host_is_linear
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_host_mhc_allele
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_host_product
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_exon__vir_host_resp_freq
    ON public.epitope_1335626_exon USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NendoU
-- 1335626 can be replaced with the virus taxon id, while nendou can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nendou (
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

ALTER TABLE public.epitope_1335626_nendou
    OWNER TO geco;

CREATE INDEX epi_1335626_nendou__cell_type
    ON public.epitope_1335626_nendou USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_an_start
    ON public.epitope_1335626_nendou USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_an_nstop
    ON public.epitope_1335626_nendou USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_frag_an_start
    ON public.epitope_1335626_nendou USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__epi_frag_an_stop
    ON public.epitope_1335626_nendou USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__host_tax_id
    ON public.epitope_1335626_nendou USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__host_tax_name
    ON public.epitope_1335626_nendou USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__iedb_id
    ON public.epitope_1335626_nendou USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__is_linear
    ON public.epitope_1335626_nendou USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__mhc_allele
    ON public.epitope_1335626_nendou USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__mhc_class_lower
    ON public.epitope_1335626_nendou USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__product_lower
    ON public.epitope_1335626_nendou USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__response_freq_pos
    ON public.epitope_1335626_nendou USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__seq_aa_alt
    ON public.epitope_1335626_nendou USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__seq_aa_orig
    ON public.epitope_1335626_nendou USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__start_aa_orig
    ON public.epitope_1335626_nendou USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__taxon_id
    ON public.epitope_1335626_nendou USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__taxon_name_lower
    ON public.epitope_1335626_nendou USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__variant_aa_length
    ON public.epitope_1335626_nendou USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__variant_aa_type
    ON public.epitope_1335626_nendou USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_n_host_tax_id
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_host_cell_type
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_host_epi_start
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_epi_stop
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__virus_host_is_linear
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_host_mhc_allele
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_host_product
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nendou__vir_host_resp_freq
    ON public.epitope_1335626_nendou USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 2'-O-methyltransferase
-- 1335626 can be replaced with the virus taxon id, while 2_o_methyltransferase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_2_o_methyltransferase (
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

ALTER TABLE public.epitope_1335626_2_o_methyltransferase
    OWNER TO geco;

CREATE INDEX epi_1335626_2_o_methyltransferase__cell_type
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__epi_an_start
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__epi_an_nstop
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__epi_frag_an_start
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__epi_frag_an_stop
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__host_tax_id
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__host_tax_name
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__iedb_id
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__is_linear
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__mhc_allele
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__mhc_class_lower
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__product_lower
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__response_freq_pos
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__seq_aa_alt
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__seq_aa_orig
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__start_aa_orig
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__taxon_id
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__taxon_name_lower
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__variant_aa_length
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__variant_aa_type
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_n_host_tax_id
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_host_cell_type
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_host_epi_start
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__virus_host_epi_stop
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__virus_host_is_linear
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_host_mhc_allele
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_host_product
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_2_o_methyltransferase__vir_host_resp_freq
    ON public.epitope_1335626_2_o_methyltransferase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT 1A polyprotein
-- 1335626 can be replaced with the virus taxon id, while 1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_1a_polyprotein (
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

ALTER TABLE public.epitope_1335626_1a_polyprotein
    OWNER TO geco;

CREATE INDEX epi_1335626_1a_polyprotein__cell_type
    ON public.epitope_1335626_1a_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__epi_an_start
    ON public.epitope_1335626_1a_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__epi_an_nstop
    ON public.epitope_1335626_1a_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__epi_frag_an_start
    ON public.epitope_1335626_1a_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__epi_frag_an_stop
    ON public.epitope_1335626_1a_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__host_tax_id
    ON public.epitope_1335626_1a_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__host_tax_name
    ON public.epitope_1335626_1a_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__iedb_id
    ON public.epitope_1335626_1a_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__is_linear
    ON public.epitope_1335626_1a_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__mhc_allele
    ON public.epitope_1335626_1a_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__mhc_class_lower
    ON public.epitope_1335626_1a_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__product_lower
    ON public.epitope_1335626_1a_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__response_freq_pos
    ON public.epitope_1335626_1a_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__seq_aa_alt
    ON public.epitope_1335626_1a_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__seq_aa_orig
    ON public.epitope_1335626_1a_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__start_aa_orig
    ON public.epitope_1335626_1a_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__taxon_id
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__taxon_name_lower
    ON public.epitope_1335626_1a_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__variant_aa_length
    ON public.epitope_1335626_1a_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__variant_aa_type
    ON public.epitope_1335626_1a_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_n_host_tax_id
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_host_cell_type
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_host_epi_start
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__virus_host_epi_stop
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__virus_host_is_linear
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_host_mhc_allele
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_host_product
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_1a_polyprotein__vir_host_resp_freq
    ON public.epitope_1335626_1a_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp1 protein
-- 1335626 can be replaced with the virus taxon id, while nsp1_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp1_protein (
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

ALTER TABLE public.epitope_1335626_nsp1_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp1_protein__cell_type
    ON public.epitope_1335626_nsp1_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__epi_an_start
    ON public.epitope_1335626_nsp1_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__epi_an_nstop
    ON public.epitope_1335626_nsp1_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp1_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp1_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__host_tax_id
    ON public.epitope_1335626_nsp1_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__host_tax_name
    ON public.epitope_1335626_nsp1_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__iedb_id
    ON public.epitope_1335626_nsp1_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__is_linear
    ON public.epitope_1335626_nsp1_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__mhc_allele
    ON public.epitope_1335626_nsp1_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__mhc_class_lower
    ON public.epitope_1335626_nsp1_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__product_lower
    ON public.epitope_1335626_nsp1_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__response_freq_pos
    ON public.epitope_1335626_nsp1_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__seq_aa_alt
    ON public.epitope_1335626_nsp1_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__seq_aa_orig
    ON public.epitope_1335626_nsp1_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__start_aa_orig
    ON public.epitope_1335626_nsp1_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__taxon_id
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__taxon_name_lower
    ON public.epitope_1335626_nsp1_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__variant_aa_length
    ON public.epitope_1335626_nsp1_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__variant_aa_type
    ON public.epitope_1335626_nsp1_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_host_product
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp1_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp1_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp2 protein
-- 1335626 can be replaced with the virus taxon id, while nsp2_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp2_protein (
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

ALTER TABLE public.epitope_1335626_nsp2_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp2_protein__cell_type
    ON public.epitope_1335626_nsp2_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__epi_an_start
    ON public.epitope_1335626_nsp2_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__epi_an_nstop
    ON public.epitope_1335626_nsp2_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp2_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp2_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__host_tax_id
    ON public.epitope_1335626_nsp2_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__host_tax_name
    ON public.epitope_1335626_nsp2_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__iedb_id
    ON public.epitope_1335626_nsp2_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__is_linear
    ON public.epitope_1335626_nsp2_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__mhc_allele
    ON public.epitope_1335626_nsp2_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__mhc_class_lower
    ON public.epitope_1335626_nsp2_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__product_lower
    ON public.epitope_1335626_nsp2_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__response_freq_pos
    ON public.epitope_1335626_nsp2_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__seq_aa_alt
    ON public.epitope_1335626_nsp2_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__seq_aa_orig
    ON public.epitope_1335626_nsp2_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__start_aa_orig
    ON public.epitope_1335626_nsp2_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__taxon_id
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__taxon_name_lower
    ON public.epitope_1335626_nsp2_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__variant_aa_length
    ON public.epitope_1335626_nsp2_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__variant_aa_type
    ON public.epitope_1335626_nsp2_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_host_product
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp2_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp2_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp3 protein
-- 1335626 can be replaced with the virus taxon id, while nsp3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp3_protein (
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

ALTER TABLE public.epitope_1335626_nsp3_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp3_protein__cell_type
    ON public.epitope_1335626_nsp3_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__epi_an_start
    ON public.epitope_1335626_nsp3_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__epi_an_nstop
    ON public.epitope_1335626_nsp3_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp3_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp3_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__host_tax_id
    ON public.epitope_1335626_nsp3_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__host_tax_name
    ON public.epitope_1335626_nsp3_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__iedb_id
    ON public.epitope_1335626_nsp3_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__is_linear
    ON public.epitope_1335626_nsp3_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__mhc_allele
    ON public.epitope_1335626_nsp3_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__mhc_class_lower
    ON public.epitope_1335626_nsp3_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__product_lower
    ON public.epitope_1335626_nsp3_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__response_freq_pos
    ON public.epitope_1335626_nsp3_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__seq_aa_alt
    ON public.epitope_1335626_nsp3_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__seq_aa_orig
    ON public.epitope_1335626_nsp3_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__start_aa_orig
    ON public.epitope_1335626_nsp3_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__taxon_id
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__taxon_name_lower
    ON public.epitope_1335626_nsp3_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__variant_aa_length
    ON public.epitope_1335626_nsp3_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__variant_aa_type
    ON public.epitope_1335626_nsp3_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_host_product
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp3_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp3_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp4 protein
-- 1335626 can be replaced with the virus taxon id, while nsp4_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp4_protein (
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

ALTER TABLE public.epitope_1335626_nsp4_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp4_protein__cell_type
    ON public.epitope_1335626_nsp4_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__epi_an_start
    ON public.epitope_1335626_nsp4_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__epi_an_nstop
    ON public.epitope_1335626_nsp4_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp4_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp4_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__host_tax_id
    ON public.epitope_1335626_nsp4_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__host_tax_name
    ON public.epitope_1335626_nsp4_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__iedb_id
    ON public.epitope_1335626_nsp4_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__is_linear
    ON public.epitope_1335626_nsp4_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__mhc_allele
    ON public.epitope_1335626_nsp4_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__mhc_class_lower
    ON public.epitope_1335626_nsp4_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__product_lower
    ON public.epitope_1335626_nsp4_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__response_freq_pos
    ON public.epitope_1335626_nsp4_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__seq_aa_alt
    ON public.epitope_1335626_nsp4_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__seq_aa_orig
    ON public.epitope_1335626_nsp4_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__start_aa_orig
    ON public.epitope_1335626_nsp4_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__taxon_id
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__taxon_name_lower
    ON public.epitope_1335626_nsp4_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__variant_aa_length
    ON public.epitope_1335626_nsp4_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__variant_aa_type
    ON public.epitope_1335626_nsp4_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_host_product
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp4_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp4_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp5 protein
-- 1335626 can be replaced with the virus taxon id, while nsp5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp5_protein (
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

ALTER TABLE public.epitope_1335626_nsp5_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp5_protein__cell_type
    ON public.epitope_1335626_nsp5_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__epi_an_start
    ON public.epitope_1335626_nsp5_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__epi_an_nstop
    ON public.epitope_1335626_nsp5_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp5_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp5_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__host_tax_id
    ON public.epitope_1335626_nsp5_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__host_tax_name
    ON public.epitope_1335626_nsp5_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__iedb_id
    ON public.epitope_1335626_nsp5_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__is_linear
    ON public.epitope_1335626_nsp5_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__mhc_allele
    ON public.epitope_1335626_nsp5_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__mhc_class_lower
    ON public.epitope_1335626_nsp5_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__product_lower
    ON public.epitope_1335626_nsp5_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__response_freq_pos
    ON public.epitope_1335626_nsp5_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__seq_aa_alt
    ON public.epitope_1335626_nsp5_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__seq_aa_orig
    ON public.epitope_1335626_nsp5_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__start_aa_orig
    ON public.epitope_1335626_nsp5_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__taxon_id
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__taxon_name_lower
    ON public.epitope_1335626_nsp5_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__variant_aa_length
    ON public.epitope_1335626_nsp5_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__variant_aa_type
    ON public.epitope_1335626_nsp5_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_host_product
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp5_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp5_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp6 protein
-- 1335626 can be replaced with the virus taxon id, while nsp6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp6_protein (
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

ALTER TABLE public.epitope_1335626_nsp6_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp6_protein__cell_type
    ON public.epitope_1335626_nsp6_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__epi_an_start
    ON public.epitope_1335626_nsp6_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__epi_an_nstop
    ON public.epitope_1335626_nsp6_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp6_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp6_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__host_tax_id
    ON public.epitope_1335626_nsp6_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__host_tax_name
    ON public.epitope_1335626_nsp6_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__iedb_id
    ON public.epitope_1335626_nsp6_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__is_linear
    ON public.epitope_1335626_nsp6_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__mhc_allele
    ON public.epitope_1335626_nsp6_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__mhc_class_lower
    ON public.epitope_1335626_nsp6_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__product_lower
    ON public.epitope_1335626_nsp6_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__response_freq_pos
    ON public.epitope_1335626_nsp6_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__seq_aa_alt
    ON public.epitope_1335626_nsp6_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__seq_aa_orig
    ON public.epitope_1335626_nsp6_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__start_aa_orig
    ON public.epitope_1335626_nsp6_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__taxon_id
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__taxon_name_lower
    ON public.epitope_1335626_nsp6_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__variant_aa_length
    ON public.epitope_1335626_nsp6_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__variant_aa_type
    ON public.epitope_1335626_nsp6_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_host_product
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp6_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp6_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp7 protein
-- 1335626 can be replaced with the virus taxon id, while nsp7_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp7_protein (
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

ALTER TABLE public.epitope_1335626_nsp7_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp7_protein__cell_type
    ON public.epitope_1335626_nsp7_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__epi_an_start
    ON public.epitope_1335626_nsp7_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__epi_an_nstop
    ON public.epitope_1335626_nsp7_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp7_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp7_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__host_tax_id
    ON public.epitope_1335626_nsp7_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__host_tax_name
    ON public.epitope_1335626_nsp7_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__iedb_id
    ON public.epitope_1335626_nsp7_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__is_linear
    ON public.epitope_1335626_nsp7_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__mhc_allele
    ON public.epitope_1335626_nsp7_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__mhc_class_lower
    ON public.epitope_1335626_nsp7_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__product_lower
    ON public.epitope_1335626_nsp7_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__response_freq_pos
    ON public.epitope_1335626_nsp7_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__seq_aa_alt
    ON public.epitope_1335626_nsp7_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__seq_aa_orig
    ON public.epitope_1335626_nsp7_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__start_aa_orig
    ON public.epitope_1335626_nsp7_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__taxon_id
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__taxon_name_lower
    ON public.epitope_1335626_nsp7_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__variant_aa_length
    ON public.epitope_1335626_nsp7_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__variant_aa_type
    ON public.epitope_1335626_nsp7_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_host_product
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp7_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp7_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp8 protein
-- 1335626 can be replaced with the virus taxon id, while nsp8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp8_protein (
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

ALTER TABLE public.epitope_1335626_nsp8_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp8_protein__cell_type
    ON public.epitope_1335626_nsp8_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__epi_an_start
    ON public.epitope_1335626_nsp8_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__epi_an_nstop
    ON public.epitope_1335626_nsp8_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp8_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp8_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__host_tax_id
    ON public.epitope_1335626_nsp8_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__host_tax_name
    ON public.epitope_1335626_nsp8_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__iedb_id
    ON public.epitope_1335626_nsp8_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__is_linear
    ON public.epitope_1335626_nsp8_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__mhc_allele
    ON public.epitope_1335626_nsp8_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__mhc_class_lower
    ON public.epitope_1335626_nsp8_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__product_lower
    ON public.epitope_1335626_nsp8_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__response_freq_pos
    ON public.epitope_1335626_nsp8_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__seq_aa_alt
    ON public.epitope_1335626_nsp8_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__seq_aa_orig
    ON public.epitope_1335626_nsp8_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__start_aa_orig
    ON public.epitope_1335626_nsp8_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__taxon_id
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__taxon_name_lower
    ON public.epitope_1335626_nsp8_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__variant_aa_length
    ON public.epitope_1335626_nsp8_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__variant_aa_type
    ON public.epitope_1335626_nsp8_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_host_product
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp8_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp8_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp9 protein
-- 1335626 can be replaced with the virus taxon id, while nsp9_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp9_protein (
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

ALTER TABLE public.epitope_1335626_nsp9_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp9_protein__cell_type
    ON public.epitope_1335626_nsp9_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__epi_an_start
    ON public.epitope_1335626_nsp9_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__epi_an_nstop
    ON public.epitope_1335626_nsp9_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp9_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp9_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__host_tax_id
    ON public.epitope_1335626_nsp9_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__host_tax_name
    ON public.epitope_1335626_nsp9_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__iedb_id
    ON public.epitope_1335626_nsp9_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__is_linear
    ON public.epitope_1335626_nsp9_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__mhc_allele
    ON public.epitope_1335626_nsp9_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__mhc_class_lower
    ON public.epitope_1335626_nsp9_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__product_lower
    ON public.epitope_1335626_nsp9_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__response_freq_pos
    ON public.epitope_1335626_nsp9_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__seq_aa_alt
    ON public.epitope_1335626_nsp9_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__seq_aa_orig
    ON public.epitope_1335626_nsp9_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__start_aa_orig
    ON public.epitope_1335626_nsp9_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__taxon_id
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__taxon_name_lower
    ON public.epitope_1335626_nsp9_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__variant_aa_length
    ON public.epitope_1335626_nsp9_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__variant_aa_type
    ON public.epitope_1335626_nsp9_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_host_product
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp9_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp9_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp10 protein
-- 1335626 can be replaced with the virus taxon id, while nsp10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp10_protein (
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

ALTER TABLE public.epitope_1335626_nsp10_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp10_protein__cell_type
    ON public.epitope_1335626_nsp10_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__epi_an_start
    ON public.epitope_1335626_nsp10_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__epi_an_nstop
    ON public.epitope_1335626_nsp10_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp10_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp10_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__host_tax_id
    ON public.epitope_1335626_nsp10_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__host_tax_name
    ON public.epitope_1335626_nsp10_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__iedb_id
    ON public.epitope_1335626_nsp10_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__is_linear
    ON public.epitope_1335626_nsp10_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__mhc_allele
    ON public.epitope_1335626_nsp10_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__mhc_class_lower
    ON public.epitope_1335626_nsp10_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__product_lower
    ON public.epitope_1335626_nsp10_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__response_freq_pos
    ON public.epitope_1335626_nsp10_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__seq_aa_alt
    ON public.epitope_1335626_nsp10_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__seq_aa_orig
    ON public.epitope_1335626_nsp10_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__start_aa_orig
    ON public.epitope_1335626_nsp10_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__taxon_id
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__taxon_name_lower
    ON public.epitope_1335626_nsp10_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__variant_aa_length
    ON public.epitope_1335626_nsp10_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__variant_aa_type
    ON public.epitope_1335626_nsp10_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_host_product
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp10_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp10_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nsp11 protein
-- 1335626 can be replaced with the virus taxon id, while nsp11_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nsp11_protein (
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

ALTER TABLE public.epitope_1335626_nsp11_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nsp11_protein__cell_type
    ON public.epitope_1335626_nsp11_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__epi_an_start
    ON public.epitope_1335626_nsp11_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__epi_an_nstop
    ON public.epitope_1335626_nsp11_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__epi_frag_an_start
    ON public.epitope_1335626_nsp11_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__epi_frag_an_stop
    ON public.epitope_1335626_nsp11_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__host_tax_id
    ON public.epitope_1335626_nsp11_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__host_tax_name
    ON public.epitope_1335626_nsp11_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__iedb_id
    ON public.epitope_1335626_nsp11_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__is_linear
    ON public.epitope_1335626_nsp11_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__mhc_allele
    ON public.epitope_1335626_nsp11_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__mhc_class_lower
    ON public.epitope_1335626_nsp11_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__product_lower
    ON public.epitope_1335626_nsp11_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__response_freq_pos
    ON public.epitope_1335626_nsp11_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__seq_aa_alt
    ON public.epitope_1335626_nsp11_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__seq_aa_orig
    ON public.epitope_1335626_nsp11_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__start_aa_orig
    ON public.epitope_1335626_nsp11_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__taxon_id
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__taxon_name_lower
    ON public.epitope_1335626_nsp11_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__variant_aa_length
    ON public.epitope_1335626_nsp11_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__variant_aa_type
    ON public.epitope_1335626_nsp11_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_host_cell_type
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_host_epi_start
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__virus_host_epi_stop
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__virus_host_is_linear
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_host_product
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nsp11_protein__vir_host_resp_freq
    ON public.epitope_1335626_nsp11_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT spike protein
-- 1335626 can be replaced with the virus taxon id, while spike_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_spike_protein (
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

ALTER TABLE public.epitope_1335626_spike_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_spike_protein__cell_type
    ON public.epitope_1335626_spike_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__epi_an_start
    ON public.epitope_1335626_spike_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__epi_an_nstop
    ON public.epitope_1335626_spike_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__epi_frag_an_start
    ON public.epitope_1335626_spike_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__epi_frag_an_stop
    ON public.epitope_1335626_spike_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__host_tax_id
    ON public.epitope_1335626_spike_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__host_tax_name
    ON public.epitope_1335626_spike_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__iedb_id
    ON public.epitope_1335626_spike_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__is_linear
    ON public.epitope_1335626_spike_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__mhc_allele
    ON public.epitope_1335626_spike_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__mhc_class_lower
    ON public.epitope_1335626_spike_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__product_lower
    ON public.epitope_1335626_spike_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__response_freq_pos
    ON public.epitope_1335626_spike_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__seq_aa_alt
    ON public.epitope_1335626_spike_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__seq_aa_orig
    ON public.epitope_1335626_spike_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__start_aa_orig
    ON public.epitope_1335626_spike_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__taxon_id
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__taxon_name_lower
    ON public.epitope_1335626_spike_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__variant_aa_length
    ON public.epitope_1335626_spike_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__variant_aa_type
    ON public.epitope_1335626_spike_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_n_host_tax_id
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_host_cell_type
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_host_epi_start
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__virus_host_epi_stop
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__virus_host_is_linear
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_host_mhc_allele
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_host_product
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_spike_protein__vir_host_resp_freq
    ON public.epitope_1335626_spike_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS3 protein
-- 1335626 can be replaced with the virus taxon id, while ns3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_ns3_protein (
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

ALTER TABLE public.epitope_1335626_ns3_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns3_protein__cell_type
    ON public.epitope_1335626_ns3_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_an_start
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_an_nstop
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_frag_an_start
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__epi_frag_an_stop
    ON public.epitope_1335626_ns3_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__host_tax_id
    ON public.epitope_1335626_ns3_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__host_tax_name
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__iedb_id
    ON public.epitope_1335626_ns3_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__is_linear
    ON public.epitope_1335626_ns3_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__mhc_allele
    ON public.epitope_1335626_ns3_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__mhc_class_lower
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__product_lower
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__response_freq_pos
    ON public.epitope_1335626_ns3_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__seq_aa_alt
    ON public.epitope_1335626_ns3_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__seq_aa_orig
    ON public.epitope_1335626_ns3_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__start_aa_orig
    ON public.epitope_1335626_ns3_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__taxon_id
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__taxon_name_lower
    ON public.epitope_1335626_ns3_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__variant_aa_length
    ON public.epitope_1335626_ns3_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__variant_aa_type
    ON public.epitope_1335626_ns3_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_n_host_tax_id
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_host_cell_type
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_host_epi_start
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_epi_stop
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__virus_host_is_linear
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_host_mhc_allele
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_host_product
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns3_protein__vir_host_resp_freq
    ON public.epitope_1335626_ns3_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS4A protein
-- 1335626 can be replaced with the virus taxon id, while ns4a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_ns4a_protein (
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

ALTER TABLE public.epitope_1335626_ns4a_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns4a_protein__cell_type
    ON public.epitope_1335626_ns4a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__epi_an_start
    ON public.epitope_1335626_ns4a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__epi_an_nstop
    ON public.epitope_1335626_ns4a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__epi_frag_an_start
    ON public.epitope_1335626_ns4a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__epi_frag_an_stop
    ON public.epitope_1335626_ns4a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__host_tax_id
    ON public.epitope_1335626_ns4a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__host_tax_name
    ON public.epitope_1335626_ns4a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__iedb_id
    ON public.epitope_1335626_ns4a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__is_linear
    ON public.epitope_1335626_ns4a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__mhc_allele
    ON public.epitope_1335626_ns4a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__mhc_class_lower
    ON public.epitope_1335626_ns4a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__product_lower
    ON public.epitope_1335626_ns4a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__response_freq_pos
    ON public.epitope_1335626_ns4a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__seq_aa_alt
    ON public.epitope_1335626_ns4a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__seq_aa_orig
    ON public.epitope_1335626_ns4a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__start_aa_orig
    ON public.epitope_1335626_ns4a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__taxon_id
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__taxon_name_lower
    ON public.epitope_1335626_ns4a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__variant_aa_length
    ON public.epitope_1335626_ns4a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__variant_aa_type
    ON public.epitope_1335626_ns4a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_n_host_tax_id
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_host_cell_type
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_host_epi_start
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__virus_host_epi_stop
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__virus_host_is_linear
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_host_mhc_allele
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_host_product
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4a_protein__vir_host_resp_freq
    ON public.epitope_1335626_ns4a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS4B protein
-- 1335626 can be replaced with the virus taxon id, while ns4b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_ns4b_protein (
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

ALTER TABLE public.epitope_1335626_ns4b_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns4b_protein__cell_type
    ON public.epitope_1335626_ns4b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__epi_an_start
    ON public.epitope_1335626_ns4b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__epi_an_nstop
    ON public.epitope_1335626_ns4b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__epi_frag_an_start
    ON public.epitope_1335626_ns4b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__epi_frag_an_stop
    ON public.epitope_1335626_ns4b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__host_tax_id
    ON public.epitope_1335626_ns4b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__host_tax_name
    ON public.epitope_1335626_ns4b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__iedb_id
    ON public.epitope_1335626_ns4b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__is_linear
    ON public.epitope_1335626_ns4b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__mhc_allele
    ON public.epitope_1335626_ns4b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__mhc_class_lower
    ON public.epitope_1335626_ns4b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__product_lower
    ON public.epitope_1335626_ns4b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__response_freq_pos
    ON public.epitope_1335626_ns4b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__seq_aa_alt
    ON public.epitope_1335626_ns4b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__seq_aa_orig
    ON public.epitope_1335626_ns4b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__start_aa_orig
    ON public.epitope_1335626_ns4b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__taxon_id
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__taxon_name_lower
    ON public.epitope_1335626_ns4b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__variant_aa_length
    ON public.epitope_1335626_ns4b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__variant_aa_type
    ON public.epitope_1335626_ns4b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_n_host_tax_id
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_host_cell_type
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_host_epi_start
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__virus_host_epi_stop
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__virus_host_is_linear
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_host_mhc_allele
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_host_product
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns4b_protein__vir_host_resp_freq
    ON public.epitope_1335626_ns4b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT NS5 protein
-- 1335626 can be replaced with the virus taxon id, while ns5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_ns5_protein (
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

ALTER TABLE public.epitope_1335626_ns5_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_ns5_protein__cell_type
    ON public.epitope_1335626_ns5_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_an_start
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_an_nstop
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_frag_an_start
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__epi_frag_an_stop
    ON public.epitope_1335626_ns5_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__host_tax_id
    ON public.epitope_1335626_ns5_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__host_tax_name
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__iedb_id
    ON public.epitope_1335626_ns5_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__is_linear
    ON public.epitope_1335626_ns5_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__mhc_allele
    ON public.epitope_1335626_ns5_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__mhc_class_lower
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__product_lower
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__response_freq_pos
    ON public.epitope_1335626_ns5_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__seq_aa_alt
    ON public.epitope_1335626_ns5_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__seq_aa_orig
    ON public.epitope_1335626_ns5_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__start_aa_orig
    ON public.epitope_1335626_ns5_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__taxon_id
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__taxon_name_lower
    ON public.epitope_1335626_ns5_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__variant_aa_length
    ON public.epitope_1335626_ns5_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__variant_aa_type
    ON public.epitope_1335626_ns5_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_n_host_tax_id
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_host_cell_type
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_host_epi_start
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_epi_stop
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__virus_host_is_linear
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_host_mhc_allele
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_host_product
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_ns5_protein__vir_host_resp_freq
    ON public.epitope_1335626_ns5_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT envelope protein
-- 1335626 can be replaced with the virus taxon id, while envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_envelope_protein (
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

ALTER TABLE public.epitope_1335626_envelope_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_envelope_protein__cell_type
    ON public.epitope_1335626_envelope_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__epi_an_start
    ON public.epitope_1335626_envelope_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__epi_an_nstop
    ON public.epitope_1335626_envelope_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__epi_frag_an_start
    ON public.epitope_1335626_envelope_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__epi_frag_an_stop
    ON public.epitope_1335626_envelope_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__host_tax_id
    ON public.epitope_1335626_envelope_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__host_tax_name
    ON public.epitope_1335626_envelope_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__iedb_id
    ON public.epitope_1335626_envelope_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__is_linear
    ON public.epitope_1335626_envelope_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__mhc_allele
    ON public.epitope_1335626_envelope_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__mhc_class_lower
    ON public.epitope_1335626_envelope_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__product_lower
    ON public.epitope_1335626_envelope_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__response_freq_pos
    ON public.epitope_1335626_envelope_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__seq_aa_alt
    ON public.epitope_1335626_envelope_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__seq_aa_orig
    ON public.epitope_1335626_envelope_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__start_aa_orig
    ON public.epitope_1335626_envelope_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__taxon_id
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__taxon_name_lower
    ON public.epitope_1335626_envelope_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__variant_aa_length
    ON public.epitope_1335626_envelope_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__variant_aa_type
    ON public.epitope_1335626_envelope_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_n_host_tax_id
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_host_cell_type
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_host_epi_start
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__virus_host_epi_stop
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__virus_host_is_linear
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_host_mhc_allele
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_host_product
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_envelope_protein__vir_host_resp_freq
    ON public.epitope_1335626_envelope_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT membrane protein
-- 1335626 can be replaced with the virus taxon id, while membrane_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_membrane_protein (
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

ALTER TABLE public.epitope_1335626_membrane_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_membrane_protein__cell_type
    ON public.epitope_1335626_membrane_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__epi_an_start
    ON public.epitope_1335626_membrane_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__epi_an_nstop
    ON public.epitope_1335626_membrane_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__epi_frag_an_start
    ON public.epitope_1335626_membrane_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__epi_frag_an_stop
    ON public.epitope_1335626_membrane_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__host_tax_id
    ON public.epitope_1335626_membrane_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__host_tax_name
    ON public.epitope_1335626_membrane_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__iedb_id
    ON public.epitope_1335626_membrane_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__is_linear
    ON public.epitope_1335626_membrane_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__mhc_allele
    ON public.epitope_1335626_membrane_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__mhc_class_lower
    ON public.epitope_1335626_membrane_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__product_lower
    ON public.epitope_1335626_membrane_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__response_freq_pos
    ON public.epitope_1335626_membrane_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__seq_aa_alt
    ON public.epitope_1335626_membrane_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__seq_aa_orig
    ON public.epitope_1335626_membrane_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__start_aa_orig
    ON public.epitope_1335626_membrane_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__taxon_id
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__taxon_name_lower
    ON public.epitope_1335626_membrane_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__variant_aa_length
    ON public.epitope_1335626_membrane_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__variant_aa_type
    ON public.epitope_1335626_membrane_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_n_host_tax_id
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_host_cell_type
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_host_epi_start
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__virus_host_epi_stop
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__virus_host_is_linear
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_host_mhc_allele
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_host_product
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_membrane_protein__vir_host_resp_freq
    ON public.epitope_1335626_membrane_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT nucleocapsid protein
-- 1335626 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_nucleocapsid_protein (
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

ALTER TABLE public.epitope_1335626_nucleocapsid_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_nucleocapsid_protein__cell_type
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__epi_an_start
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__epi_an_nstop
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__epi_frag_an_start
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__epi_frag_an_stop
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__host_tax_id
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__host_tax_name
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__iedb_id
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__is_linear
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__mhc_allele
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__mhc_class_lower
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__product_lower
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__response_freq_pos
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__seq_aa_alt
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__seq_aa_orig
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__start_aa_orig
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__taxon_id
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__taxon_name_lower
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__variant_aa_length
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__variant_aa_type
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_n_host_tax_id
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_host_cell_type
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_host_epi_start
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__virus_host_epi_stop
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__virus_host_is_linear
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_host_mhc_allele
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_host_product
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_nucleocapsid_protein__vir_host_resp_freq
    ON public.epitope_1335626_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR mers and PROT ORF8b protein
-- 1335626 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_1335626_orf8b_protein (
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

ALTER TABLE public.epitope_1335626_orf8b_protein
    OWNER TO geco;

CREATE INDEX epi_1335626_orf8b_protein__cell_type
    ON public.epitope_1335626_orf8b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__epi_an_start
    ON public.epitope_1335626_orf8b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__epi_an_nstop
    ON public.epitope_1335626_orf8b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__epi_frag_an_start
    ON public.epitope_1335626_orf8b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__epi_frag_an_stop
    ON public.epitope_1335626_orf8b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__host_tax_id
    ON public.epitope_1335626_orf8b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__host_tax_name
    ON public.epitope_1335626_orf8b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__iedb_id
    ON public.epitope_1335626_orf8b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__is_linear
    ON public.epitope_1335626_orf8b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__mhc_allele
    ON public.epitope_1335626_orf8b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__mhc_class_lower
    ON public.epitope_1335626_orf8b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__product_lower
    ON public.epitope_1335626_orf8b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__response_freq_pos
    ON public.epitope_1335626_orf8b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__seq_aa_alt
    ON public.epitope_1335626_orf8b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__seq_aa_orig
    ON public.epitope_1335626_orf8b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__start_aa_orig
    ON public.epitope_1335626_orf8b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__taxon_id
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__taxon_name_lower
    ON public.epitope_1335626_orf8b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__variant_aa_length
    ON public.epitope_1335626_orf8b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__variant_aa_type
    ON public.epitope_1335626_orf8b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_n_host_tax_id
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_host_cell_type
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_host_epi_start
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__virus_host_epi_stop
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__virus_host_is_linear
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_host_mhc_allele
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_host_product
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_1335626_orf8b_protein__vir_host_resp_freq
    ON public.epitope_1335626_orf8b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


