-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1ab polyprotein
-- 694009 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf1ab_polyprotein (
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

ALTER TABLE public.epitope_694009_orf1ab_polyprotein
    OWNER TO geco;

CREATE INDEX epi_694009_orf1ab_polyprotein__cell_type
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__epi_an_start
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__epi_an_nstop
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__epi_frag_an_start
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__epi_frag_an_stop
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__host_tax_id
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__host_tax_name
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__iedb_id
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__is_linear
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__mhc_allele
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__mhc_class_lower
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__product_lower
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__response_freq_pos
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__seq_aa_alt
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__seq_aa_orig
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__start_aa_orig
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__taxon_id
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__taxon_name_lower
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__variant_aa_length
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__variant_aa_type
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_n_host_tax_id
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_host_cell_type
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_host_epi_start
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__virus_host_epi_stop
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__virus_host_is_linear
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_host_mhc_allele
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_host_product
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1ab_polyprotein__vir_host_resp_freq
    ON public.epitope_694009_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT RNA-dependent RNA polymerase
-- 694009 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_rna_dependent_rna_polymerase (
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

ALTER TABLE public.epitope_694009_rna_dependent_rna_polymerase
    OWNER TO geco;

CREATE INDEX epi_694009_rna_dependent_rna_polymerase__cell_type
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__epi_an_start
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__epi_an_nstop
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__epi_frag_an_start
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__epi_frag_an_stop
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__host_tax_id
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__host_tax_name
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__iedb_id
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__is_linear
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__mhc_allele
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__mhc_class_lower
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__product_lower
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__response_freq_pos
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__seq_aa_alt
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__seq_aa_orig
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__start_aa_orig
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__taxon_id
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__taxon_name_lower
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__variant_aa_length
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__variant_aa_type
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_n_host_tax_id
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_host_cell_type
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_host_epi_start
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__virus_host_epi_stop
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__virus_host_is_linear
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_host_mhc_allele
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_host_product
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_rna_dependent_rna_polymerase__vir_host_resp_freq
    ON public.epitope_694009_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT helicase/NTPase
-- 694009 can be replaced with the virus taxon id, while helicase_ntpase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_helicase_ntpase (
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

ALTER TABLE public.epitope_694009_helicase_ntpase
    OWNER TO geco;

CREATE INDEX epi_694009_helicase_ntpase__cell_type
    ON public.epitope_694009_helicase_ntpase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__epi_an_start
    ON public.epitope_694009_helicase_ntpase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__epi_an_nstop
    ON public.epitope_694009_helicase_ntpase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__epi_frag_an_start
    ON public.epitope_694009_helicase_ntpase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__epi_frag_an_stop
    ON public.epitope_694009_helicase_ntpase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__host_tax_id
    ON public.epitope_694009_helicase_ntpase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__host_tax_name
    ON public.epitope_694009_helicase_ntpase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__iedb_id
    ON public.epitope_694009_helicase_ntpase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__is_linear
    ON public.epitope_694009_helicase_ntpase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__mhc_allele
    ON public.epitope_694009_helicase_ntpase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__mhc_class_lower
    ON public.epitope_694009_helicase_ntpase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__product_lower
    ON public.epitope_694009_helicase_ntpase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__response_freq_pos
    ON public.epitope_694009_helicase_ntpase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__seq_aa_alt
    ON public.epitope_694009_helicase_ntpase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__seq_aa_orig
    ON public.epitope_694009_helicase_ntpase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__start_aa_orig
    ON public.epitope_694009_helicase_ntpase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__taxon_id
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__taxon_name_lower
    ON public.epitope_694009_helicase_ntpase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__variant_aa_length
    ON public.epitope_694009_helicase_ntpase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__variant_aa_type
    ON public.epitope_694009_helicase_ntpase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_n_host_tax_id
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_host_cell_type
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_host_epi_start
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__virus_host_epi_stop
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__virus_host_is_linear
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_host_mhc_allele
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_host_product
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_helicase_ntpase__vir_host_resp_freq
    ON public.epitope_694009_helicase_ntpase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3' to 5' exonuclease
-- 694009 can be replaced with the virus taxon id, while 3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_3_to_5_exonuclease (
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

ALTER TABLE public.epitope_694009_3_to_5_exonuclease
    OWNER TO geco;

CREATE INDEX epi_694009_3_to_5_exonuclease__cell_type
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__epi_an_start
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__epi_an_nstop
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__epi_frag_an_start
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__epi_frag_an_stop
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__host_tax_id
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__host_tax_name
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__iedb_id
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__is_linear
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__mhc_allele
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__mhc_class_lower
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__product_lower
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__response_freq_pos
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__seq_aa_alt
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__seq_aa_orig
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__start_aa_orig
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__taxon_id
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__taxon_name_lower
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__variant_aa_length
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__variant_aa_type
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_n_host_tax_id
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_host_cell_type
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_host_epi_start
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__virus_host_epi_stop
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__virus_host_is_linear
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_host_mhc_allele
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_host_product
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3_to_5_exonuclease__vir_host_resp_freq
    ON public.epitope_694009_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT endoribonuclease
-- 694009 can be replaced with the virus taxon id, while endoribonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_endoribonuclease (
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

ALTER TABLE public.epitope_694009_endoribonuclease
    OWNER TO geco;

CREATE INDEX epi_694009_endoribonuclease__cell_type
    ON public.epitope_694009_endoribonuclease USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__epi_an_start
    ON public.epitope_694009_endoribonuclease USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__epi_an_nstop
    ON public.epitope_694009_endoribonuclease USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__epi_frag_an_start
    ON public.epitope_694009_endoribonuclease USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__epi_frag_an_stop
    ON public.epitope_694009_endoribonuclease USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__host_tax_id
    ON public.epitope_694009_endoribonuclease USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__host_tax_name
    ON public.epitope_694009_endoribonuclease USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__iedb_id
    ON public.epitope_694009_endoribonuclease USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__is_linear
    ON public.epitope_694009_endoribonuclease USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__mhc_allele
    ON public.epitope_694009_endoribonuclease USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__mhc_class_lower
    ON public.epitope_694009_endoribonuclease USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__product_lower
    ON public.epitope_694009_endoribonuclease USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__response_freq_pos
    ON public.epitope_694009_endoribonuclease USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__seq_aa_alt
    ON public.epitope_694009_endoribonuclease USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__seq_aa_orig
    ON public.epitope_694009_endoribonuclease USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__start_aa_orig
    ON public.epitope_694009_endoribonuclease USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__taxon_id
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__taxon_name_lower
    ON public.epitope_694009_endoribonuclease USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__variant_aa_length
    ON public.epitope_694009_endoribonuclease USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__variant_aa_type
    ON public.epitope_694009_endoribonuclease USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_n_host_tax_id
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_host_cell_type
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_host_epi_start
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__virus_host_epi_stop
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__virus_host_is_linear
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_host_mhc_allele
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_host_product
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_endoribonuclease__vir_host_resp_freq
    ON public.epitope_694009_endoribonuclease USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 2'-O-MTase
-- 694009 can be replaced with the virus taxon id, while 2_o_mtase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_2_o_mtase (
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

ALTER TABLE public.epitope_694009_2_o_mtase
    OWNER TO geco;

CREATE INDEX epi_694009_2_o_mtase__cell_type
    ON public.epitope_694009_2_o_mtase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_an_start
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_an_nstop
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_frag_an_start
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__epi_frag_an_stop
    ON public.epitope_694009_2_o_mtase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__host_tax_id
    ON public.epitope_694009_2_o_mtase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__host_tax_name
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__iedb_id
    ON public.epitope_694009_2_o_mtase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__is_linear
    ON public.epitope_694009_2_o_mtase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__mhc_allele
    ON public.epitope_694009_2_o_mtase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__mhc_class_lower
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__product_lower
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__response_freq_pos
    ON public.epitope_694009_2_o_mtase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__seq_aa_alt
    ON public.epitope_694009_2_o_mtase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__seq_aa_orig
    ON public.epitope_694009_2_o_mtase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__start_aa_orig
    ON public.epitope_694009_2_o_mtase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__taxon_id
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__taxon_name_lower
    ON public.epitope_694009_2_o_mtase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__variant_aa_length
    ON public.epitope_694009_2_o_mtase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__variant_aa_type
    ON public.epitope_694009_2_o_mtase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_n_host_tax_id
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_host_cell_type
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_host_epi_start
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_epi_stop
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__virus_host_is_linear
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_host_mhc_allele
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_host_product
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_2_o_mtase__vir_host_resp_freq
    ON public.epitope_694009_2_o_mtase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1a polyprotein
-- 694009 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf1a_polyprotein (
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

ALTER TABLE public.epitope_694009_orf1a_polyprotein
    OWNER TO geco;

CREATE INDEX epi_694009_orf1a_polyprotein__cell_type
    ON public.epitope_694009_orf1a_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__epi_an_start
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__epi_an_nstop
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__epi_frag_an_start
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__epi_frag_an_stop
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__host_tax_id
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__host_tax_name
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__iedb_id
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__is_linear
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__mhc_allele
    ON public.epitope_694009_orf1a_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__mhc_class_lower
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__product_lower
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__response_freq_pos
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__seq_aa_alt
    ON public.epitope_694009_orf1a_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__seq_aa_orig
    ON public.epitope_694009_orf1a_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__start_aa_orig
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__taxon_id
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__taxon_name_lower
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__variant_aa_length
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__variant_aa_type
    ON public.epitope_694009_orf1a_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_n_host_tax_id
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_host_cell_type
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_host_epi_start
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__virus_host_epi_stop
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__virus_host_is_linear
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_host_mhc_allele
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_host_product
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf1a_polyprotein__vir_host_resp_freq
    ON public.epitope_694009_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp1
-- 694009 can be replaced with the virus taxon id, while nsp1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp1 (
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

ALTER TABLE public.epitope_694009_nsp1
    OWNER TO geco;

CREATE INDEX epi_694009_nsp1__cell_type
    ON public.epitope_694009_nsp1 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_an_start
    ON public.epitope_694009_nsp1 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_an_nstop
    ON public.epitope_694009_nsp1 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_frag_an_start
    ON public.epitope_694009_nsp1 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__epi_frag_an_stop
    ON public.epitope_694009_nsp1 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__host_tax_id
    ON public.epitope_694009_nsp1 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__host_tax_name
    ON public.epitope_694009_nsp1 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__iedb_id
    ON public.epitope_694009_nsp1 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__is_linear
    ON public.epitope_694009_nsp1 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__mhc_allele
    ON public.epitope_694009_nsp1 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__mhc_class_lower
    ON public.epitope_694009_nsp1 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__product_lower
    ON public.epitope_694009_nsp1 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__response_freq_pos
    ON public.epitope_694009_nsp1 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__seq_aa_alt
    ON public.epitope_694009_nsp1 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__seq_aa_orig
    ON public.epitope_694009_nsp1 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__start_aa_orig
    ON public.epitope_694009_nsp1 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__taxon_id
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__taxon_name_lower
    ON public.epitope_694009_nsp1 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__variant_aa_length
    ON public.epitope_694009_nsp1 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__variant_aa_type
    ON public.epitope_694009_nsp1 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_n_host_tax_id
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_host_cell_type
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_host_epi_start
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_epi_stop
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__virus_host_is_linear
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_host_mhc_allele
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_host_product
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp1__vir_host_resp_freq
    ON public.epitope_694009_nsp1 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp2
-- 694009 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp2 (
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

ALTER TABLE public.epitope_694009_nsp2
    OWNER TO geco;

CREATE INDEX epi_694009_nsp2__cell_type
    ON public.epitope_694009_nsp2 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_an_start
    ON public.epitope_694009_nsp2 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_an_nstop
    ON public.epitope_694009_nsp2 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_frag_an_start
    ON public.epitope_694009_nsp2 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__epi_frag_an_stop
    ON public.epitope_694009_nsp2 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__host_tax_id
    ON public.epitope_694009_nsp2 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__host_tax_name
    ON public.epitope_694009_nsp2 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__iedb_id
    ON public.epitope_694009_nsp2 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__is_linear
    ON public.epitope_694009_nsp2 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__mhc_allele
    ON public.epitope_694009_nsp2 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__mhc_class_lower
    ON public.epitope_694009_nsp2 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__product_lower
    ON public.epitope_694009_nsp2 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__response_freq_pos
    ON public.epitope_694009_nsp2 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__seq_aa_alt
    ON public.epitope_694009_nsp2 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__seq_aa_orig
    ON public.epitope_694009_nsp2 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__start_aa_orig
    ON public.epitope_694009_nsp2 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__taxon_id
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__taxon_name_lower
    ON public.epitope_694009_nsp2 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__variant_aa_length
    ON public.epitope_694009_nsp2 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__variant_aa_type
    ON public.epitope_694009_nsp2 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_n_host_tax_id
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_host_cell_type
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_host_epi_start
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_epi_stop
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__virus_host_is_linear
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_host_mhc_allele
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_host_product
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp2__vir_host_resp_freq
    ON public.epitope_694009_nsp2 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp3
-- 694009 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp3 (
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

ALTER TABLE public.epitope_694009_nsp3
    OWNER TO geco;

CREATE INDEX epi_694009_nsp3__cell_type
    ON public.epitope_694009_nsp3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_an_start
    ON public.epitope_694009_nsp3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_an_nstop
    ON public.epitope_694009_nsp3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_frag_an_start
    ON public.epitope_694009_nsp3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__epi_frag_an_stop
    ON public.epitope_694009_nsp3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__host_tax_id
    ON public.epitope_694009_nsp3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__host_tax_name
    ON public.epitope_694009_nsp3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__iedb_id
    ON public.epitope_694009_nsp3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__is_linear
    ON public.epitope_694009_nsp3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__mhc_allele
    ON public.epitope_694009_nsp3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__mhc_class_lower
    ON public.epitope_694009_nsp3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__product_lower
    ON public.epitope_694009_nsp3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__response_freq_pos
    ON public.epitope_694009_nsp3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__seq_aa_alt
    ON public.epitope_694009_nsp3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__seq_aa_orig
    ON public.epitope_694009_nsp3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__start_aa_orig
    ON public.epitope_694009_nsp3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__taxon_id
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__taxon_name_lower
    ON public.epitope_694009_nsp3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__variant_aa_length
    ON public.epitope_694009_nsp3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__variant_aa_type
    ON public.epitope_694009_nsp3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_n_host_tax_id
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_host_cell_type
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_host_epi_start
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_epi_stop
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__virus_host_is_linear
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_host_mhc_allele
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_host_product
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp3__vir_host_resp_freq
    ON public.epitope_694009_nsp3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp4
-- 694009 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp4 (
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

ALTER TABLE public.epitope_694009_nsp4
    OWNER TO geco;

CREATE INDEX epi_694009_nsp4__cell_type
    ON public.epitope_694009_nsp4 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_an_start
    ON public.epitope_694009_nsp4 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_an_nstop
    ON public.epitope_694009_nsp4 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_frag_an_start
    ON public.epitope_694009_nsp4 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__epi_frag_an_stop
    ON public.epitope_694009_nsp4 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__host_tax_id
    ON public.epitope_694009_nsp4 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__host_tax_name
    ON public.epitope_694009_nsp4 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__iedb_id
    ON public.epitope_694009_nsp4 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__is_linear
    ON public.epitope_694009_nsp4 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__mhc_allele
    ON public.epitope_694009_nsp4 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__mhc_class_lower
    ON public.epitope_694009_nsp4 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__product_lower
    ON public.epitope_694009_nsp4 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__response_freq_pos
    ON public.epitope_694009_nsp4 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__seq_aa_alt
    ON public.epitope_694009_nsp4 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__seq_aa_orig
    ON public.epitope_694009_nsp4 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__start_aa_orig
    ON public.epitope_694009_nsp4 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__taxon_id
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__taxon_name_lower
    ON public.epitope_694009_nsp4 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__variant_aa_length
    ON public.epitope_694009_nsp4 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__variant_aa_type
    ON public.epitope_694009_nsp4 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_n_host_tax_id
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_host_cell_type
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_host_epi_start
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_epi_stop
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__virus_host_is_linear
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_host_mhc_allele
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_host_product
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp4__vir_host_resp_freq
    ON public.epitope_694009_nsp4 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3C-like protease
-- 694009 can be replaced with the virus taxon id, while 3c_like_protease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_3c_like_protease (
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

ALTER TABLE public.epitope_694009_3c_like_protease
    OWNER TO geco;

CREATE INDEX epi_694009_3c_like_protease__cell_type
    ON public.epitope_694009_3c_like_protease USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__epi_an_start
    ON public.epitope_694009_3c_like_protease USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__epi_an_nstop
    ON public.epitope_694009_3c_like_protease USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__epi_frag_an_start
    ON public.epitope_694009_3c_like_protease USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__epi_frag_an_stop
    ON public.epitope_694009_3c_like_protease USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__host_tax_id
    ON public.epitope_694009_3c_like_protease USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__host_tax_name
    ON public.epitope_694009_3c_like_protease USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__iedb_id
    ON public.epitope_694009_3c_like_protease USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__is_linear
    ON public.epitope_694009_3c_like_protease USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__mhc_allele
    ON public.epitope_694009_3c_like_protease USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__mhc_class_lower
    ON public.epitope_694009_3c_like_protease USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__product_lower
    ON public.epitope_694009_3c_like_protease USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__response_freq_pos
    ON public.epitope_694009_3c_like_protease USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__seq_aa_alt
    ON public.epitope_694009_3c_like_protease USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__seq_aa_orig
    ON public.epitope_694009_3c_like_protease USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__start_aa_orig
    ON public.epitope_694009_3c_like_protease USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__taxon_id
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__taxon_name_lower
    ON public.epitope_694009_3c_like_protease USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__variant_aa_length
    ON public.epitope_694009_3c_like_protease USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__variant_aa_type
    ON public.epitope_694009_3c_like_protease USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_n_host_tax_id
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_host_cell_type
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_host_epi_start
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__virus_host_epi_stop
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__virus_host_is_linear
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_host_mhc_allele
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_host_product
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_3c_like_protease__vir_host_resp_freq
    ON public.epitope_694009_3c_like_protease USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp6
-- 694009 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp6 (
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

ALTER TABLE public.epitope_694009_nsp6
    OWNER TO geco;

CREATE INDEX epi_694009_nsp6__cell_type
    ON public.epitope_694009_nsp6 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_an_start
    ON public.epitope_694009_nsp6 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_an_nstop
    ON public.epitope_694009_nsp6 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_frag_an_start
    ON public.epitope_694009_nsp6 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__epi_frag_an_stop
    ON public.epitope_694009_nsp6 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__host_tax_id
    ON public.epitope_694009_nsp6 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__host_tax_name
    ON public.epitope_694009_nsp6 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__iedb_id
    ON public.epitope_694009_nsp6 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__is_linear
    ON public.epitope_694009_nsp6 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__mhc_allele
    ON public.epitope_694009_nsp6 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__mhc_class_lower
    ON public.epitope_694009_nsp6 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__product_lower
    ON public.epitope_694009_nsp6 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__response_freq_pos
    ON public.epitope_694009_nsp6 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__seq_aa_alt
    ON public.epitope_694009_nsp6 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__seq_aa_orig
    ON public.epitope_694009_nsp6 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__start_aa_orig
    ON public.epitope_694009_nsp6 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__taxon_id
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__taxon_name_lower
    ON public.epitope_694009_nsp6 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__variant_aa_length
    ON public.epitope_694009_nsp6 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__variant_aa_type
    ON public.epitope_694009_nsp6 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_n_host_tax_id
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_host_cell_type
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_host_epi_start
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_epi_stop
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__virus_host_is_linear
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_host_mhc_allele
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_host_product
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp6__vir_host_resp_freq
    ON public.epitope_694009_nsp6 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp7
-- 694009 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp7 (
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

ALTER TABLE public.epitope_694009_nsp7
    OWNER TO geco;

CREATE INDEX epi_694009_nsp7__cell_type
    ON public.epitope_694009_nsp7 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_an_start
    ON public.epitope_694009_nsp7 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_an_nstop
    ON public.epitope_694009_nsp7 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_frag_an_start
    ON public.epitope_694009_nsp7 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__epi_frag_an_stop
    ON public.epitope_694009_nsp7 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__host_tax_id
    ON public.epitope_694009_nsp7 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__host_tax_name
    ON public.epitope_694009_nsp7 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__iedb_id
    ON public.epitope_694009_nsp7 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__is_linear
    ON public.epitope_694009_nsp7 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__mhc_allele
    ON public.epitope_694009_nsp7 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__mhc_class_lower
    ON public.epitope_694009_nsp7 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__product_lower
    ON public.epitope_694009_nsp7 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__response_freq_pos
    ON public.epitope_694009_nsp7 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__seq_aa_alt
    ON public.epitope_694009_nsp7 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__seq_aa_orig
    ON public.epitope_694009_nsp7 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__start_aa_orig
    ON public.epitope_694009_nsp7 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__taxon_id
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__taxon_name_lower
    ON public.epitope_694009_nsp7 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__variant_aa_length
    ON public.epitope_694009_nsp7 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__variant_aa_type
    ON public.epitope_694009_nsp7 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_n_host_tax_id
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_host_cell_type
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_host_epi_start
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_epi_stop
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__virus_host_is_linear
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_host_mhc_allele
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_host_product
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp7__vir_host_resp_freq
    ON public.epitope_694009_nsp7 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp8
-- 694009 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp8 (
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

ALTER TABLE public.epitope_694009_nsp8
    OWNER TO geco;

CREATE INDEX epi_694009_nsp8__cell_type
    ON public.epitope_694009_nsp8 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_an_start
    ON public.epitope_694009_nsp8 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_an_nstop
    ON public.epitope_694009_nsp8 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_frag_an_start
    ON public.epitope_694009_nsp8 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__epi_frag_an_stop
    ON public.epitope_694009_nsp8 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__host_tax_id
    ON public.epitope_694009_nsp8 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__host_tax_name
    ON public.epitope_694009_nsp8 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__iedb_id
    ON public.epitope_694009_nsp8 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__is_linear
    ON public.epitope_694009_nsp8 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__mhc_allele
    ON public.epitope_694009_nsp8 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__mhc_class_lower
    ON public.epitope_694009_nsp8 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__product_lower
    ON public.epitope_694009_nsp8 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__response_freq_pos
    ON public.epitope_694009_nsp8 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__seq_aa_alt
    ON public.epitope_694009_nsp8 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__seq_aa_orig
    ON public.epitope_694009_nsp8 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__start_aa_orig
    ON public.epitope_694009_nsp8 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__taxon_id
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__taxon_name_lower
    ON public.epitope_694009_nsp8 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__variant_aa_length
    ON public.epitope_694009_nsp8 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__variant_aa_type
    ON public.epitope_694009_nsp8 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_n_host_tax_id
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_host_cell_type
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_host_epi_start
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_epi_stop
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__virus_host_is_linear
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_host_mhc_allele
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_host_product
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp8__vir_host_resp_freq
    ON public.epitope_694009_nsp8 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp9
-- 694009 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp9 (
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

ALTER TABLE public.epitope_694009_nsp9
    OWNER TO geco;

CREATE INDEX epi_694009_nsp9__cell_type
    ON public.epitope_694009_nsp9 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_an_start
    ON public.epitope_694009_nsp9 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_an_nstop
    ON public.epitope_694009_nsp9 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_frag_an_start
    ON public.epitope_694009_nsp9 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__epi_frag_an_stop
    ON public.epitope_694009_nsp9 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__host_tax_id
    ON public.epitope_694009_nsp9 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__host_tax_name
    ON public.epitope_694009_nsp9 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__iedb_id
    ON public.epitope_694009_nsp9 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__is_linear
    ON public.epitope_694009_nsp9 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__mhc_allele
    ON public.epitope_694009_nsp9 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__mhc_class_lower
    ON public.epitope_694009_nsp9 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__product_lower
    ON public.epitope_694009_nsp9 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__response_freq_pos
    ON public.epitope_694009_nsp9 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__seq_aa_alt
    ON public.epitope_694009_nsp9 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__seq_aa_orig
    ON public.epitope_694009_nsp9 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__start_aa_orig
    ON public.epitope_694009_nsp9 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__taxon_id
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__taxon_name_lower
    ON public.epitope_694009_nsp9 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__variant_aa_length
    ON public.epitope_694009_nsp9 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__variant_aa_type
    ON public.epitope_694009_nsp9 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_n_host_tax_id
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_host_cell_type
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_host_epi_start
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_epi_stop
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__virus_host_is_linear
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_host_mhc_allele
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_host_product
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp9__vir_host_resp_freq
    ON public.epitope_694009_nsp9 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp10
-- 694009 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nsp10 (
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

ALTER TABLE public.epitope_694009_nsp10
    OWNER TO geco;

CREATE INDEX epi_694009_nsp10__cell_type
    ON public.epitope_694009_nsp10 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_an_start
    ON public.epitope_694009_nsp10 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_an_nstop
    ON public.epitope_694009_nsp10 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_frag_an_start
    ON public.epitope_694009_nsp10 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__epi_frag_an_stop
    ON public.epitope_694009_nsp10 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__host_tax_id
    ON public.epitope_694009_nsp10 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__host_tax_name
    ON public.epitope_694009_nsp10 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__iedb_id
    ON public.epitope_694009_nsp10 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__is_linear
    ON public.epitope_694009_nsp10 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__mhc_allele
    ON public.epitope_694009_nsp10 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__mhc_class_lower
    ON public.epitope_694009_nsp10 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__product_lower
    ON public.epitope_694009_nsp10 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__response_freq_pos
    ON public.epitope_694009_nsp10 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__seq_aa_alt
    ON public.epitope_694009_nsp10 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__seq_aa_orig
    ON public.epitope_694009_nsp10 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__start_aa_orig
    ON public.epitope_694009_nsp10 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__taxon_id
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__taxon_name_lower
    ON public.epitope_694009_nsp10 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__variant_aa_length
    ON public.epitope_694009_nsp10 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__variant_aa_type
    ON public.epitope_694009_nsp10 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_n_host_tax_id
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_host_cell_type
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_host_epi_start
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_epi_stop
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__virus_host_is_linear
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_host_mhc_allele
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_host_product
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nsp10__vir_host_resp_freq
    ON public.epitope_694009_nsp10 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ndp11
-- 694009 can be replaced with the virus taxon id, while ndp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_ndp11 (
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

ALTER TABLE public.epitope_694009_ndp11
    OWNER TO geco;

CREATE INDEX epi_694009_ndp11__cell_type
    ON public.epitope_694009_ndp11 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_an_start
    ON public.epitope_694009_ndp11 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_an_nstop
    ON public.epitope_694009_ndp11 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_frag_an_start
    ON public.epitope_694009_ndp11 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__epi_frag_an_stop
    ON public.epitope_694009_ndp11 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__host_tax_id
    ON public.epitope_694009_ndp11 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__host_tax_name
    ON public.epitope_694009_ndp11 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__iedb_id
    ON public.epitope_694009_ndp11 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__is_linear
    ON public.epitope_694009_ndp11 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__mhc_allele
    ON public.epitope_694009_ndp11 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__mhc_class_lower
    ON public.epitope_694009_ndp11 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__product_lower
    ON public.epitope_694009_ndp11 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__response_freq_pos
    ON public.epitope_694009_ndp11 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__seq_aa_alt
    ON public.epitope_694009_ndp11 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__seq_aa_orig
    ON public.epitope_694009_ndp11 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__start_aa_orig
    ON public.epitope_694009_ndp11 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__taxon_id
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__taxon_name_lower
    ON public.epitope_694009_ndp11 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__variant_aa_length
    ON public.epitope_694009_ndp11 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__variant_aa_type
    ON public.epitope_694009_ndp11 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_n_host_tax_id
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_host_cell_type
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_host_epi_start
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_epi_stop
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__virus_host_is_linear
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_host_mhc_allele
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_host_product
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_ndp11__vir_host_resp_freq
    ON public.epitope_694009_ndp11 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT spike glycoprotein
-- 694009 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_spike_glycoprotein (
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

ALTER TABLE public.epitope_694009_spike_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_694009_spike_glycoprotein__cell_type
    ON public.epitope_694009_spike_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__epi_an_start
    ON public.epitope_694009_spike_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__epi_an_nstop
    ON public.epitope_694009_spike_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__epi_frag_an_start
    ON public.epitope_694009_spike_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__epi_frag_an_stop
    ON public.epitope_694009_spike_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__host_tax_id
    ON public.epitope_694009_spike_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__host_tax_name
    ON public.epitope_694009_spike_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__iedb_id
    ON public.epitope_694009_spike_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__is_linear
    ON public.epitope_694009_spike_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__mhc_allele
    ON public.epitope_694009_spike_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__mhc_class_lower
    ON public.epitope_694009_spike_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__product_lower
    ON public.epitope_694009_spike_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__response_freq_pos
    ON public.epitope_694009_spike_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__seq_aa_alt
    ON public.epitope_694009_spike_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__seq_aa_orig
    ON public.epitope_694009_spike_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__start_aa_orig
    ON public.epitope_694009_spike_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__taxon_id
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__taxon_name_lower
    ON public.epitope_694009_spike_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__variant_aa_length
    ON public.epitope_694009_spike_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__variant_aa_type
    ON public.epitope_694009_spike_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_n_host_tax_id
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_host_cell_type
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_host_epi_start
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__virus_host_epi_stop
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__virus_host_is_linear
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_host_mhc_allele
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_host_product
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_spike_glycoprotein__vir_host_resp_freq
    ON public.epitope_694009_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3a protein
-- 694009 can be replaced with the virus taxon id, while orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf3a_protein (
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

ALTER TABLE public.epitope_694009_orf3a_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf3a_protein__cell_type
    ON public.epitope_694009_orf3a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__epi_an_start
    ON public.epitope_694009_orf3a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__epi_an_nstop
    ON public.epitope_694009_orf3a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__epi_frag_an_start
    ON public.epitope_694009_orf3a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__epi_frag_an_stop
    ON public.epitope_694009_orf3a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__host_tax_id
    ON public.epitope_694009_orf3a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__host_tax_name
    ON public.epitope_694009_orf3a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__iedb_id
    ON public.epitope_694009_orf3a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__is_linear
    ON public.epitope_694009_orf3a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__mhc_allele
    ON public.epitope_694009_orf3a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__mhc_class_lower
    ON public.epitope_694009_orf3a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__product_lower
    ON public.epitope_694009_orf3a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__response_freq_pos
    ON public.epitope_694009_orf3a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__seq_aa_alt
    ON public.epitope_694009_orf3a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__seq_aa_orig
    ON public.epitope_694009_orf3a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__start_aa_orig
    ON public.epitope_694009_orf3a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__taxon_id
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__taxon_name_lower
    ON public.epitope_694009_orf3a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__variant_aa_length
    ON public.epitope_694009_orf3a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__variant_aa_type
    ON public.epitope_694009_orf3a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_host_cell_type
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_host_epi_start
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__virus_host_epi_stop
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__virus_host_is_linear
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_host_product
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3a_protein__vir_host_resp_freq
    ON public.epitope_694009_orf3a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3b protein
-- 694009 can be replaced with the virus taxon id, while orf3b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf3b_protein (
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

ALTER TABLE public.epitope_694009_orf3b_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf3b_protein__cell_type
    ON public.epitope_694009_orf3b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__epi_an_start
    ON public.epitope_694009_orf3b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__epi_an_nstop
    ON public.epitope_694009_orf3b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__epi_frag_an_start
    ON public.epitope_694009_orf3b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__epi_frag_an_stop
    ON public.epitope_694009_orf3b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__host_tax_id
    ON public.epitope_694009_orf3b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__host_tax_name
    ON public.epitope_694009_orf3b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__iedb_id
    ON public.epitope_694009_orf3b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__is_linear
    ON public.epitope_694009_orf3b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__mhc_allele
    ON public.epitope_694009_orf3b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__mhc_class_lower
    ON public.epitope_694009_orf3b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__product_lower
    ON public.epitope_694009_orf3b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__response_freq_pos
    ON public.epitope_694009_orf3b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__seq_aa_alt
    ON public.epitope_694009_orf3b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__seq_aa_orig
    ON public.epitope_694009_orf3b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__start_aa_orig
    ON public.epitope_694009_orf3b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__taxon_id
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__taxon_name_lower
    ON public.epitope_694009_orf3b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__variant_aa_length
    ON public.epitope_694009_orf3b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__variant_aa_type
    ON public.epitope_694009_orf3b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_host_cell_type
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_host_epi_start
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__virus_host_epi_stop
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__virus_host_is_linear
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_host_product
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf3b_protein__vir_host_resp_freq
    ON public.epitope_694009_orf3b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT small envelope protein
-- 694009 can be replaced with the virus taxon id, while small_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_small_envelope_protein (
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

ALTER TABLE public.epitope_694009_small_envelope_protein
    OWNER TO geco;

CREATE INDEX epi_694009_small_envelope_protein__cell_type
    ON public.epitope_694009_small_envelope_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__epi_an_start
    ON public.epitope_694009_small_envelope_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__epi_an_nstop
    ON public.epitope_694009_small_envelope_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__epi_frag_an_start
    ON public.epitope_694009_small_envelope_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__epi_frag_an_stop
    ON public.epitope_694009_small_envelope_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__host_tax_id
    ON public.epitope_694009_small_envelope_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__host_tax_name
    ON public.epitope_694009_small_envelope_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__iedb_id
    ON public.epitope_694009_small_envelope_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__is_linear
    ON public.epitope_694009_small_envelope_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__mhc_allele
    ON public.epitope_694009_small_envelope_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__mhc_class_lower
    ON public.epitope_694009_small_envelope_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__product_lower
    ON public.epitope_694009_small_envelope_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__response_freq_pos
    ON public.epitope_694009_small_envelope_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__seq_aa_alt
    ON public.epitope_694009_small_envelope_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__seq_aa_orig
    ON public.epitope_694009_small_envelope_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__start_aa_orig
    ON public.epitope_694009_small_envelope_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__taxon_id
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__taxon_name_lower
    ON public.epitope_694009_small_envelope_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__variant_aa_length
    ON public.epitope_694009_small_envelope_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__variant_aa_type
    ON public.epitope_694009_small_envelope_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_n_host_tax_id
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_host_cell_type
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_host_epi_start
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__virus_host_epi_stop
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__virus_host_is_linear
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_host_mhc_allele
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_host_product
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_small_envelope_protein__vir_host_resp_freq
    ON public.epitope_694009_small_envelope_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT membrane glycoprotein M
-- 694009 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_membrane_glycoprotein_m (
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

ALTER TABLE public.epitope_694009_membrane_glycoprotein_m
    OWNER TO geco;

CREATE INDEX epi_694009_membrane_glycoprotein_m__cell_type
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__epi_an_start
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__epi_an_nstop
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__epi_frag_an_start
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__epi_frag_an_stop
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__host_tax_id
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__host_tax_name
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__iedb_id
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__is_linear
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__mhc_allele
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__mhc_class_lower
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__product_lower
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__response_freq_pos
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__seq_aa_alt
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__seq_aa_orig
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__start_aa_orig
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__taxon_id
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__taxon_name_lower
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__variant_aa_length
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__variant_aa_type
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_n_host_tax_id
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_host_cell_type
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_host_epi_start
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__virus_host_epi_stop
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__virus_host_is_linear
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_host_mhc_allele
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_host_product
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_membrane_glycoprotein_m__vir_host_resp_freq
    ON public.epitope_694009_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF6 protein
-- 694009 can be replaced with the virus taxon id, while orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf6_protein (
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

ALTER TABLE public.epitope_694009_orf6_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf6_protein__cell_type
    ON public.epitope_694009_orf6_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__epi_an_start
    ON public.epitope_694009_orf6_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__epi_an_nstop
    ON public.epitope_694009_orf6_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__epi_frag_an_start
    ON public.epitope_694009_orf6_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__epi_frag_an_stop
    ON public.epitope_694009_orf6_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__host_tax_id
    ON public.epitope_694009_orf6_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__host_tax_name
    ON public.epitope_694009_orf6_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__iedb_id
    ON public.epitope_694009_orf6_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__is_linear
    ON public.epitope_694009_orf6_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__mhc_allele
    ON public.epitope_694009_orf6_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__mhc_class_lower
    ON public.epitope_694009_orf6_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__product_lower
    ON public.epitope_694009_orf6_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__response_freq_pos
    ON public.epitope_694009_orf6_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__seq_aa_alt
    ON public.epitope_694009_orf6_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__seq_aa_orig
    ON public.epitope_694009_orf6_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__start_aa_orig
    ON public.epitope_694009_orf6_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__taxon_id
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__taxon_name_lower
    ON public.epitope_694009_orf6_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__variant_aa_length
    ON public.epitope_694009_orf6_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__variant_aa_type
    ON public.epitope_694009_orf6_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_host_cell_type
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_host_epi_start
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__virus_host_epi_stop
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__virus_host_is_linear
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_host_product
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf6_protein__vir_host_resp_freq
    ON public.epitope_694009_orf6_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7a protein
-- 694009 can be replaced with the virus taxon id, while orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf7a_protein (
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

ALTER TABLE public.epitope_694009_orf7a_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf7a_protein__cell_type
    ON public.epitope_694009_orf7a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__epi_an_start
    ON public.epitope_694009_orf7a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__epi_an_nstop
    ON public.epitope_694009_orf7a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__epi_frag_an_start
    ON public.epitope_694009_orf7a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__epi_frag_an_stop
    ON public.epitope_694009_orf7a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__host_tax_id
    ON public.epitope_694009_orf7a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__host_tax_name
    ON public.epitope_694009_orf7a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__iedb_id
    ON public.epitope_694009_orf7a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__is_linear
    ON public.epitope_694009_orf7a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__mhc_allele
    ON public.epitope_694009_orf7a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__mhc_class_lower
    ON public.epitope_694009_orf7a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__product_lower
    ON public.epitope_694009_orf7a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__response_freq_pos
    ON public.epitope_694009_orf7a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__seq_aa_alt
    ON public.epitope_694009_orf7a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__seq_aa_orig
    ON public.epitope_694009_orf7a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__start_aa_orig
    ON public.epitope_694009_orf7a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__taxon_id
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__taxon_name_lower
    ON public.epitope_694009_orf7a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__variant_aa_length
    ON public.epitope_694009_orf7a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__variant_aa_type
    ON public.epitope_694009_orf7a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_host_cell_type
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_host_epi_start
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__virus_host_epi_stop
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__virus_host_is_linear
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_host_product
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7a_protein__vir_host_resp_freq
    ON public.epitope_694009_orf7a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7b protein
-- 694009 can be replaced with the virus taxon id, while orf7b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf7b_protein (
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

ALTER TABLE public.epitope_694009_orf7b_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf7b_protein__cell_type
    ON public.epitope_694009_orf7b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__epi_an_start
    ON public.epitope_694009_orf7b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__epi_an_nstop
    ON public.epitope_694009_orf7b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__epi_frag_an_start
    ON public.epitope_694009_orf7b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__epi_frag_an_stop
    ON public.epitope_694009_orf7b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__host_tax_id
    ON public.epitope_694009_orf7b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__host_tax_name
    ON public.epitope_694009_orf7b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__iedb_id
    ON public.epitope_694009_orf7b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__is_linear
    ON public.epitope_694009_orf7b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__mhc_allele
    ON public.epitope_694009_orf7b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__mhc_class_lower
    ON public.epitope_694009_orf7b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__product_lower
    ON public.epitope_694009_orf7b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__response_freq_pos
    ON public.epitope_694009_orf7b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__seq_aa_alt
    ON public.epitope_694009_orf7b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__seq_aa_orig
    ON public.epitope_694009_orf7b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__start_aa_orig
    ON public.epitope_694009_orf7b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__taxon_id
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__taxon_name_lower
    ON public.epitope_694009_orf7b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__variant_aa_length
    ON public.epitope_694009_orf7b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__variant_aa_type
    ON public.epitope_694009_orf7b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_host_cell_type
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_host_epi_start
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__virus_host_epi_stop
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__virus_host_is_linear
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_host_product
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf7b_protein__vir_host_resp_freq
    ON public.epitope_694009_orf7b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8a protein
-- 694009 can be replaced with the virus taxon id, while orf8a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf8a_protein (
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

ALTER TABLE public.epitope_694009_orf8a_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf8a_protein__cell_type
    ON public.epitope_694009_orf8a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__epi_an_start
    ON public.epitope_694009_orf8a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__epi_an_nstop
    ON public.epitope_694009_orf8a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__epi_frag_an_start
    ON public.epitope_694009_orf8a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__epi_frag_an_stop
    ON public.epitope_694009_orf8a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__host_tax_id
    ON public.epitope_694009_orf8a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__host_tax_name
    ON public.epitope_694009_orf8a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__iedb_id
    ON public.epitope_694009_orf8a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__is_linear
    ON public.epitope_694009_orf8a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__mhc_allele
    ON public.epitope_694009_orf8a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__mhc_class_lower
    ON public.epitope_694009_orf8a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__product_lower
    ON public.epitope_694009_orf8a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__response_freq_pos
    ON public.epitope_694009_orf8a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__seq_aa_alt
    ON public.epitope_694009_orf8a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__seq_aa_orig
    ON public.epitope_694009_orf8a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__start_aa_orig
    ON public.epitope_694009_orf8a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__taxon_id
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__taxon_name_lower
    ON public.epitope_694009_orf8a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__variant_aa_length
    ON public.epitope_694009_orf8a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__variant_aa_type
    ON public.epitope_694009_orf8a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_host_cell_type
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_host_epi_start
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__virus_host_epi_stop
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__virus_host_is_linear
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_host_product
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8a_protein__vir_host_resp_freq
    ON public.epitope_694009_orf8a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8b protein
-- 694009 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf8b_protein (
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

ALTER TABLE public.epitope_694009_orf8b_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf8b_protein__cell_type
    ON public.epitope_694009_orf8b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__epi_an_start
    ON public.epitope_694009_orf8b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__epi_an_nstop
    ON public.epitope_694009_orf8b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__epi_frag_an_start
    ON public.epitope_694009_orf8b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__epi_frag_an_stop
    ON public.epitope_694009_orf8b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__host_tax_id
    ON public.epitope_694009_orf8b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__host_tax_name
    ON public.epitope_694009_orf8b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__iedb_id
    ON public.epitope_694009_orf8b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__is_linear
    ON public.epitope_694009_orf8b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__mhc_allele
    ON public.epitope_694009_orf8b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__mhc_class_lower
    ON public.epitope_694009_orf8b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__product_lower
    ON public.epitope_694009_orf8b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__response_freq_pos
    ON public.epitope_694009_orf8b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__seq_aa_alt
    ON public.epitope_694009_orf8b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__seq_aa_orig
    ON public.epitope_694009_orf8b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__start_aa_orig
    ON public.epitope_694009_orf8b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__taxon_id
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__taxon_name_lower
    ON public.epitope_694009_orf8b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__variant_aa_length
    ON public.epitope_694009_orf8b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__variant_aa_type
    ON public.epitope_694009_orf8b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_host_cell_type
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_host_epi_start
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__virus_host_epi_stop
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__virus_host_is_linear
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_host_product
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf8b_protein__vir_host_resp_freq
    ON public.epitope_694009_orf8b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nucleocapsid protein
-- 694009 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_nucleocapsid_protein (
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

ALTER TABLE public.epitope_694009_nucleocapsid_protein
    OWNER TO geco;

CREATE INDEX epi_694009_nucleocapsid_protein__cell_type
    ON public.epitope_694009_nucleocapsid_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__epi_an_start
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__epi_an_nstop
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__epi_frag_an_start
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__epi_frag_an_stop
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__host_tax_id
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__host_tax_name
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__iedb_id
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__is_linear
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__mhc_allele
    ON public.epitope_694009_nucleocapsid_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__mhc_class_lower
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__product_lower
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__response_freq_pos
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__seq_aa_alt
    ON public.epitope_694009_nucleocapsid_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__seq_aa_orig
    ON public.epitope_694009_nucleocapsid_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__start_aa_orig
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__taxon_id
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__taxon_name_lower
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__variant_aa_length
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__variant_aa_type
    ON public.epitope_694009_nucleocapsid_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_n_host_tax_id
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_host_cell_type
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_host_epi_start
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__virus_host_epi_stop
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__virus_host_is_linear
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_host_mhc_allele
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_host_product
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_nucleocapsid_protein__vir_host_resp_freq
    ON public.epitope_694009_nucleocapsid_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9b protein
-- 694009 can be replaced with the virus taxon id, while orf9b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf9b_protein (
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

ALTER TABLE public.epitope_694009_orf9b_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf9b_protein__cell_type
    ON public.epitope_694009_orf9b_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__epi_an_start
    ON public.epitope_694009_orf9b_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__epi_an_nstop
    ON public.epitope_694009_orf9b_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__epi_frag_an_start
    ON public.epitope_694009_orf9b_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__epi_frag_an_stop
    ON public.epitope_694009_orf9b_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__host_tax_id
    ON public.epitope_694009_orf9b_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__host_tax_name
    ON public.epitope_694009_orf9b_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__iedb_id
    ON public.epitope_694009_orf9b_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__is_linear
    ON public.epitope_694009_orf9b_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__mhc_allele
    ON public.epitope_694009_orf9b_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__mhc_class_lower
    ON public.epitope_694009_orf9b_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__product_lower
    ON public.epitope_694009_orf9b_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__response_freq_pos
    ON public.epitope_694009_orf9b_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__seq_aa_alt
    ON public.epitope_694009_orf9b_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__seq_aa_orig
    ON public.epitope_694009_orf9b_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__start_aa_orig
    ON public.epitope_694009_orf9b_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__taxon_id
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__taxon_name_lower
    ON public.epitope_694009_orf9b_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__variant_aa_length
    ON public.epitope_694009_orf9b_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__variant_aa_type
    ON public.epitope_694009_orf9b_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_host_cell_type
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_host_epi_start
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__virus_host_epi_stop
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__virus_host_is_linear
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_host_product
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9b_protein__vir_host_resp_freq
    ON public.epitope_694009_orf9b_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9a protein
-- 694009 can be replaced with the virus taxon id, while orf9a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_694009_orf9a_protein (
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

ALTER TABLE public.epitope_694009_orf9a_protein
    OWNER TO geco;

CREATE INDEX epi_694009_orf9a_protein__cell_type
    ON public.epitope_694009_orf9a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__epi_an_start
    ON public.epitope_694009_orf9a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__epi_an_nstop
    ON public.epitope_694009_orf9a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__epi_frag_an_start
    ON public.epitope_694009_orf9a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__epi_frag_an_stop
    ON public.epitope_694009_orf9a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__host_tax_id
    ON public.epitope_694009_orf9a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__host_tax_name
    ON public.epitope_694009_orf9a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__iedb_id
    ON public.epitope_694009_orf9a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__is_linear
    ON public.epitope_694009_orf9a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__mhc_allele
    ON public.epitope_694009_orf9a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__mhc_class_lower
    ON public.epitope_694009_orf9a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__product_lower
    ON public.epitope_694009_orf9a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__response_freq_pos
    ON public.epitope_694009_orf9a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__seq_aa_alt
    ON public.epitope_694009_orf9a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__seq_aa_orig
    ON public.epitope_694009_orf9a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__start_aa_orig
    ON public.epitope_694009_orf9a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__taxon_id
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__taxon_name_lower
    ON public.epitope_694009_orf9a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__variant_aa_length
    ON public.epitope_694009_orf9a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__variant_aa_type
    ON public.epitope_694009_orf9a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_n_host_tax_id
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_host_cell_type
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_host_epi_start
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__virus_host_epi_stop
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__virus_host_is_linear
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_host_mhc_allele
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_host_product
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_694009_orf9a_protein__vir_host_resp_freq
    ON public.epitope_694009_orf9a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


