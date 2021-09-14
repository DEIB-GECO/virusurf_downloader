-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT polyprotein
-- 11070 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_polyprotein (
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

ALTER TABLE public.epitope_11070_polyprotein
    OWNER TO geco;

CREATE INDEX epi_11070_polyprotein__cell_type
    ON public.epitope_11070_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_an_start
    ON public.epitope_11070_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_an_nstop
    ON public.epitope_11070_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_frag_an_start
    ON public.epitope_11070_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__epi_frag_an_stop
    ON public.epitope_11070_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__host_tax_id
    ON public.epitope_11070_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__host_tax_name
    ON public.epitope_11070_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__iedb_id
    ON public.epitope_11070_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__is_linear
    ON public.epitope_11070_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__mhc_allele
    ON public.epitope_11070_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__mhc_class_lower
    ON public.epitope_11070_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__product_lower
    ON public.epitope_11070_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__response_freq_pos
    ON public.epitope_11070_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__seq_aa_alt
    ON public.epitope_11070_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__seq_aa_orig
    ON public.epitope_11070_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__start_aa_orig
    ON public.epitope_11070_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__taxon_id
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__taxon_name_lower
    ON public.epitope_11070_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__variant_aa_length
    ON public.epitope_11070_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__variant_aa_type
    ON public.epitope_11070_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_n_host_tax_id
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_host_cell_type
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_host_epi_start
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_epi_stop
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__virus_host_is_linear
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_host_mhc_allele
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_host_product
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_polyprotein__vir_host_resp_freq
    ON public.epitope_11070_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT anchored capsid protein ancC
-- 11070 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_anchored_capsid_protein_ancc (
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

ALTER TABLE public.epitope_11070_anchored_capsid_protein_ancc
    OWNER TO geco;

CREATE INDEX epi_11070_anchored_capsid_protein_ancc__cell_type
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__epi_an_start
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__epi_an_nstop
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__epi_frag_an_start
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__epi_frag_an_stop
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__host_tax_id
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__host_tax_name
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__iedb_id
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__is_linear
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__mhc_allele
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__mhc_class_lower
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__product_lower
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__response_freq_pos
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__seq_aa_alt
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__seq_aa_orig
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__start_aa_orig
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__taxon_id
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__taxon_name_lower
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__variant_aa_length
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__variant_aa_type
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_n_host_tax_id
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_host_cell_type
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_host_epi_start
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__virus_host_epi_stop
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__virus_host_is_linear
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_host_mhc_allele
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_host_product
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_anchored_capsid_protein_ancc__vir_host_resp_freq
    ON public.epitope_11070_anchored_capsid_protein_ancc USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT capsid protein C
-- 11070 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_capsid_protein_c (
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

ALTER TABLE public.epitope_11070_capsid_protein_c
    OWNER TO geco;

CREATE INDEX epi_11070_capsid_protein_c__cell_type
    ON public.epitope_11070_capsid_protein_c USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__epi_an_start
    ON public.epitope_11070_capsid_protein_c USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__epi_an_nstop
    ON public.epitope_11070_capsid_protein_c USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__epi_frag_an_start
    ON public.epitope_11070_capsid_protein_c USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__epi_frag_an_stop
    ON public.epitope_11070_capsid_protein_c USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__host_tax_id
    ON public.epitope_11070_capsid_protein_c USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__host_tax_name
    ON public.epitope_11070_capsid_protein_c USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__iedb_id
    ON public.epitope_11070_capsid_protein_c USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__is_linear
    ON public.epitope_11070_capsid_protein_c USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__mhc_allele
    ON public.epitope_11070_capsid_protein_c USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__mhc_class_lower
    ON public.epitope_11070_capsid_protein_c USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__product_lower
    ON public.epitope_11070_capsid_protein_c USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__response_freq_pos
    ON public.epitope_11070_capsid_protein_c USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__seq_aa_alt
    ON public.epitope_11070_capsid_protein_c USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__seq_aa_orig
    ON public.epitope_11070_capsid_protein_c USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__start_aa_orig
    ON public.epitope_11070_capsid_protein_c USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__taxon_id
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__taxon_name_lower
    ON public.epitope_11070_capsid_protein_c USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__variant_aa_length
    ON public.epitope_11070_capsid_protein_c USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__variant_aa_type
    ON public.epitope_11070_capsid_protein_c USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_n_host_tax_id
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_host_cell_type
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_host_epi_start
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__virus_host_epi_stop
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__virus_host_is_linear
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_host_mhc_allele
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_host_product
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_capsid_protein_c__vir_host_resp_freq
    ON public.epitope_11070_capsid_protein_c USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein precursor prM
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_membrane_glycoprotein_precur (
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

ALTER TABLE public.epitope_11070_membrane_glycoprotein_precur
    OWNER TO geco;

CREATE INDEX epi_11070_membrane_glycoprotein_precur__cell_type
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__epi_an_start
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__epi_an_nstop
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__epi_frag_an_start
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__epi_frag_an_stop
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__host_tax_id
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__host_tax_name
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__iedb_id
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__is_linear
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__mhc_allele
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__mhc_class_lower
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__product_lower
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__response_freq_pos
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__seq_aa_alt
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__seq_aa_orig
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__start_aa_orig
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__taxon_id
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__taxon_name_lower
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__variant_aa_length
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__variant_aa_type
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_n_host_tax_id
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_host_cell_type
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_host_epi_start
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__virus_host_epi_stop
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__virus_host_is_linear
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_host_mhc_allele
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_host_product
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_precur__vir_host_resp_freq
    ON public.epitope_11070_membrane_glycoprotein_precur USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT protein pr
-- 11070 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_protein_pr (
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

ALTER TABLE public.epitope_11070_protein_pr
    OWNER TO geco;

CREATE INDEX epi_11070_protein_pr__cell_type
    ON public.epitope_11070_protein_pr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_an_start
    ON public.epitope_11070_protein_pr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_an_nstop
    ON public.epitope_11070_protein_pr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_frag_an_start
    ON public.epitope_11070_protein_pr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__epi_frag_an_stop
    ON public.epitope_11070_protein_pr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__host_tax_id
    ON public.epitope_11070_protein_pr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__host_tax_name
    ON public.epitope_11070_protein_pr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__iedb_id
    ON public.epitope_11070_protein_pr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__is_linear
    ON public.epitope_11070_protein_pr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__mhc_allele
    ON public.epitope_11070_protein_pr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__mhc_class_lower
    ON public.epitope_11070_protein_pr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__product_lower
    ON public.epitope_11070_protein_pr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__response_freq_pos
    ON public.epitope_11070_protein_pr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__seq_aa_alt
    ON public.epitope_11070_protein_pr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__seq_aa_orig
    ON public.epitope_11070_protein_pr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__start_aa_orig
    ON public.epitope_11070_protein_pr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__taxon_id
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__taxon_name_lower
    ON public.epitope_11070_protein_pr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__variant_aa_length
    ON public.epitope_11070_protein_pr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__variant_aa_type
    ON public.epitope_11070_protein_pr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_n_host_tax_id
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_host_cell_type
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_host_epi_start
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_epi_stop
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__virus_host_is_linear
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_host_mhc_allele
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_host_product
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_pr__vir_host_resp_freq
    ON public.epitope_11070_protein_pr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein M
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_membrane_glycoprotein_m (
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

ALTER TABLE public.epitope_11070_membrane_glycoprotein_m
    OWNER TO geco;

CREATE INDEX epi_11070_membrane_glycoprotein_m__cell_type
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__epi_an_start
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__epi_an_nstop
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__epi_frag_an_start
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__epi_frag_an_stop
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__host_tax_id
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__host_tax_name
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__iedb_id
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__is_linear
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__mhc_allele
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__mhc_class_lower
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__product_lower
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__response_freq_pos
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__seq_aa_alt
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__seq_aa_orig
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__start_aa_orig
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__taxon_id
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__taxon_name_lower
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__variant_aa_length
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__variant_aa_type
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_n_host_tax_id
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_host_cell_type
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_host_epi_start
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__virus_host_epi_stop
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__virus_host_is_linear
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_host_mhc_allele
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_host_product
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_membrane_glycoprotein_m__vir_host_resp_freq
    ON public.epitope_11070_membrane_glycoprotein_m USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT envelope protein E
-- 11070 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_envelope_protein_e (
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

ALTER TABLE public.epitope_11070_envelope_protein_e
    OWNER TO geco;

CREATE INDEX epi_11070_envelope_protein_e__cell_type
    ON public.epitope_11070_envelope_protein_e USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__epi_an_start
    ON public.epitope_11070_envelope_protein_e USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__epi_an_nstop
    ON public.epitope_11070_envelope_protein_e USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__epi_frag_an_start
    ON public.epitope_11070_envelope_protein_e USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__epi_frag_an_stop
    ON public.epitope_11070_envelope_protein_e USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__host_tax_id
    ON public.epitope_11070_envelope_protein_e USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__host_tax_name
    ON public.epitope_11070_envelope_protein_e USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__iedb_id
    ON public.epitope_11070_envelope_protein_e USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__is_linear
    ON public.epitope_11070_envelope_protein_e USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__mhc_allele
    ON public.epitope_11070_envelope_protein_e USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__mhc_class_lower
    ON public.epitope_11070_envelope_protein_e USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__product_lower
    ON public.epitope_11070_envelope_protein_e USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__response_freq_pos
    ON public.epitope_11070_envelope_protein_e USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__seq_aa_alt
    ON public.epitope_11070_envelope_protein_e USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__seq_aa_orig
    ON public.epitope_11070_envelope_protein_e USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__start_aa_orig
    ON public.epitope_11070_envelope_protein_e USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__taxon_id
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__taxon_name_lower
    ON public.epitope_11070_envelope_protein_e USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__variant_aa_length
    ON public.epitope_11070_envelope_protein_e USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__variant_aa_type
    ON public.epitope_11070_envelope_protein_e USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_n_host_tax_id
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_host_cell_type
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_host_epi_start
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__virus_host_epi_stop
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__virus_host_is_linear
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_host_mhc_allele
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_host_product
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_envelope_protein_e__vir_host_resp_freq
    ON public.epitope_11070_envelope_protein_e USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS1
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns1 (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns1
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns1__cell_type
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__is_linear
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__product_lower
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns1__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns1 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2A
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns2a (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns2a
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns2a__cell_type
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__is_linear
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__product_lower
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2a__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns2a USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2B
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns2b (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns2b
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns2b__cell_type
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__is_linear
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__product_lower
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns2b__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns2b USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS3
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns3 (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns3
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns3__cell_type
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__is_linear
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__product_lower
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns3__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4A
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns4a (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns4a
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns4a__cell_type
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__is_linear
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__product_lower
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4a__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns4a USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT protein 2K
-- 11070 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_protein_2k (
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

ALTER TABLE public.epitope_11070_protein_2k
    OWNER TO geco;

CREATE INDEX epi_11070_protein_2k__cell_type
    ON public.epitope_11070_protein_2k USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_an_start
    ON public.epitope_11070_protein_2k USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_an_nstop
    ON public.epitope_11070_protein_2k USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_frag_an_start
    ON public.epitope_11070_protein_2k USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__epi_frag_an_stop
    ON public.epitope_11070_protein_2k USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__host_tax_id
    ON public.epitope_11070_protein_2k USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__host_tax_name
    ON public.epitope_11070_protein_2k USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__iedb_id
    ON public.epitope_11070_protein_2k USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__is_linear
    ON public.epitope_11070_protein_2k USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__mhc_allele
    ON public.epitope_11070_protein_2k USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__mhc_class_lower
    ON public.epitope_11070_protein_2k USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__product_lower
    ON public.epitope_11070_protein_2k USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__response_freq_pos
    ON public.epitope_11070_protein_2k USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__seq_aa_alt
    ON public.epitope_11070_protein_2k USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__seq_aa_orig
    ON public.epitope_11070_protein_2k USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__start_aa_orig
    ON public.epitope_11070_protein_2k USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__taxon_id
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__taxon_name_lower
    ON public.epitope_11070_protein_2k USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__variant_aa_length
    ON public.epitope_11070_protein_2k USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__variant_aa_type
    ON public.epitope_11070_protein_2k USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_n_host_tax_id
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_host_cell_type
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_host_epi_start
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_epi_stop
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__virus_host_is_linear
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_host_mhc_allele
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_host_product
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_protein_2k__vir_host_resp_freq
    ON public.epitope_11070_protein_2k USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4B
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_nonstructural_protein_ns4b (
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

ALTER TABLE public.epitope_11070_nonstructural_protein_ns4b
    OWNER TO geco;

CREATE INDEX epi_11070_nonstructural_protein_ns4b__cell_type
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__epi_an_start
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__epi_an_nstop
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__epi_frag_an_start
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__epi_frag_an_stop
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__host_tax_name
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__iedb_id
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__is_linear
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__mhc_class_lower
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__product_lower
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__response_freq_pos
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__seq_aa_alt
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__seq_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__start_aa_orig
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__taxon_id
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__taxon_name_lower
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__variant_aa_length
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__variant_aa_type
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_n_host_tax_id
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_host_cell_type
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_host_epi_start
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__virus_host_epi_stop
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__virus_host_is_linear
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_host_mhc_allele
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_host_product
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_nonstructural_protein_ns4b__vir_host_resp_freq
    ON public.epitope_11070_nonstructural_protein_ns4b USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR dengue_4 and PROT RNA-dependent RNA polymerase NS5
-- 11070 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_11070_rna_dependent_rna_polymerase (
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

ALTER TABLE public.epitope_11070_rna_dependent_rna_polymerase
    OWNER TO geco;

CREATE INDEX epi_11070_rna_dependent_rna_polymerase__cell_type
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__epi_an_start
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__epi_an_nstop
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__epi_frag_an_start
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__epi_frag_an_stop
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__host_tax_id
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__host_tax_name
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__iedb_id
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__is_linear
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__mhc_allele
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__mhc_class_lower
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__product_lower
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__response_freq_pos
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__seq_aa_alt
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__seq_aa_orig
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__start_aa_orig
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__taxon_id
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__taxon_name_lower
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__variant_aa_length
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__variant_aa_type
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_n_host_tax_id
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_host_cell_type
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_host_epi_start
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__virus_host_epi_stop
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__virus_host_is_linear
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_host_mhc_allele
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_host_product
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_11070_rna_dependent_rna_polymerase__vir_host_resp_freq
    ON public.epitope_11070_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


