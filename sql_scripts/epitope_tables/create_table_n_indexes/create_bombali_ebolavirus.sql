-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT nucleoprotein
-- 2010960 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_nucleoprotein (
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

ALTER TABLE public.epitope_2010960_nucleoprotein
    OWNER TO geco;

CREATE INDEX epi_2010960_nucleoprotein__cell_type
    ON public.epitope_2010960_nucleoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__epi_an_start
    ON public.epitope_2010960_nucleoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__epi_an_nstop
    ON public.epitope_2010960_nucleoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__epi_frag_an_start
    ON public.epitope_2010960_nucleoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__epi_frag_an_stop
    ON public.epitope_2010960_nucleoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__host_tax_id
    ON public.epitope_2010960_nucleoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__host_tax_name
    ON public.epitope_2010960_nucleoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__iedb_id
    ON public.epitope_2010960_nucleoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__is_linear
    ON public.epitope_2010960_nucleoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__mhc_allele
    ON public.epitope_2010960_nucleoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__mhc_class_lower
    ON public.epitope_2010960_nucleoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__product_lower
    ON public.epitope_2010960_nucleoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__response_freq_pos
    ON public.epitope_2010960_nucleoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__seq_aa_alt
    ON public.epitope_2010960_nucleoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__seq_aa_orig
    ON public.epitope_2010960_nucleoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__start_aa_orig
    ON public.epitope_2010960_nucleoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__taxon_id
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__taxon_name_lower
    ON public.epitope_2010960_nucleoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__variant_aa_length
    ON public.epitope_2010960_nucleoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__variant_aa_type
    ON public.epitope_2010960_nucleoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_n_host_tax_id
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_host_cell_type
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_host_epi_start
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__virus_host_epi_stop
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__virus_host_is_linear
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_host_mhc_allele
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_host_product
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_nucleoprotein__vir_host_resp_freq
    ON public.epitope_2010960_nucleoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT polymerase complex protein
-- 2010960 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_polymerase_complex_protein (
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

ALTER TABLE public.epitope_2010960_polymerase_complex_protein
    OWNER TO geco;

CREATE INDEX epi_2010960_polymerase_complex_protein__cell_type
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__epi_an_start
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__epi_an_nstop
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__epi_frag_an_start
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__epi_frag_an_stop
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__host_tax_id
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__host_tax_name
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__iedb_id
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__is_linear
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__mhc_allele
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__mhc_class_lower
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__product_lower
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__response_freq_pos
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__seq_aa_alt
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__seq_aa_orig
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__start_aa_orig
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__taxon_id
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__taxon_name_lower
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__variant_aa_length
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__variant_aa_type
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_n_host_tax_id
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_host_cell_type
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_host_epi_start
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__virus_host_epi_stop
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__virus_host_is_linear
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_host_mhc_allele
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_host_product
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_polymerase_complex_protein__vir_host_resp_freq
    ON public.epitope_2010960_polymerase_complex_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT matrix protein
-- 2010960 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_matrix_protein (
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

ALTER TABLE public.epitope_2010960_matrix_protein
    OWNER TO geco;

CREATE INDEX epi_2010960_matrix_protein__cell_type
    ON public.epitope_2010960_matrix_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__epi_an_start
    ON public.epitope_2010960_matrix_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__epi_an_nstop
    ON public.epitope_2010960_matrix_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__epi_frag_an_start
    ON public.epitope_2010960_matrix_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__epi_frag_an_stop
    ON public.epitope_2010960_matrix_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__host_tax_id
    ON public.epitope_2010960_matrix_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__host_tax_name
    ON public.epitope_2010960_matrix_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__iedb_id
    ON public.epitope_2010960_matrix_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__is_linear
    ON public.epitope_2010960_matrix_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__mhc_allele
    ON public.epitope_2010960_matrix_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__mhc_class_lower
    ON public.epitope_2010960_matrix_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__product_lower
    ON public.epitope_2010960_matrix_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__response_freq_pos
    ON public.epitope_2010960_matrix_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__seq_aa_alt
    ON public.epitope_2010960_matrix_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__seq_aa_orig
    ON public.epitope_2010960_matrix_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__start_aa_orig
    ON public.epitope_2010960_matrix_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__taxon_id
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__taxon_name_lower
    ON public.epitope_2010960_matrix_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__variant_aa_length
    ON public.epitope_2010960_matrix_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__variant_aa_type
    ON public.epitope_2010960_matrix_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_n_host_tax_id
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_host_cell_type
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_host_epi_start
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__virus_host_epi_stop
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__virus_host_is_linear
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_host_mhc_allele
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_host_product
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_matrix_protein__vir_host_resp_freq
    ON public.epitope_2010960_matrix_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT spike glycoprotein
-- 2010960 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_spike_glycoprotein (
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

ALTER TABLE public.epitope_2010960_spike_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_2010960_spike_glycoprotein__cell_type
    ON public.epitope_2010960_spike_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__epi_an_start
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__epi_an_nstop
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__epi_frag_an_start
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__epi_frag_an_stop
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__host_tax_id
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__host_tax_name
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__iedb_id
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__is_linear
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__mhc_allele
    ON public.epitope_2010960_spike_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__mhc_class_lower
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__product_lower
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__response_freq_pos
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__seq_aa_alt
    ON public.epitope_2010960_spike_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__seq_aa_orig
    ON public.epitope_2010960_spike_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__start_aa_orig
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__taxon_id
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__taxon_name_lower
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__variant_aa_length
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__variant_aa_type
    ON public.epitope_2010960_spike_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_n_host_tax_id
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_host_cell_type
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_host_epi_start
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__virus_host_epi_stop
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__virus_host_is_linear
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_host_mhc_allele
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_host_product
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_spike_glycoprotein__vir_host_resp_freq
    ON public.epitope_2010960_spike_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT small secreted glycoprotein
-- 2010960 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_small_secreted_glycoprotein (
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

ALTER TABLE public.epitope_2010960_small_secreted_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_2010960_small_secreted_glycoprotein__cell_type
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__epi_an_start
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__epi_an_nstop
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__epi_frag_an_start
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__epi_frag_an_stop
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__host_tax_id
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__host_tax_name
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__iedb_id
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__is_linear
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__mhc_allele
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__mhc_class_lower
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__product_lower
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__response_freq_pos
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__seq_aa_alt
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__seq_aa_orig
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__start_aa_orig
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__taxon_id
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__taxon_name_lower
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__variant_aa_length
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__variant_aa_type
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_n_host_tax_id
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_host_cell_type
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_host_epi_start
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__virus_host_epi_stop
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__virus_host_is_linear
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_host_mhc_allele
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_host_product
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_small_secreted_glycoprotein__vir_host_resp_freq
    ON public.epitope_2010960_small_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT second secreted glycoprotein
-- 2010960 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_second_secreted_glycoprotein (
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

ALTER TABLE public.epitope_2010960_second_secreted_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_2010960_second_secreted_glycoprotein__cell_type
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__epi_an_start
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__epi_an_nstop
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__epi_frag_an_start
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__epi_frag_an_stop
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__host_tax_id
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__host_tax_name
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__iedb_id
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__is_linear
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__mhc_allele
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__mhc_class_lower
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__product_lower
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__response_freq_pos
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__seq_aa_alt
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__seq_aa_orig
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__start_aa_orig
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__taxon_id
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__taxon_name_lower
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__variant_aa_length
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__variant_aa_type
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_n_host_tax_id
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_host_cell_type
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_host_epi_start
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__virus_host_epi_stop
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__virus_host_is_linear
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_host_mhc_allele
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_host_product
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_second_secreted_glycoprotein__vir_host_resp_freq
    ON public.epitope_2010960_second_secreted_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT minor nucleoprotein
-- 2010960 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_minor_nucleoprotein (
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

ALTER TABLE public.epitope_2010960_minor_nucleoprotein
    OWNER TO geco;

CREATE INDEX epi_2010960_minor_nucleoprotein__cell_type
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__epi_an_start
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__epi_an_nstop
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__epi_frag_an_start
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__epi_frag_an_stop
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__host_tax_id
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__host_tax_name
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__iedb_id
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__is_linear
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__mhc_allele
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__mhc_class_lower
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__product_lower
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__response_freq_pos
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__seq_aa_alt
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__seq_aa_orig
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__start_aa_orig
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__taxon_id
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__taxon_name_lower
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__variant_aa_length
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__variant_aa_type
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_n_host_tax_id
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_host_cell_type
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_host_epi_start
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__virus_host_epi_stop
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__virus_host_is_linear
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_host_mhc_allele
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_host_product
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_minor_nucleoprotein__vir_host_resp_freq
    ON public.epitope_2010960_minor_nucleoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT membrane-associated protein
-- 2010960 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_membrane_associated_protein (
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

ALTER TABLE public.epitope_2010960_membrane_associated_protein
    OWNER TO geco;

CREATE INDEX epi_2010960_membrane_associated_protein__cell_type
    ON public.epitope_2010960_membrane_associated_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__epi_an_start
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__epi_an_nstop
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__epi_frag_an_start
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__epi_frag_an_stop
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__host_tax_id
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__host_tax_name
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__iedb_id
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__is_linear
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__mhc_allele
    ON public.epitope_2010960_membrane_associated_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__mhc_class_lower
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__product_lower
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__response_freq_pos
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__seq_aa_alt
    ON public.epitope_2010960_membrane_associated_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__seq_aa_orig
    ON public.epitope_2010960_membrane_associated_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__start_aa_orig
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__taxon_id
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__taxon_name_lower
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__variant_aa_length
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__variant_aa_type
    ON public.epitope_2010960_membrane_associated_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_n_host_tax_id
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_host_cell_type
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_host_epi_start
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__virus_host_epi_stop
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__virus_host_is_linear
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_host_mhc_allele
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_host_product
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_membrane_associated_protein__vir_host_resp_freq
    ON public.epitope_2010960_membrane_associated_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT RNA-dependent RNA polymerase
-- 2010960 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2010960_rna_dependent_rna_polymerase (
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

ALTER TABLE public.epitope_2010960_rna_dependent_rna_polymerase
    OWNER TO geco;

CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__cell_type
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__epi_an_start
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__epi_an_nstop
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__epi_frag_an_start
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__epi_frag_an_stop
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__host_tax_id
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__host_tax_name
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__iedb_id
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__is_linear
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__mhc_allele
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__mhc_class_lower
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__product_lower
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__response_freq_pos
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__seq_aa_alt
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__seq_aa_orig
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__start_aa_orig
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__taxon_id
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__taxon_name_lower
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__variant_aa_length
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__variant_aa_type
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_n_host_tax_id
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_host_cell_type
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_host_epi_start
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__virus_host_epi_stop
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__virus_host_is_linear
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_host_mhc_allele
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_host_product
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2010960_rna_dependent_rna_polymerase__vir_host_resp_freq
    ON public.epitope_2010960_rna_dependent_rna_polymerase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


