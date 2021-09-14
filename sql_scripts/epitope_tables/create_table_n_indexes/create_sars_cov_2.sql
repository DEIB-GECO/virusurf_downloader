-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1ab polyprotein
-- 2697049 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_orf1ab_polyprotein (
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

ALTER TABLE public.epitope_2697049_orf1ab_polyprotein
    OWNER TO geco;

CREATE INDEX epi_2697049_orf1ab_polyprotein__cell_type
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_an_start
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_an_nstop
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_frag_an_start
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_frag_an_stop
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__host_tax_id
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__host_tax_name
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__iedb_id
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__is_linear
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__mhc_allele
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__mhc_class_lower
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__product_lower
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__response_freq_pos
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__seq_aa_alt
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__seq_aa_orig
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__start_aa_orig
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__taxon_id
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__taxon_name_lower
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__variant_aa_length
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__variant_aa_type
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_n_host_tax_id
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_cell_type
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_epi_start
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__virus_host_epi_stop
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__virus_host_is_linear
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_mhc_allele
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_product
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_resp_freq
    ON public.epitope_2697049_orf1ab_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP12 (RNA-dependent RNA polymerase)
-- 2697049 can be replaced with the virus taxon id, while nsp12_rna_dependent_rna_poly can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp12_rna_dependent_rna_poly (
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

ALTER TABLE public.epitope_2697049_nsp12_rna_dependent_rna_poly
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__cell_type
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_start
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_nstop
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_start
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_stop
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_id
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_name
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__iedb_id
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__is_linear
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__mhc_allele
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__mhc_class_lower
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__product_lower
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__response_freq_pos
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_alt
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_orig
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__start_aa_orig
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__taxon_id
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__taxon_name_lower
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_length
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_type
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_n_host_tax_id
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_cell_type
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_epi_start
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_epi_stop
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_is_linear
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_mhc_allele
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_product
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_resp_freq
    ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP13 (helicase)
-- 2697049 can be replaced with the virus taxon id, while nsp13_helicase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp13_helicase (
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

ALTER TABLE public.epitope_2697049_nsp13_helicase
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp13_helicase__cell_type
    ON public.epitope_2697049_nsp13_helicase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__epi_an_start
    ON public.epitope_2697049_nsp13_helicase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__epi_an_nstop
    ON public.epitope_2697049_nsp13_helicase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__epi_frag_an_start
    ON public.epitope_2697049_nsp13_helicase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__epi_frag_an_stop
    ON public.epitope_2697049_nsp13_helicase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__host_tax_id
    ON public.epitope_2697049_nsp13_helicase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__host_tax_name
    ON public.epitope_2697049_nsp13_helicase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__iedb_id
    ON public.epitope_2697049_nsp13_helicase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__is_linear
    ON public.epitope_2697049_nsp13_helicase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__mhc_allele
    ON public.epitope_2697049_nsp13_helicase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__mhc_class_lower
    ON public.epitope_2697049_nsp13_helicase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__product_lower
    ON public.epitope_2697049_nsp13_helicase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__response_freq_pos
    ON public.epitope_2697049_nsp13_helicase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__seq_aa_alt
    ON public.epitope_2697049_nsp13_helicase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__seq_aa_orig
    ON public.epitope_2697049_nsp13_helicase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__start_aa_orig
    ON public.epitope_2697049_nsp13_helicase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__taxon_id
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__taxon_name_lower
    ON public.epitope_2697049_nsp13_helicase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__variant_aa_length
    ON public.epitope_2697049_nsp13_helicase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__variant_aa_type
    ON public.epitope_2697049_nsp13_helicase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_n_host_tax_id
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_host_cell_type
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_host_epi_start
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__virus_host_epi_stop
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__virus_host_is_linear
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_host_mhc_allele
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_host_product
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp13_helicase__vir_host_resp_freq
    ON public.epitope_2697049_nsp13_helicase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP14 (3'-to-5' exonuclease)
-- 2697049 can be replaced with the virus taxon id, while nsp14_3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp14_3_to_5_exonuclease (
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

ALTER TABLE public.epitope_2697049_nsp14_3_to_5_exonuclease
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__cell_type
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_an_start
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_an_nstop
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_start
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_stop
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__host_tax_id
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__host_tax_name
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__iedb_id
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__is_linear
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__mhc_allele
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__mhc_class_lower
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__product_lower
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__response_freq_pos
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_alt
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_orig
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__start_aa_orig
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__taxon_id
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__taxon_name_lower
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_length
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_type
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_n_host_tax_id
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_cell_type
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_epi_start
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__virus_host_epi_stop
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__virus_host_is_linear
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_mhc_allele
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_product
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_resp_freq
    ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP15 (endoRNAse)
-- 2697049 can be replaced with the virus taxon id, while nsp15_endornase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp15_endornase (
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

ALTER TABLE public.epitope_2697049_nsp15_endornase
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp15_endornase__cell_type
    ON public.epitope_2697049_nsp15_endornase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__epi_an_start
    ON public.epitope_2697049_nsp15_endornase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__epi_an_nstop
    ON public.epitope_2697049_nsp15_endornase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__epi_frag_an_start
    ON public.epitope_2697049_nsp15_endornase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__epi_frag_an_stop
    ON public.epitope_2697049_nsp15_endornase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__host_tax_id
    ON public.epitope_2697049_nsp15_endornase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__host_tax_name
    ON public.epitope_2697049_nsp15_endornase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__iedb_id
    ON public.epitope_2697049_nsp15_endornase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__is_linear
    ON public.epitope_2697049_nsp15_endornase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__mhc_allele
    ON public.epitope_2697049_nsp15_endornase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__mhc_class_lower
    ON public.epitope_2697049_nsp15_endornase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__product_lower
    ON public.epitope_2697049_nsp15_endornase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__response_freq_pos
    ON public.epitope_2697049_nsp15_endornase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__seq_aa_alt
    ON public.epitope_2697049_nsp15_endornase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__seq_aa_orig
    ON public.epitope_2697049_nsp15_endornase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__start_aa_orig
    ON public.epitope_2697049_nsp15_endornase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__taxon_id
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__taxon_name_lower
    ON public.epitope_2697049_nsp15_endornase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__variant_aa_length
    ON public.epitope_2697049_nsp15_endornase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__variant_aa_type
    ON public.epitope_2697049_nsp15_endornase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_n_host_tax_id
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_host_cell_type
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_host_epi_start
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__virus_host_epi_stop
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__virus_host_is_linear
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_host_mhc_allele
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_host_product
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp15_endornase__vir_host_resp_freq
    ON public.epitope_2697049_nsp15_endornase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP16 (2'-O-ribose methyltransferase)
-- 2697049 can be replaced with the virus taxon id, while nsp16_2_o_ribose_methyltrans can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp16_2_o_ribose_methyltrans (
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

ALTER TABLE public.epitope_2697049_nsp16_2_o_ribose_methyltrans
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__cell_type
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_start
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_nstop
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_start
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_stop
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_id
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_name
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__iedb_id
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__is_linear
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_allele
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_class_lower
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__product_lower
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__response_freq_pos
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_alt
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_orig
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__start_aa_orig
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_id
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_name_lower
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_length
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_type
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_n_host_tax_id
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_cell_type
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_epi_start
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_epi_stop
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_is_linear
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_mhc_allele
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_product
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_resp_freq
    ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1a polyprotein
-- 2697049 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_orf1a_polyprotein (
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

ALTER TABLE public.epitope_2697049_orf1a_polyprotein
    OWNER TO geco;

CREATE INDEX epi_2697049_orf1a_polyprotein__cell_type
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__epi_an_start
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__epi_an_nstop
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__epi_frag_an_start
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__epi_frag_an_stop
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__host_tax_id
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__host_tax_name
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__iedb_id
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__is_linear
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__mhc_allele
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__mhc_class_lower
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__product_lower
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__response_freq_pos
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__seq_aa_alt
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__seq_aa_orig
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__start_aa_orig
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__taxon_id
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__taxon_name_lower
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__variant_aa_length
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__variant_aa_type
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_n_host_tax_id
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_cell_type
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_epi_start
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__virus_host_epi_stop
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__virus_host_is_linear
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_mhc_allele
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_product
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_resp_freq
    ON public.epitope_2697049_orf1a_polyprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP1 (leader protein)
-- 2697049 can be replaced with the virus taxon id, while nsp1_leader_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp1_leader_protein (
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

ALTER TABLE public.epitope_2697049_nsp1_leader_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp1_leader_protein__cell_type
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__epi_an_start
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__epi_an_nstop
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__epi_frag_an_start
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__epi_frag_an_stop
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__host_tax_id
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__host_tax_name
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__iedb_id
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__is_linear
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__mhc_allele
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__mhc_class_lower
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__product_lower
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__response_freq_pos
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__seq_aa_alt
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__seq_aa_orig
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__start_aa_orig
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__taxon_id
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__taxon_name_lower
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__variant_aa_length
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__variant_aa_type
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_n_host_tax_id
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_cell_type
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_epi_start
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__virus_host_epi_stop
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__virus_host_is_linear
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_mhc_allele
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_product
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_resp_freq
    ON public.epitope_2697049_nsp1_leader_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP2
-- 2697049 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp2 (
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

ALTER TABLE public.epitope_2697049_nsp2
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp2__cell_type
    ON public.epitope_2697049_nsp2 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_an_start
    ON public.epitope_2697049_nsp2 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_an_nstop
    ON public.epitope_2697049_nsp2 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_frag_an_start
    ON public.epitope_2697049_nsp2 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__epi_frag_an_stop
    ON public.epitope_2697049_nsp2 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__host_tax_id
    ON public.epitope_2697049_nsp2 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__host_tax_name
    ON public.epitope_2697049_nsp2 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__iedb_id
    ON public.epitope_2697049_nsp2 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__is_linear
    ON public.epitope_2697049_nsp2 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__mhc_allele
    ON public.epitope_2697049_nsp2 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__mhc_class_lower
    ON public.epitope_2697049_nsp2 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__product_lower
    ON public.epitope_2697049_nsp2 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__response_freq_pos
    ON public.epitope_2697049_nsp2 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__seq_aa_alt
    ON public.epitope_2697049_nsp2 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__seq_aa_orig
    ON public.epitope_2697049_nsp2 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__start_aa_orig
    ON public.epitope_2697049_nsp2 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__taxon_id
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__taxon_name_lower
    ON public.epitope_2697049_nsp2 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__variant_aa_length
    ON public.epitope_2697049_nsp2 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__variant_aa_type
    ON public.epitope_2697049_nsp2 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_n_host_tax_id
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_host_cell_type
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_host_epi_start
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_epi_stop
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__virus_host_is_linear
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_host_mhc_allele
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_host_product
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp2__vir_host_resp_freq
    ON public.epitope_2697049_nsp2 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP3
-- 2697049 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp3 (
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

ALTER TABLE public.epitope_2697049_nsp3
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp3__cell_type
    ON public.epitope_2697049_nsp3 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_an_start
    ON public.epitope_2697049_nsp3 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_an_nstop
    ON public.epitope_2697049_nsp3 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_frag_an_start
    ON public.epitope_2697049_nsp3 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__epi_frag_an_stop
    ON public.epitope_2697049_nsp3 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__host_tax_id
    ON public.epitope_2697049_nsp3 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__host_tax_name
    ON public.epitope_2697049_nsp3 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__iedb_id
    ON public.epitope_2697049_nsp3 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__is_linear
    ON public.epitope_2697049_nsp3 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__mhc_allele
    ON public.epitope_2697049_nsp3 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__mhc_class_lower
    ON public.epitope_2697049_nsp3 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__product_lower
    ON public.epitope_2697049_nsp3 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__response_freq_pos
    ON public.epitope_2697049_nsp3 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__seq_aa_alt
    ON public.epitope_2697049_nsp3 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__seq_aa_orig
    ON public.epitope_2697049_nsp3 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__start_aa_orig
    ON public.epitope_2697049_nsp3 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__taxon_id
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__taxon_name_lower
    ON public.epitope_2697049_nsp3 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__variant_aa_length
    ON public.epitope_2697049_nsp3 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__variant_aa_type
    ON public.epitope_2697049_nsp3 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_n_host_tax_id
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_host_cell_type
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_host_epi_start
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_epi_stop
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__virus_host_is_linear
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_host_mhc_allele
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_host_product
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp3__vir_host_resp_freq
    ON public.epitope_2697049_nsp3 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP4
-- 2697049 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp4 (
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

ALTER TABLE public.epitope_2697049_nsp4
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp4__cell_type
    ON public.epitope_2697049_nsp4 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_an_start
    ON public.epitope_2697049_nsp4 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_an_nstop
    ON public.epitope_2697049_nsp4 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_frag_an_start
    ON public.epitope_2697049_nsp4 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__epi_frag_an_stop
    ON public.epitope_2697049_nsp4 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__host_tax_id
    ON public.epitope_2697049_nsp4 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__host_tax_name
    ON public.epitope_2697049_nsp4 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__iedb_id
    ON public.epitope_2697049_nsp4 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__is_linear
    ON public.epitope_2697049_nsp4 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__mhc_allele
    ON public.epitope_2697049_nsp4 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__mhc_class_lower
    ON public.epitope_2697049_nsp4 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__product_lower
    ON public.epitope_2697049_nsp4 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__response_freq_pos
    ON public.epitope_2697049_nsp4 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__seq_aa_alt
    ON public.epitope_2697049_nsp4 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__seq_aa_orig
    ON public.epitope_2697049_nsp4 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__start_aa_orig
    ON public.epitope_2697049_nsp4 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__taxon_id
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__taxon_name_lower
    ON public.epitope_2697049_nsp4 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__variant_aa_length
    ON public.epitope_2697049_nsp4 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__variant_aa_type
    ON public.epitope_2697049_nsp4 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_n_host_tax_id
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_host_cell_type
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_host_epi_start
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_epi_stop
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__virus_host_is_linear
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_host_mhc_allele
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_host_product
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp4__vir_host_resp_freq
    ON public.epitope_2697049_nsp4 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP5 (3C-like proteinase)
-- 2697049 can be replaced with the virus taxon id, while nsp5_3c_like_proteinase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp5_3c_like_proteinase (
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

ALTER TABLE public.epitope_2697049_nsp5_3c_like_proteinase
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__cell_type
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_an_start
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_an_nstop
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_start
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_stop
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__host_tax_id
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__host_tax_name
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__iedb_id
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__is_linear
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__mhc_allele
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__mhc_class_lower
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__product_lower
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__response_freq_pos
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__seq_aa_alt
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__seq_aa_orig
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__start_aa_orig
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__taxon_id
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__taxon_name_lower
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__variant_aa_length
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__variant_aa_type
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_n_host_tax_id
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_cell_type
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_epi_start
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__virus_host_epi_stop
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__virus_host_is_linear
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_mhc_allele
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_product
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_resp_freq
    ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP6
-- 2697049 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp6 (
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

ALTER TABLE public.epitope_2697049_nsp6
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp6__cell_type
    ON public.epitope_2697049_nsp6 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_an_start
    ON public.epitope_2697049_nsp6 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_an_nstop
    ON public.epitope_2697049_nsp6 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_frag_an_start
    ON public.epitope_2697049_nsp6 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__epi_frag_an_stop
    ON public.epitope_2697049_nsp6 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__host_tax_id
    ON public.epitope_2697049_nsp6 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__host_tax_name
    ON public.epitope_2697049_nsp6 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__iedb_id
    ON public.epitope_2697049_nsp6 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__is_linear
    ON public.epitope_2697049_nsp6 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__mhc_allele
    ON public.epitope_2697049_nsp6 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__mhc_class_lower
    ON public.epitope_2697049_nsp6 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__product_lower
    ON public.epitope_2697049_nsp6 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__response_freq_pos
    ON public.epitope_2697049_nsp6 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__seq_aa_alt
    ON public.epitope_2697049_nsp6 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__seq_aa_orig
    ON public.epitope_2697049_nsp6 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__start_aa_orig
    ON public.epitope_2697049_nsp6 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__taxon_id
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__taxon_name_lower
    ON public.epitope_2697049_nsp6 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__variant_aa_length
    ON public.epitope_2697049_nsp6 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__variant_aa_type
    ON public.epitope_2697049_nsp6 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_n_host_tax_id
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_host_cell_type
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_host_epi_start
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_epi_stop
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__virus_host_is_linear
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_host_mhc_allele
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_host_product
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp6__vir_host_resp_freq
    ON public.epitope_2697049_nsp6 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP7
-- 2697049 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp7 (
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

ALTER TABLE public.epitope_2697049_nsp7
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp7__cell_type
    ON public.epitope_2697049_nsp7 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_an_start
    ON public.epitope_2697049_nsp7 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_an_nstop
    ON public.epitope_2697049_nsp7 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_frag_an_start
    ON public.epitope_2697049_nsp7 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__epi_frag_an_stop
    ON public.epitope_2697049_nsp7 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__host_tax_id
    ON public.epitope_2697049_nsp7 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__host_tax_name
    ON public.epitope_2697049_nsp7 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__iedb_id
    ON public.epitope_2697049_nsp7 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__is_linear
    ON public.epitope_2697049_nsp7 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__mhc_allele
    ON public.epitope_2697049_nsp7 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__mhc_class_lower
    ON public.epitope_2697049_nsp7 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__product_lower
    ON public.epitope_2697049_nsp7 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__response_freq_pos
    ON public.epitope_2697049_nsp7 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__seq_aa_alt
    ON public.epitope_2697049_nsp7 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__seq_aa_orig
    ON public.epitope_2697049_nsp7 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__start_aa_orig
    ON public.epitope_2697049_nsp7 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__taxon_id
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__taxon_name_lower
    ON public.epitope_2697049_nsp7 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__variant_aa_length
    ON public.epitope_2697049_nsp7 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__variant_aa_type
    ON public.epitope_2697049_nsp7 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_n_host_tax_id
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_host_cell_type
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_host_epi_start
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_epi_stop
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__virus_host_is_linear
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_host_mhc_allele
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_host_product
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp7__vir_host_resp_freq
    ON public.epitope_2697049_nsp7 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP8
-- 2697049 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp8 (
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

ALTER TABLE public.epitope_2697049_nsp8
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp8__cell_type
    ON public.epitope_2697049_nsp8 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_an_start
    ON public.epitope_2697049_nsp8 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_an_nstop
    ON public.epitope_2697049_nsp8 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_frag_an_start
    ON public.epitope_2697049_nsp8 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__epi_frag_an_stop
    ON public.epitope_2697049_nsp8 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__host_tax_id
    ON public.epitope_2697049_nsp8 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__host_tax_name
    ON public.epitope_2697049_nsp8 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__iedb_id
    ON public.epitope_2697049_nsp8 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__is_linear
    ON public.epitope_2697049_nsp8 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__mhc_allele
    ON public.epitope_2697049_nsp8 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__mhc_class_lower
    ON public.epitope_2697049_nsp8 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__product_lower
    ON public.epitope_2697049_nsp8 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__response_freq_pos
    ON public.epitope_2697049_nsp8 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__seq_aa_alt
    ON public.epitope_2697049_nsp8 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__seq_aa_orig
    ON public.epitope_2697049_nsp8 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__start_aa_orig
    ON public.epitope_2697049_nsp8 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__taxon_id
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__taxon_name_lower
    ON public.epitope_2697049_nsp8 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__variant_aa_length
    ON public.epitope_2697049_nsp8 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__variant_aa_type
    ON public.epitope_2697049_nsp8 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_n_host_tax_id
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_host_cell_type
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_host_epi_start
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_epi_stop
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__virus_host_is_linear
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_host_mhc_allele
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_host_product
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp8__vir_host_resp_freq
    ON public.epitope_2697049_nsp8 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP9
-- 2697049 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp9 (
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

ALTER TABLE public.epitope_2697049_nsp9
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp9__cell_type
    ON public.epitope_2697049_nsp9 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_an_start
    ON public.epitope_2697049_nsp9 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_an_nstop
    ON public.epitope_2697049_nsp9 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_frag_an_start
    ON public.epitope_2697049_nsp9 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__epi_frag_an_stop
    ON public.epitope_2697049_nsp9 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__host_tax_id
    ON public.epitope_2697049_nsp9 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__host_tax_name
    ON public.epitope_2697049_nsp9 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__iedb_id
    ON public.epitope_2697049_nsp9 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__is_linear
    ON public.epitope_2697049_nsp9 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__mhc_allele
    ON public.epitope_2697049_nsp9 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__mhc_class_lower
    ON public.epitope_2697049_nsp9 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__product_lower
    ON public.epitope_2697049_nsp9 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__response_freq_pos
    ON public.epitope_2697049_nsp9 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__seq_aa_alt
    ON public.epitope_2697049_nsp9 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__seq_aa_orig
    ON public.epitope_2697049_nsp9 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__start_aa_orig
    ON public.epitope_2697049_nsp9 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__taxon_id
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__taxon_name_lower
    ON public.epitope_2697049_nsp9 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__variant_aa_length
    ON public.epitope_2697049_nsp9 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__variant_aa_type
    ON public.epitope_2697049_nsp9 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_n_host_tax_id
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_host_cell_type
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_host_epi_start
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_epi_stop
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__virus_host_is_linear
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_host_mhc_allele
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_host_product
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp9__vir_host_resp_freq
    ON public.epitope_2697049_nsp9 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP10
-- 2697049 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp10 (
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

ALTER TABLE public.epitope_2697049_nsp10
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp10__cell_type
    ON public.epitope_2697049_nsp10 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_an_start
    ON public.epitope_2697049_nsp10 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_an_nstop
    ON public.epitope_2697049_nsp10 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_frag_an_start
    ON public.epitope_2697049_nsp10 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__epi_frag_an_stop
    ON public.epitope_2697049_nsp10 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__host_tax_id
    ON public.epitope_2697049_nsp10 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__host_tax_name
    ON public.epitope_2697049_nsp10 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__iedb_id
    ON public.epitope_2697049_nsp10 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__is_linear
    ON public.epitope_2697049_nsp10 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__mhc_allele
    ON public.epitope_2697049_nsp10 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__mhc_class_lower
    ON public.epitope_2697049_nsp10 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__product_lower
    ON public.epitope_2697049_nsp10 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__response_freq_pos
    ON public.epitope_2697049_nsp10 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__seq_aa_alt
    ON public.epitope_2697049_nsp10 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__seq_aa_orig
    ON public.epitope_2697049_nsp10 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__start_aa_orig
    ON public.epitope_2697049_nsp10 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__taxon_id
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__taxon_name_lower
    ON public.epitope_2697049_nsp10 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__variant_aa_length
    ON public.epitope_2697049_nsp10 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__variant_aa_type
    ON public.epitope_2697049_nsp10 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_n_host_tax_id
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_host_cell_type
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_host_epi_start
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_epi_stop
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__virus_host_is_linear
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_host_mhc_allele
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_host_product
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp10__vir_host_resp_freq
    ON public.epitope_2697049_nsp10 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP11
-- 2697049 can be replaced with the virus taxon id, while nsp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_nsp11 (
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

ALTER TABLE public.epitope_2697049_nsp11
    OWNER TO geco;

CREATE INDEX epi_2697049_nsp11__cell_type
    ON public.epitope_2697049_nsp11 USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_an_start
    ON public.epitope_2697049_nsp11 USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_an_nstop
    ON public.epitope_2697049_nsp11 USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_frag_an_start
    ON public.epitope_2697049_nsp11 USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__epi_frag_an_stop
    ON public.epitope_2697049_nsp11 USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__host_tax_id
    ON public.epitope_2697049_nsp11 USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__host_tax_name
    ON public.epitope_2697049_nsp11 USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__iedb_id
    ON public.epitope_2697049_nsp11 USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__is_linear
    ON public.epitope_2697049_nsp11 USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__mhc_allele
    ON public.epitope_2697049_nsp11 USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__mhc_class_lower
    ON public.epitope_2697049_nsp11 USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__product_lower
    ON public.epitope_2697049_nsp11 USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__response_freq_pos
    ON public.epitope_2697049_nsp11 USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__seq_aa_alt
    ON public.epitope_2697049_nsp11 USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__seq_aa_orig
    ON public.epitope_2697049_nsp11 USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__start_aa_orig
    ON public.epitope_2697049_nsp11 USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__taxon_id
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__taxon_name_lower
    ON public.epitope_2697049_nsp11 USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__variant_aa_length
    ON public.epitope_2697049_nsp11 USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__variant_aa_type
    ON public.epitope_2697049_nsp11 USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_n_host_tax_id
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_host_cell_type
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_host_epi_start
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_epi_stop
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__virus_host_is_linear
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_host_mhc_allele
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_host_product
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_nsp11__vir_host_resp_freq
    ON public.epitope_2697049_nsp11 USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT Spike (surface glycoprotein)
-- 2697049 can be replaced with the virus taxon id, while spike_surface_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_spike_surface_glycoprotein (
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

ALTER TABLE public.epitope_2697049_spike_surface_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_2697049_spike_surface_glycoprotein__cell_type
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_an_start
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_an_nstop
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_frag_an_start
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_frag_an_stop
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__host_tax_id
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__host_tax_name
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__iedb_id
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__is_linear
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__mhc_allele
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__mhc_class_lower
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__product_lower
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__response_freq_pos
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__seq_aa_alt
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__seq_aa_orig
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__start_aa_orig
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__taxon_id
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__taxon_name_lower
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__variant_aa_length
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__variant_aa_type
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_n_host_tax_id
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_cell_type
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_epi_start
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__virus_host_epi_stop
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__virus_host_is_linear
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_mhc_allele
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_product
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_resp_freq
    ON public.epitope_2697049_spike_surface_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS3 (ORF3a protein)
-- 2697049 can be replaced with the virus taxon id, while ns3_orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_ns3_orf3a_protein (
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

ALTER TABLE public.epitope_2697049_ns3_orf3a_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_ns3_orf3a_protein__cell_type
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_an_start
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_an_nstop
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_frag_an_start
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_frag_an_stop
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__host_tax_id
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__host_tax_name
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__iedb_id
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__is_linear
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__mhc_allele
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__mhc_class_lower
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__product_lower
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__response_freq_pos
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__seq_aa_alt
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__seq_aa_orig
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__start_aa_orig
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__taxon_id
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__taxon_name_lower
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__variant_aa_length
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__variant_aa_type
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_n_host_tax_id
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_cell_type
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_epi_start
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__virus_host_epi_stop
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__virus_host_is_linear
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_mhc_allele
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_product
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_resp_freq
    ON public.epitope_2697049_ns3_orf3a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT E (envelope protein)
-- 2697049 can be replaced with the virus taxon id, while e_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_e_envelope_protein (
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

ALTER TABLE public.epitope_2697049_e_envelope_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_e_envelope_protein__cell_type
    ON public.epitope_2697049_e_envelope_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__epi_an_start
    ON public.epitope_2697049_e_envelope_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__epi_an_nstop
    ON public.epitope_2697049_e_envelope_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__epi_frag_an_start
    ON public.epitope_2697049_e_envelope_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__epi_frag_an_stop
    ON public.epitope_2697049_e_envelope_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__host_tax_id
    ON public.epitope_2697049_e_envelope_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__host_tax_name
    ON public.epitope_2697049_e_envelope_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__iedb_id
    ON public.epitope_2697049_e_envelope_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__is_linear
    ON public.epitope_2697049_e_envelope_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__mhc_allele
    ON public.epitope_2697049_e_envelope_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__mhc_class_lower
    ON public.epitope_2697049_e_envelope_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__product_lower
    ON public.epitope_2697049_e_envelope_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__response_freq_pos
    ON public.epitope_2697049_e_envelope_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__seq_aa_alt
    ON public.epitope_2697049_e_envelope_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__seq_aa_orig
    ON public.epitope_2697049_e_envelope_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__start_aa_orig
    ON public.epitope_2697049_e_envelope_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__taxon_id
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__taxon_name_lower
    ON public.epitope_2697049_e_envelope_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__variant_aa_length
    ON public.epitope_2697049_e_envelope_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__variant_aa_type
    ON public.epitope_2697049_e_envelope_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_n_host_tax_id
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_host_cell_type
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_host_epi_start
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__virus_host_epi_stop
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__virus_host_is_linear
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_host_mhc_allele
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_host_product
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_e_envelope_protein__vir_host_resp_freq
    ON public.epitope_2697049_e_envelope_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT M (membrane glycoprotein)
-- 2697049 can be replaced with the virus taxon id, while m_membrane_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_m_membrane_glycoprotein (
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

ALTER TABLE public.epitope_2697049_m_membrane_glycoprotein
    OWNER TO geco;

CREATE INDEX epi_2697049_m_membrane_glycoprotein__cell_type
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_an_start
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_an_nstop
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_frag_an_start
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_frag_an_stop
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__host_tax_id
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__host_tax_name
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__iedb_id
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__is_linear
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__mhc_allele
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__mhc_class_lower
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__product_lower
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__response_freq_pos
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__seq_aa_alt
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__seq_aa_orig
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__start_aa_orig
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__taxon_id
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__taxon_name_lower
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__variant_aa_length
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__variant_aa_type
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_n_host_tax_id
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_cell_type
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_epi_start
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__virus_host_epi_stop
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__virus_host_is_linear
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_mhc_allele
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_product
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_resp_freq
    ON public.epitope_2697049_m_membrane_glycoprotein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS6 (ORF6 protein)
-- 2697049 can be replaced with the virus taxon id, while ns6_orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_ns6_orf6_protein (
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

ALTER TABLE public.epitope_2697049_ns6_orf6_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_ns6_orf6_protein__cell_type
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__epi_an_start
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__epi_an_nstop
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__epi_frag_an_start
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__epi_frag_an_stop
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__host_tax_id
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__host_tax_name
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__iedb_id
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__is_linear
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__mhc_allele
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__mhc_class_lower
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__product_lower
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__response_freq_pos
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__seq_aa_alt
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__seq_aa_orig
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__start_aa_orig
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__taxon_id
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__taxon_name_lower
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__variant_aa_length
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__variant_aa_type
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_n_host_tax_id
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_cell_type
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_epi_start
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__virus_host_epi_stop
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__virus_host_is_linear
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_mhc_allele
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_product
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_resp_freq
    ON public.epitope_2697049_ns6_orf6_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7a (ORF7a protein)
-- 2697049 can be replaced with the virus taxon id, while ns7a_orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_ns7a_orf7a_protein (
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

ALTER TABLE public.epitope_2697049_ns7a_orf7a_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_ns7a_orf7a_protein__cell_type
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_an_start
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_an_nstop
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_frag_an_start
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_frag_an_stop
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__host_tax_id
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__host_tax_name
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__iedb_id
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__is_linear
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__mhc_allele
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__mhc_class_lower
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__product_lower
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__response_freq_pos
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__seq_aa_alt
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__seq_aa_orig
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__start_aa_orig
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__taxon_id
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__taxon_name_lower
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__variant_aa_length
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__variant_aa_type
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_n_host_tax_id
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_cell_type
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_epi_start
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__virus_host_epi_stop
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__virus_host_is_linear
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_mhc_allele
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_product
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_resp_freq
    ON public.epitope_2697049_ns7a_orf7a_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7b (ORF7b)
-- 2697049 can be replaced with the virus taxon id, while ns7b_orf7b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_ns7b_orf7b (
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

ALTER TABLE public.epitope_2697049_ns7b_orf7b
    OWNER TO geco;

CREATE INDEX epi_2697049_ns7b_orf7b__cell_type
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_an_start
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_an_nstop
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_an_start
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_an_stop
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__host_tax_id
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__host_tax_name
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__iedb_id
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__is_linear
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__mhc_allele
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__mhc_class_lower
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__product_lower
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__response_freq_pos
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__seq_aa_alt
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__seq_aa_orig
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__start_aa_orig
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__taxon_id
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__taxon_name_lower
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_length
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_type
    ON public.epitope_2697049_ns7b_orf7b USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_n_host_tax_id
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_cell_type
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_epi_start
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_stop
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_is_linear
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_mhc_allele
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_product
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_resp_freq
    ON public.epitope_2697049_ns7b_orf7b USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS8 (ORF8 protein)
-- 2697049 can be replaced with the virus taxon id, while ns8_orf8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_ns8_orf8_protein (
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

ALTER TABLE public.epitope_2697049_ns8_orf8_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_ns8_orf8_protein__cell_type
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__epi_an_start
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__epi_an_nstop
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__epi_frag_an_start
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__epi_frag_an_stop
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__host_tax_id
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__host_tax_name
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__iedb_id
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__is_linear
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__mhc_allele
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__mhc_class_lower
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__product_lower
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__response_freq_pos
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__seq_aa_alt
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__seq_aa_orig
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__start_aa_orig
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__taxon_id
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__taxon_name_lower
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__variant_aa_length
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__variant_aa_type
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_n_host_tax_id
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_cell_type
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_epi_start
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__virus_host_epi_stop
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__virus_host_is_linear
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_mhc_allele
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_product
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_resp_freq
    ON public.epitope_2697049_ns8_orf8_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT N (nucleocapsid phosphoprotein)
-- 2697049 can be replaced with the virus taxon id, while n_nucleocapsid_phosphoprotei can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_n_nucleocapsid_phosphoprotei (
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

ALTER TABLE public.epitope_2697049_n_nucleocapsid_phosphoprotei
    OWNER TO geco;

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__cell_type
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_start
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_nstop
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_start
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_stop
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_id
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_name
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__iedb_id
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__is_linear
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__mhc_allele
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__mhc_class_lower
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__product_lower
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__response_freq_pos
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_alt
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_orig
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__start_aa_orig
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__taxon_id
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__taxon_name_lower
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_length
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_type
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_n_host_tax_id
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_cell_type
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_epi_start
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_epi_stop
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_is_linear
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_mhc_allele
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_product
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_resp_freq
    ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- CREATE TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF10 protein
-- 2697049 can be replaced with the virus taxon id, while orf10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
CREATE TABLE public.epitope_2697049_orf10_protein (
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

ALTER TABLE public.epitope_2697049_orf10_protein
    OWNER TO geco;

CREATE INDEX epi_2697049_orf10_protein__cell_type
    ON public.epitope_2697049_orf10_protein USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__epi_an_start
    ON public.epitope_2697049_orf10_protein USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__epi_an_nstop
    ON public.epitope_2697049_orf10_protein USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__epi_frag_an_start
    ON public.epitope_2697049_orf10_protein USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__epi_frag_an_stop
    ON public.epitope_2697049_orf10_protein USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__host_tax_id
    ON public.epitope_2697049_orf10_protein USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__host_tax_name
    ON public.epitope_2697049_orf10_protein USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__iedb_id
    ON public.epitope_2697049_orf10_protein USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__is_linear
    ON public.epitope_2697049_orf10_protein USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__mhc_allele
    ON public.epitope_2697049_orf10_protein USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__mhc_class_lower
    ON public.epitope_2697049_orf10_protein USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__product_lower
    ON public.epitope_2697049_orf10_protein USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__response_freq_pos
    ON public.epitope_2697049_orf10_protein USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__seq_aa_alt
    ON public.epitope_2697049_orf10_protein USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__seq_aa_orig
    ON public.epitope_2697049_orf10_protein USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__start_aa_orig
    ON public.epitope_2697049_orf10_protein USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__taxon_id
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__taxon_name_lower
    ON public.epitope_2697049_orf10_protein USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__variant_aa_length
    ON public.epitope_2697049_orf10_protein USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__variant_aa_type
    ON public.epitope_2697049_orf10_protein USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_n_host_tax_id
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_host_cell_type
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_host_epi_start
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__virus_host_epi_stop
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__virus_host_is_linear
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_host_mhc_allele
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_host_product
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_2697049_orf10_protein__vir_host_resp_freq
    ON public.epitope_2697049_orf10_protein USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


