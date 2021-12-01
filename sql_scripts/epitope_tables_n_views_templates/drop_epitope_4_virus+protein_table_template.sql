-- $virus_id can be replaced with the virus taxon id, while $short_prot_name can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_$virus_id_$short_prot_name;

DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__cell_type;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__epi_an_start;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__epi_an_nstop;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__epi_frag_an_start;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__host_tax_id;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__host_tax_name;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__iedb_id;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__is_linear;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__mhc_allele;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__mhc_class_lower;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__product_lower;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__response_freq_pos;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__seq_aa_alt;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__seq_aa_orig;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__start_aa_orig;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__taxon_id;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__taxon_name_lower;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__variant_aa_length;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__variant_aa_type;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_host_cell_type;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_host_epi_start;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__virus_host_is_linear;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_host_product;
DROP INDEX IF EXISTS epi_$virus_id_$short_prot_name__vir_host_resp_freq;
