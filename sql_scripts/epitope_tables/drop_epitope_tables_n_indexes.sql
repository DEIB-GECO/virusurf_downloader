-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1ab polyprotein
-- 2697049 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_orf1ab_polyprotein;

DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_orf1ab_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP12 (RNA-dependent RNA polymerase)
-- 2697049 can be replaced with the virus taxon id, while nsp12_rna_dependent_rna_poly can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp12_rna_dependent_rna_poly;

DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP13 (helicase)
-- 2697049 can be replaced with the virus taxon id, while nsp13_helicase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp13_helicase;

DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp13_helicase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP14 (3'-to-5' exonuclease)
-- 2697049 can be replaced with the virus taxon id, while nsp14_3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp14_3_to_5_exonuclease;

DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp14_3_to_5_exonuclease__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP15 (endoRNAse)
-- 2697049 can be replaced with the virus taxon id, while nsp15_endornase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp15_endornase;

DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp15_endornase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP16 (2'-O-ribose methyltransferase)
-- 2697049 can be replaced with the virus taxon id, while nsp16_2_o_ribose_methyltrans can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp16_2_o_ribose_methyltrans;

DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF1a polyprotein
-- 2697049 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_orf1a_polyprotein;

DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_orf1a_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP1 (leader protein)
-- 2697049 can be replaced with the virus taxon id, while nsp1_leader_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp1_leader_protein;

DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp1_leader_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP2
-- 2697049 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp2;

DROP INDEX IF EXISTS epi_2697049_nsp2__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp2__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp2__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp2__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp2__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp2__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp2__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp2__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp2__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp2__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp2__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp2__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp2__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp2__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp2__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp2__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp2__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp2__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp2__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp2__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp2__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp2__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp2__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP3
-- 2697049 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp3;

DROP INDEX IF EXISTS epi_2697049_nsp3__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp3__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp3__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp3__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp3__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp3__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp3__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp3__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp3__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp3__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp3__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp3__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp3__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp3__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp3__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp3__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp3__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp3__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP4
-- 2697049 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp4;

DROP INDEX IF EXISTS epi_2697049_nsp4__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp4__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp4__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp4__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp4__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp4__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp4__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp4__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp4__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp4__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp4__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp4__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp4__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp4__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp4__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp4__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp4__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp4__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp4__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp4__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp4__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp4__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp4__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP5 (3C-like proteinase)
-- 2697049 can be replaced with the virus taxon id, while nsp5_3c_like_proteinase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp5_3c_like_proteinase;

DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp5_3c_like_proteinase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP6
-- 2697049 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp6;

DROP INDEX IF EXISTS epi_2697049_nsp6__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp6__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp6__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp6__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp6__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp6__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp6__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp6__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp6__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp6__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp6__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp6__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp6__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp6__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp6__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp6__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp6__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp6__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp6__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp6__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp6__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp6__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp6__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP7
-- 2697049 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp7;

DROP INDEX IF EXISTS epi_2697049_nsp7__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp7__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp7__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp7__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp7__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp7__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp7__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp7__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp7__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp7__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp7__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp7__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp7__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp7__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp7__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp7__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp7__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp7__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp7__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp7__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp7__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp7__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp7__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP8
-- 2697049 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp8;

DROP INDEX IF EXISTS epi_2697049_nsp8__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp8__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp8__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp8__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp8__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp8__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp8__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp8__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp8__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp8__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp8__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp8__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp8__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp8__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp8__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp8__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp8__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp8__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp8__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp8__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp8__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp8__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp8__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP9
-- 2697049 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp9;

DROP INDEX IF EXISTS epi_2697049_nsp9__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp9__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp9__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp9__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp9__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp9__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp9__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp9__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp9__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp9__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp9__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp9__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp9__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp9__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp9__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp9__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp9__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp9__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp9__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp9__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp9__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp9__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp9__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP10
-- 2697049 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp10;

DROP INDEX IF EXISTS epi_2697049_nsp10__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp10__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp10__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp10__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp10__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp10__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp10__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp10__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp10__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp10__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp10__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp10__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp10__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp10__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp10__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp10__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp10__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp10__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp10__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp10__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp10__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp10__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp10__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NSP11
-- 2697049 can be replaced with the virus taxon id, while nsp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_nsp11;

DROP INDEX IF EXISTS epi_2697049_nsp11__cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp11__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp11__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_nsp11__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_nsp11__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_nsp11__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp11__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_nsp11__iedb_id;
DROP INDEX IF EXISTS epi_2697049_nsp11__is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp11__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp11__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_nsp11__product_lower;
DROP INDEX IF EXISTS epi_2697049_nsp11__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_nsp11__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_nsp11__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp11__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_nsp11__taxon_id;
DROP INDEX IF EXISTS epi_2697049_nsp11__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_nsp11__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_nsp11__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_nsp11__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_nsp11__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_nsp11__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT Spike (surface glycoprotein)
-- 2697049 can be replaced with the virus taxon id, while spike_surface_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_spike_surface_glycoprotein;

DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_spike_surface_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS3 (ORF3a protein)
-- 2697049 can be replaced with the virus taxon id, while ns3_orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_ns3_orf3a_protein;

DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_ns3_orf3a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT E (envelope protein)
-- 2697049 can be replaced with the virus taxon id, while e_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_e_envelope_protein;

DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_e_envelope_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT M (membrane glycoprotein)
-- 2697049 can be replaced with the virus taxon id, while m_membrane_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_m_membrane_glycoprotein;

DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_m_membrane_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS6 (ORF6 protein)
-- 2697049 can be replaced with the virus taxon id, while ns6_orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_ns6_orf6_protein;

DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_ns6_orf6_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7a (ORF7a protein)
-- 2697049 can be replaced with the virus taxon id, while ns7a_orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_ns7a_orf7a_protein;

DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_ns7a_orf7a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS7b (ORF7b)
-- 2697049 can be replaced with the virus taxon id, while ns7b_orf7b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_ns7b_orf7b;

DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__cell_type;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__iedb_id;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__is_linear;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__product_lower;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__taxon_id;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_ns7b_orf7b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT NS8 (ORF8 protein)
-- 2697049 can be replaced with the virus taxon id, while ns8_orf8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_ns8_orf8_protein;

DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_ns8_orf8_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT N (nucleocapsid phosphoprotein)
-- 2697049 can be replaced with the virus taxon id, while n_nucleocapsid_phosphoprotei can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_n_nucleocapsid_phosphoprotei;

DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__cell_type;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__iedb_id;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__is_linear;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__product_lower;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__taxon_id;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_2 and PROT ORF10 protein
-- 2697049 can be replaced with the virus taxon id, while orf10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2697049_orf10_protein;

DROP INDEX IF EXISTS epi_2697049_orf10_protein__cell_type;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__iedb_id;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__is_linear;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__product_lower;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__taxon_id;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2697049_orf10_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1ab polyprotein
-- 694009 can be replaced with the virus taxon id, while orf1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf1ab_polyprotein;

DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf1ab_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT RNA-dependent RNA polymerase
-- 694009 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_694009_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT helicase/NTPase
-- 694009 can be replaced with the virus taxon id, while helicase_ntpase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_helicase_ntpase;

DROP INDEX IF EXISTS epi_694009_helicase_ntpase__cell_type;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__epi_an_start;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__host_tax_id;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__host_tax_name;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__iedb_id;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__is_linear;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__mhc_allele;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__product_lower;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__taxon_id;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_host_product;
DROP INDEX IF EXISTS epi_694009_helicase_ntpase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3' to 5' exonuclease
-- 694009 can be replaced with the virus taxon id, while 3_to_5_exonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_3_to_5_exonuclease;

DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__cell_type;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__epi_an_start;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__host_tax_id;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__host_tax_name;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__iedb_id;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__is_linear;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__mhc_allele;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__product_lower;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__taxon_id;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_host_product;
DROP INDEX IF EXISTS epi_694009_3_to_5_exonuclease__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT endoribonuclease
-- 694009 can be replaced with the virus taxon id, while endoribonuclease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_endoribonuclease;

DROP INDEX IF EXISTS epi_694009_endoribonuclease__cell_type;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__epi_an_start;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__host_tax_id;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__host_tax_name;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__iedb_id;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__is_linear;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__mhc_allele;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__product_lower;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__taxon_id;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_host_product;
DROP INDEX IF EXISTS epi_694009_endoribonuclease__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 2'-O-MTase
-- 694009 can be replaced with the virus taxon id, while 2_o_mtase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_2_o_mtase;

DROP INDEX IF EXISTS epi_694009_2_o_mtase__cell_type;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__epi_an_start;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__host_tax_id;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__host_tax_name;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__iedb_id;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__is_linear;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__mhc_allele;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__product_lower;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__taxon_id;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_host_product;
DROP INDEX IF EXISTS epi_694009_2_o_mtase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF1a polyprotein
-- 694009 can be replaced with the virus taxon id, while orf1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf1a_polyprotein;

DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf1a_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp1
-- 694009 can be replaced with the virus taxon id, while nsp1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp1;

DROP INDEX IF EXISTS epi_694009_nsp1__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp1__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp1__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp1__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp1__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp1__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp1__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp1__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp1__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp1__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp1__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp1__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp1__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp1__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp1__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp1__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp1__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp1__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp1__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp1__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp1__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp1__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp1__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp2
-- 694009 can be replaced with the virus taxon id, while nsp2 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp2;

DROP INDEX IF EXISTS epi_694009_nsp2__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp2__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp2__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp2__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp2__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp2__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp2__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp2__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp2__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp2__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp2__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp2__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp2__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp2__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp2__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp2__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp2__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp2__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp2__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp2__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp2__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp2__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp2__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp3
-- 694009 can be replaced with the virus taxon id, while nsp3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp3;

DROP INDEX IF EXISTS epi_694009_nsp3__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp3__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp3__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp3__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp3__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp3__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp3__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp3__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp3__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp3__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp3__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp3__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp3__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp3__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp3__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp3__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp3__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp3__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp4
-- 694009 can be replaced with the virus taxon id, while nsp4 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp4;

DROP INDEX IF EXISTS epi_694009_nsp4__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp4__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp4__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp4__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp4__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp4__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp4__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp4__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp4__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp4__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp4__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp4__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp4__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp4__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp4__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp4__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp4__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp4__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp4__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp4__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp4__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp4__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp4__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT 3C-like protease
-- 694009 can be replaced with the virus taxon id, while 3c_like_protease can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_3c_like_protease;

DROP INDEX IF EXISTS epi_694009_3c_like_protease__cell_type;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__epi_an_start;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__host_tax_id;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__host_tax_name;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__iedb_id;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__is_linear;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__mhc_allele;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__product_lower;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__taxon_id;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_host_product;
DROP INDEX IF EXISTS epi_694009_3c_like_protease__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp6
-- 694009 can be replaced with the virus taxon id, while nsp6 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp6;

DROP INDEX IF EXISTS epi_694009_nsp6__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp6__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp6__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp6__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp6__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp6__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp6__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp6__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp6__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp6__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp6__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp6__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp6__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp6__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp6__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp6__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp6__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp6__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp6__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp6__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp6__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp6__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp6__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp7
-- 694009 can be replaced with the virus taxon id, while nsp7 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp7;

DROP INDEX IF EXISTS epi_694009_nsp7__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp7__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp7__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp7__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp7__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp7__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp7__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp7__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp7__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp7__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp7__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp7__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp7__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp7__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp7__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp7__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp7__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp7__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp7__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp7__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp7__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp7__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp7__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp8
-- 694009 can be replaced with the virus taxon id, while nsp8 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp8;

DROP INDEX IF EXISTS epi_694009_nsp8__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp8__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp8__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp8__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp8__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp8__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp8__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp8__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp8__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp8__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp8__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp8__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp8__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp8__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp8__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp8__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp8__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp8__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp8__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp8__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp8__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp8__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp8__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp9
-- 694009 can be replaced with the virus taxon id, while nsp9 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp9;

DROP INDEX IF EXISTS epi_694009_nsp9__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp9__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp9__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp9__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp9__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp9__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp9__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp9__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp9__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp9__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp9__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp9__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp9__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp9__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp9__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp9__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp9__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp9__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp9__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp9__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp9__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp9__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp9__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nsp10
-- 694009 can be replaced with the virus taxon id, while nsp10 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nsp10;

DROP INDEX IF EXISTS epi_694009_nsp10__cell_type;
DROP INDEX IF EXISTS epi_694009_nsp10__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nsp10__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nsp10__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nsp10__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nsp10__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp10__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nsp10__iedb_id;
DROP INDEX IF EXISTS epi_694009_nsp10__is_linear;
DROP INDEX IF EXISTS epi_694009_nsp10__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp10__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nsp10__product_lower;
DROP INDEX IF EXISTS epi_694009_nsp10__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nsp10__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nsp10__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp10__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nsp10__taxon_id;
DROP INDEX IF EXISTS epi_694009_nsp10__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nsp10__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nsp10__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nsp10__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nsp10__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nsp10__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ndp11
-- 694009 can be replaced with the virus taxon id, while ndp11 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_ndp11;

DROP INDEX IF EXISTS epi_694009_ndp11__cell_type;
DROP INDEX IF EXISTS epi_694009_ndp11__epi_an_start;
DROP INDEX IF EXISTS epi_694009_ndp11__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_ndp11__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_ndp11__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_ndp11__host_tax_id;
DROP INDEX IF EXISTS epi_694009_ndp11__host_tax_name;
DROP INDEX IF EXISTS epi_694009_ndp11__iedb_id;
DROP INDEX IF EXISTS epi_694009_ndp11__is_linear;
DROP INDEX IF EXISTS epi_694009_ndp11__mhc_allele;
DROP INDEX IF EXISTS epi_694009_ndp11__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_ndp11__product_lower;
DROP INDEX IF EXISTS epi_694009_ndp11__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_ndp11__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_ndp11__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_ndp11__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_ndp11__taxon_id;
DROP INDEX IF EXISTS epi_694009_ndp11__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_ndp11__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_ndp11__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_ndp11__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_ndp11__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_host_product;
DROP INDEX IF EXISTS epi_694009_ndp11__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT spike glycoprotein
-- 694009 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_spike_glycoprotein;

DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3a protein
-- 694009 can be replaced with the virus taxon id, while orf3a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf3a_protein;

DROP INDEX IF EXISTS epi_694009_orf3a_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf3a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF3b protein
-- 694009 can be replaced with the virus taxon id, while orf3b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf3b_protein;

DROP INDEX IF EXISTS epi_694009_orf3b_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf3b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT small envelope protein
-- 694009 can be replaced with the virus taxon id, while small_envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_small_envelope_protein;

DROP INDEX IF EXISTS epi_694009_small_envelope_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_small_envelope_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT membrane glycoprotein M
-- 694009 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_membrane_glycoprotein_m;

DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__cell_type;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__epi_an_start;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__host_tax_id;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__host_tax_name;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__iedb_id;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__is_linear;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__mhc_allele;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__product_lower;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__taxon_id;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_host_product;
DROP INDEX IF EXISTS epi_694009_membrane_glycoprotein_m__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF6 protein
-- 694009 can be replaced with the virus taxon id, while orf6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf6_protein;

DROP INDEX IF EXISTS epi_694009_orf6_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf6_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf6_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf6_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf6_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf6_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf6_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf6_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf6_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf6_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf6_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf6_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf6_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf6_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf6_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf6_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf6_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf6_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf6_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf6_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf6_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf6_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf6_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7a protein
-- 694009 can be replaced with the virus taxon id, while orf7a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf7a_protein;

DROP INDEX IF EXISTS epi_694009_orf7a_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf7a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF7b protein
-- 694009 can be replaced with the virus taxon id, while orf7b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf7b_protein;

DROP INDEX IF EXISTS epi_694009_orf7b_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf7b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8a protein
-- 694009 can be replaced with the virus taxon id, while orf8a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf8a_protein;

DROP INDEX IF EXISTS epi_694009_orf8a_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf8a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF8b protein
-- 694009 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf8b_protein;

DROP INDEX IF EXISTS epi_694009_orf8b_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf8b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT nucleocapsid protein
-- 694009 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_nucleocapsid_protein;

DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_nucleocapsid_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9b protein
-- 694009 can be replaced with the virus taxon id, while orf9b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf9b_protein;

DROP INDEX IF EXISTS epi_694009_orf9b_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf9b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sars_cov_1 and PROT ORF9a protein
-- 694009 can be replaced with the virus taxon id, while orf9a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_694009_orf9a_protein;

DROP INDEX IF EXISTS epi_694009_orf9a_protein__cell_type;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__iedb_id;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__is_linear;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__product_lower;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__taxon_id;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_694009_orf9a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT polyprotein
-- 11053 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_polyprotein;

DROP INDEX IF EXISTS epi_11053_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_11053_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_11053_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_11053_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_11053_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_11053_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_11053_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_11053_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_11053_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_11053_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_11053_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT anchored capsid protein ancC
-- 11053 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_anchored_capsid_protein_ancc;

DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__cell_type;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__epi_an_start;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__host_tax_id;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__host_tax_name;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__iedb_id;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__is_linear;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__mhc_allele;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__product_lower;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__taxon_id;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_host_product;
DROP INDEX IF EXISTS epi_11053_anchored_capsid_protein_ancc__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT capsid protein C
-- 11053 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_capsid_protein_c;

DROP INDEX IF EXISTS epi_11053_capsid_protein_c__cell_type;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__epi_an_start;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__host_tax_id;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__host_tax_name;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__iedb_id;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__is_linear;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__mhc_allele;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__product_lower;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__taxon_id;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_host_product;
DROP INDEX IF EXISTS epi_11053_capsid_protein_c__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT membrane glycoprotein precursor prM
-- 11053 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_membrane_glycoprotein_precur;

DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__cell_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__epi_an_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__host_tax_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__host_tax_name;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__iedb_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__is_linear;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__mhc_allele;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__product_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__taxon_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_host_product;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_precur__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT protein pr
-- 11053 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_protein_pr;

DROP INDEX IF EXISTS epi_11053_protein_pr__cell_type;
DROP INDEX IF EXISTS epi_11053_protein_pr__epi_an_start;
DROP INDEX IF EXISTS epi_11053_protein_pr__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_protein_pr__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_protein_pr__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_protein_pr__host_tax_id;
DROP INDEX IF EXISTS epi_11053_protein_pr__host_tax_name;
DROP INDEX IF EXISTS epi_11053_protein_pr__iedb_id;
DROP INDEX IF EXISTS epi_11053_protein_pr__is_linear;
DROP INDEX IF EXISTS epi_11053_protein_pr__mhc_allele;
DROP INDEX IF EXISTS epi_11053_protein_pr__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_protein_pr__product_lower;
DROP INDEX IF EXISTS epi_11053_protein_pr__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_protein_pr__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_protein_pr__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_protein_pr__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_protein_pr__taxon_id;
DROP INDEX IF EXISTS epi_11053_protein_pr__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_protein_pr__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_protein_pr__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_protein_pr__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_protein_pr__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_host_product;
DROP INDEX IF EXISTS epi_11053_protein_pr__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT membrane glycoprotein M
-- 11053 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_membrane_glycoprotein_m;

DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__cell_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__epi_an_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__host_tax_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__host_tax_name;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__iedb_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__is_linear;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__mhc_allele;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__product_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__taxon_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_host_product;
DROP INDEX IF EXISTS epi_11053_membrane_glycoprotein_m__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT envelope protein E
-- 11053 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_envelope_protein_e;

DROP INDEX IF EXISTS epi_11053_envelope_protein_e__cell_type;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__epi_an_start;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__host_tax_id;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__host_tax_name;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__iedb_id;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__is_linear;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__mhc_allele;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__product_lower;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__taxon_id;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_host_product;
DROP INDEX IF EXISTS epi_11053_envelope_protein_e__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS1
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns1;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns1__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS2A
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns2a;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS2B
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns2b;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns2b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS3
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns3;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS4A
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns4a;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT protein 2K
-- 11053 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_protein_2k;

DROP INDEX IF EXISTS epi_11053_protein_2k__cell_type;
DROP INDEX IF EXISTS epi_11053_protein_2k__epi_an_start;
DROP INDEX IF EXISTS epi_11053_protein_2k__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_protein_2k__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_protein_2k__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_protein_2k__host_tax_id;
DROP INDEX IF EXISTS epi_11053_protein_2k__host_tax_name;
DROP INDEX IF EXISTS epi_11053_protein_2k__iedb_id;
DROP INDEX IF EXISTS epi_11053_protein_2k__is_linear;
DROP INDEX IF EXISTS epi_11053_protein_2k__mhc_allele;
DROP INDEX IF EXISTS epi_11053_protein_2k__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_protein_2k__product_lower;
DROP INDEX IF EXISTS epi_11053_protein_2k__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_protein_2k__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_protein_2k__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_protein_2k__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_protein_2k__taxon_id;
DROP INDEX IF EXISTS epi_11053_protein_2k__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_protein_2k__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_protein_2k__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_protein_2k__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_protein_2k__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_host_product;
DROP INDEX IF EXISTS epi_11053_protein_2k__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT nonstructural protein NS4B
-- 11053 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_nonstructural_protein_ns4b;

DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__epi_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__host_tax_name;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__iedb_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__product_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__taxon_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_host_product;
DROP INDEX IF EXISTS epi_11053_nonstructural_protein_ns4b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_1 and PROT RNA-dependent RNA polymerase NS5
-- 11053 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11053_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_11053_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT polyprotein
-- 11060 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_polyprotein;

DROP INDEX IF EXISTS epi_11060_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_11060_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_11060_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_11060_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_11060_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_11060_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_11060_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_11060_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_11060_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_11060_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_11060_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT anchored capsid protein ancC
-- 11060 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_anchored_capsid_protein_ancc;

DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__cell_type;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__epi_an_start;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__host_tax_id;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__host_tax_name;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__iedb_id;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__is_linear;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__mhc_allele;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__product_lower;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__taxon_id;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_host_product;
DROP INDEX IF EXISTS epi_11060_anchored_capsid_protein_ancc__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT capsid protein C
-- 11060 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_capsid_protein_c;

DROP INDEX IF EXISTS epi_11060_capsid_protein_c__cell_type;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__epi_an_start;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__host_tax_id;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__host_tax_name;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__iedb_id;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__is_linear;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__mhc_allele;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__product_lower;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__taxon_id;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_host_product;
DROP INDEX IF EXISTS epi_11060_capsid_protein_c__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT membrane glycoprotein precursor prM
-- 11060 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_membrane_glycoprotein_precur;

DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__cell_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__epi_an_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__host_tax_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__host_tax_name;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__iedb_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__is_linear;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__mhc_allele;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__product_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__taxon_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_host_product;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_precur__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT protein pr
-- 11060 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_protein_pr;

DROP INDEX IF EXISTS epi_11060_protein_pr__cell_type;
DROP INDEX IF EXISTS epi_11060_protein_pr__epi_an_start;
DROP INDEX IF EXISTS epi_11060_protein_pr__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_protein_pr__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_protein_pr__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_protein_pr__host_tax_id;
DROP INDEX IF EXISTS epi_11060_protein_pr__host_tax_name;
DROP INDEX IF EXISTS epi_11060_protein_pr__iedb_id;
DROP INDEX IF EXISTS epi_11060_protein_pr__is_linear;
DROP INDEX IF EXISTS epi_11060_protein_pr__mhc_allele;
DROP INDEX IF EXISTS epi_11060_protein_pr__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_protein_pr__product_lower;
DROP INDEX IF EXISTS epi_11060_protein_pr__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_protein_pr__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_protein_pr__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_protein_pr__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_protein_pr__taxon_id;
DROP INDEX IF EXISTS epi_11060_protein_pr__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_protein_pr__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_protein_pr__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_protein_pr__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_protein_pr__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_host_product;
DROP INDEX IF EXISTS epi_11060_protein_pr__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT membrane glycoprotein M
-- 11060 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_membrane_glycoprotein_m;

DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__cell_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__epi_an_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__host_tax_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__host_tax_name;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__iedb_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__is_linear;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__mhc_allele;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__product_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__taxon_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_host_product;
DROP INDEX IF EXISTS epi_11060_membrane_glycoprotein_m__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT envelope protein E
-- 11060 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_envelope_protein_e;

DROP INDEX IF EXISTS epi_11060_envelope_protein_e__cell_type;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__epi_an_start;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__host_tax_id;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__host_tax_name;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__iedb_id;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__is_linear;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__mhc_allele;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__product_lower;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__taxon_id;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_host_product;
DROP INDEX IF EXISTS epi_11060_envelope_protein_e__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS1
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns1;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns1__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS2A
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns2a;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS2B
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns2b;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns2b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS3
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns3;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS4A
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns4a;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT protein 2K
-- 11060 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_protein_2k;

DROP INDEX IF EXISTS epi_11060_protein_2k__cell_type;
DROP INDEX IF EXISTS epi_11060_protein_2k__epi_an_start;
DROP INDEX IF EXISTS epi_11060_protein_2k__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_protein_2k__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_protein_2k__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_protein_2k__host_tax_id;
DROP INDEX IF EXISTS epi_11060_protein_2k__host_tax_name;
DROP INDEX IF EXISTS epi_11060_protein_2k__iedb_id;
DROP INDEX IF EXISTS epi_11060_protein_2k__is_linear;
DROP INDEX IF EXISTS epi_11060_protein_2k__mhc_allele;
DROP INDEX IF EXISTS epi_11060_protein_2k__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_protein_2k__product_lower;
DROP INDEX IF EXISTS epi_11060_protein_2k__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_protein_2k__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_protein_2k__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_protein_2k__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_protein_2k__taxon_id;
DROP INDEX IF EXISTS epi_11060_protein_2k__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_protein_2k__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_protein_2k__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_protein_2k__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_protein_2k__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_host_product;
DROP INDEX IF EXISTS epi_11060_protein_2k__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT nonstructural protein NS4B
-- 11060 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_nonstructural_protein_ns4b;

DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__epi_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__host_tax_name;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__iedb_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__product_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__taxon_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_host_product;
DROP INDEX IF EXISTS epi_11060_nonstructural_protein_ns4b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_2 and PROT RNA-dependent RNA polymerase NS5
-- 11060 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11060_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_11060_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT polyprotein
-- 11069 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_polyprotein;

DROP INDEX IF EXISTS epi_11069_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_11069_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_11069_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_11069_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_11069_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_11069_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_11069_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_11069_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_11069_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_11069_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_11069_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT anchored capsid protein ancC
-- 11069 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_anchored_capsid_protein_ancc;

DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__cell_type;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__epi_an_start;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__host_tax_id;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__host_tax_name;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__iedb_id;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__is_linear;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__mhc_allele;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__product_lower;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__taxon_id;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_host_product;
DROP INDEX IF EXISTS epi_11069_anchored_capsid_protein_ancc__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT capsid protein C
-- 11069 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_capsid_protein_c;

DROP INDEX IF EXISTS epi_11069_capsid_protein_c__cell_type;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__epi_an_start;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__host_tax_id;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__host_tax_name;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__iedb_id;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__is_linear;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__mhc_allele;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__product_lower;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__taxon_id;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_host_product;
DROP INDEX IF EXISTS epi_11069_capsid_protein_c__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT membrane glycoprotein precursor prM
-- 11069 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_membrane_glycoprotein_precur;

DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__cell_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__epi_an_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__host_tax_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__host_tax_name;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__iedb_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__is_linear;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__mhc_allele;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__product_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__taxon_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_host_product;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_precur__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT protein pr
-- 11069 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_protein_pr;

DROP INDEX IF EXISTS epi_11069_protein_pr__cell_type;
DROP INDEX IF EXISTS epi_11069_protein_pr__epi_an_start;
DROP INDEX IF EXISTS epi_11069_protein_pr__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_protein_pr__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_protein_pr__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_protein_pr__host_tax_id;
DROP INDEX IF EXISTS epi_11069_protein_pr__host_tax_name;
DROP INDEX IF EXISTS epi_11069_protein_pr__iedb_id;
DROP INDEX IF EXISTS epi_11069_protein_pr__is_linear;
DROP INDEX IF EXISTS epi_11069_protein_pr__mhc_allele;
DROP INDEX IF EXISTS epi_11069_protein_pr__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_protein_pr__product_lower;
DROP INDEX IF EXISTS epi_11069_protein_pr__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_protein_pr__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_protein_pr__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_protein_pr__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_protein_pr__taxon_id;
DROP INDEX IF EXISTS epi_11069_protein_pr__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_protein_pr__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_protein_pr__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_protein_pr__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_protein_pr__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_host_product;
DROP INDEX IF EXISTS epi_11069_protein_pr__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT membrane glycoprotein M
-- 11069 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_membrane_glycoprotein_m;

DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__cell_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__epi_an_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__host_tax_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__host_tax_name;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__iedb_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__is_linear;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__mhc_allele;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__product_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__taxon_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_host_product;
DROP INDEX IF EXISTS epi_11069_membrane_glycoprotein_m__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT envelope protein E
-- 11069 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_envelope_protein_e;

DROP INDEX IF EXISTS epi_11069_envelope_protein_e__cell_type;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__epi_an_start;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__host_tax_id;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__host_tax_name;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__iedb_id;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__is_linear;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__mhc_allele;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__product_lower;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__taxon_id;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_host_product;
DROP INDEX IF EXISTS epi_11069_envelope_protein_e__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS1
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns1;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns1__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS2A
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns2a;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS2B
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns2b;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns2b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS3
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns3;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS4A
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns4a;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT protein 2K
-- 11069 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_protein_2k;

DROP INDEX IF EXISTS epi_11069_protein_2k__cell_type;
DROP INDEX IF EXISTS epi_11069_protein_2k__epi_an_start;
DROP INDEX IF EXISTS epi_11069_protein_2k__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_protein_2k__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_protein_2k__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_protein_2k__host_tax_id;
DROP INDEX IF EXISTS epi_11069_protein_2k__host_tax_name;
DROP INDEX IF EXISTS epi_11069_protein_2k__iedb_id;
DROP INDEX IF EXISTS epi_11069_protein_2k__is_linear;
DROP INDEX IF EXISTS epi_11069_protein_2k__mhc_allele;
DROP INDEX IF EXISTS epi_11069_protein_2k__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_protein_2k__product_lower;
DROP INDEX IF EXISTS epi_11069_protein_2k__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_protein_2k__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_protein_2k__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_protein_2k__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_protein_2k__taxon_id;
DROP INDEX IF EXISTS epi_11069_protein_2k__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_protein_2k__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_protein_2k__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_protein_2k__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_protein_2k__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_host_product;
DROP INDEX IF EXISTS epi_11069_protein_2k__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT nonstructural protein NS4B
-- 11069 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_nonstructural_protein_ns4b;

DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__epi_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__host_tax_name;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__iedb_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__product_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__taxon_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_host_product;
DROP INDEX IF EXISTS epi_11069_nonstructural_protein_ns4b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_3 and PROT RNA-dependent RNA polymerase NS5
-- 11069 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11069_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_11069_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT polyprotein
-- 11070 can be replaced with the virus taxon id, while polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_polyprotein;

DROP INDEX IF EXISTS epi_11070_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_11070_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_11070_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_11070_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_11070_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_11070_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_11070_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_11070_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_11070_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_11070_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_11070_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT anchored capsid protein ancC
-- 11070 can be replaced with the virus taxon id, while anchored_capsid_protein_ancc can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_anchored_capsid_protein_ancc;

DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__cell_type;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__epi_an_start;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__host_tax_id;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__host_tax_name;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__iedb_id;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__is_linear;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__mhc_allele;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__product_lower;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__taxon_id;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_host_product;
DROP INDEX IF EXISTS epi_11070_anchored_capsid_protein_ancc__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT capsid protein C
-- 11070 can be replaced with the virus taxon id, while capsid_protein_c can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_capsid_protein_c;

DROP INDEX IF EXISTS epi_11070_capsid_protein_c__cell_type;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__epi_an_start;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__host_tax_id;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__host_tax_name;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__iedb_id;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__is_linear;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__mhc_allele;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__product_lower;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__taxon_id;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_host_product;
DROP INDEX IF EXISTS epi_11070_capsid_protein_c__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein precursor prM
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_precur can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_membrane_glycoprotein_precur;

DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__cell_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__epi_an_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__host_tax_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__host_tax_name;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__iedb_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__is_linear;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__mhc_allele;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__product_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__taxon_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_host_product;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_precur__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT protein pr
-- 11070 can be replaced with the virus taxon id, while protein_pr can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_protein_pr;

DROP INDEX IF EXISTS epi_11070_protein_pr__cell_type;
DROP INDEX IF EXISTS epi_11070_protein_pr__epi_an_start;
DROP INDEX IF EXISTS epi_11070_protein_pr__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_protein_pr__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_protein_pr__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_protein_pr__host_tax_id;
DROP INDEX IF EXISTS epi_11070_protein_pr__host_tax_name;
DROP INDEX IF EXISTS epi_11070_protein_pr__iedb_id;
DROP INDEX IF EXISTS epi_11070_protein_pr__is_linear;
DROP INDEX IF EXISTS epi_11070_protein_pr__mhc_allele;
DROP INDEX IF EXISTS epi_11070_protein_pr__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_protein_pr__product_lower;
DROP INDEX IF EXISTS epi_11070_protein_pr__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_protein_pr__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_protein_pr__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_protein_pr__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_protein_pr__taxon_id;
DROP INDEX IF EXISTS epi_11070_protein_pr__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_protein_pr__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_protein_pr__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_protein_pr__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_protein_pr__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_host_product;
DROP INDEX IF EXISTS epi_11070_protein_pr__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT membrane glycoprotein M
-- 11070 can be replaced with the virus taxon id, while membrane_glycoprotein_m can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_membrane_glycoprotein_m;

DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__cell_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__epi_an_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__host_tax_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__host_tax_name;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__iedb_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__is_linear;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__mhc_allele;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__product_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__taxon_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_host_product;
DROP INDEX IF EXISTS epi_11070_membrane_glycoprotein_m__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT envelope protein E
-- 11070 can be replaced with the virus taxon id, while envelope_protein_e can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_envelope_protein_e;

DROP INDEX IF EXISTS epi_11070_envelope_protein_e__cell_type;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__epi_an_start;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__host_tax_id;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__host_tax_name;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__iedb_id;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__is_linear;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__mhc_allele;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__product_lower;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__taxon_id;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_host_product;
DROP INDEX IF EXISTS epi_11070_envelope_protein_e__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS1
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns1 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns1;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns1__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2A
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns2a;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS2B
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns2b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns2b;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns2b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS3
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns3 can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns3;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns3__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4A
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4a can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns4a;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4a__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT protein 2K
-- 11070 can be replaced with the virus taxon id, while protein_2k can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_protein_2k;

DROP INDEX IF EXISTS epi_11070_protein_2k__cell_type;
DROP INDEX IF EXISTS epi_11070_protein_2k__epi_an_start;
DROP INDEX IF EXISTS epi_11070_protein_2k__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_protein_2k__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_protein_2k__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_protein_2k__host_tax_id;
DROP INDEX IF EXISTS epi_11070_protein_2k__host_tax_name;
DROP INDEX IF EXISTS epi_11070_protein_2k__iedb_id;
DROP INDEX IF EXISTS epi_11070_protein_2k__is_linear;
DROP INDEX IF EXISTS epi_11070_protein_2k__mhc_allele;
DROP INDEX IF EXISTS epi_11070_protein_2k__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_protein_2k__product_lower;
DROP INDEX IF EXISTS epi_11070_protein_2k__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_protein_2k__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_protein_2k__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_protein_2k__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_protein_2k__taxon_id;
DROP INDEX IF EXISTS epi_11070_protein_2k__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_protein_2k__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_protein_2k__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_protein_2k__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_protein_2k__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_host_product;
DROP INDEX IF EXISTS epi_11070_protein_2k__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT nonstructural protein NS4B
-- 11070 can be replaced with the virus taxon id, while nonstructural_protein_ns4b can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_nonstructural_protein_ns4b;

DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__epi_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__host_tax_name;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__iedb_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__product_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__taxon_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_host_product;
DROP INDEX IF EXISTS epi_11070_nonstructural_protein_ns4b__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR dengue_4 and PROT RNA-dependent RNA polymerase NS5
-- 11070 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_11070_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_11070_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT 1AB polyprotein
-- 1335626 can be replaced with the virus taxon id, while 1ab_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_1ab_polyprotein;

DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_1ab_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT RNA-dependent RNA polymerase
-- 1335626 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT Hel
-- 1335626 can be replaced with the virus taxon id, while hel can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_hel;

DROP INDEX IF EXISTS epi_1335626_hel__cell_type;
DROP INDEX IF EXISTS epi_1335626_hel__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_hel__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_hel__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_hel__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_hel__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_hel__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_hel__iedb_id;
DROP INDEX IF EXISTS epi_1335626_hel__is_linear;
DROP INDEX IF EXISTS epi_1335626_hel__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_hel__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_hel__product_lower;
DROP INDEX IF EXISTS epi_1335626_hel__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_hel__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_hel__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_hel__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_hel__taxon_id;
DROP INDEX IF EXISTS epi_1335626_hel__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_hel__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_hel__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_hel__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_hel__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_hel__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_hel__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_hel__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_hel__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_hel__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_hel__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT ExoN
-- 1335626 can be replaced with the virus taxon id, while exon can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_exon;

DROP INDEX IF EXISTS epi_1335626_exon__cell_type;
DROP INDEX IF EXISTS epi_1335626_exon__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_exon__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_exon__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_exon__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_exon__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_exon__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_exon__iedb_id;
DROP INDEX IF EXISTS epi_1335626_exon__is_linear;
DROP INDEX IF EXISTS epi_1335626_exon__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_exon__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_exon__product_lower;
DROP INDEX IF EXISTS epi_1335626_exon__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_exon__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_exon__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_exon__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_exon__taxon_id;
DROP INDEX IF EXISTS epi_1335626_exon__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_exon__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_exon__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_exon__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_exon__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_exon__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_exon__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_exon__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_exon__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_exon__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_exon__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT NendoU
-- 1335626 can be replaced with the virus taxon id, while nendou can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nendou;

DROP INDEX IF EXISTS epi_1335626_nendou__cell_type;
DROP INDEX IF EXISTS epi_1335626_nendou__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nendou__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nendou__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nendou__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nendou__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nendou__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nendou__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nendou__is_linear;
DROP INDEX IF EXISTS epi_1335626_nendou__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nendou__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nendou__product_lower;
DROP INDEX IF EXISTS epi_1335626_nendou__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nendou__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nendou__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nendou__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nendou__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nendou__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nendou__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nendou__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nendou__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nendou__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nendou__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT 2'-O-methyltransferase
-- 1335626 can be replaced with the virus taxon id, while 2_o_methyltransferase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_2_o_methyltransferase;

DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__cell_type;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__iedb_id;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__is_linear;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__product_lower;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__taxon_id;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_2_o_methyltransferase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT 1A polyprotein
-- 1335626 can be replaced with the virus taxon id, while 1a_polyprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_1a_polyprotein;

DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__cell_type;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__is_linear;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__product_lower;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_1a_polyprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp1 protein
-- 1335626 can be replaced with the virus taxon id, while nsp1_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp1_protein;

DROP INDEX IF EXISTS epi_1335626_nsp1_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp1_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp2 protein
-- 1335626 can be replaced with the virus taxon id, while nsp2_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp2_protein;

DROP INDEX IF EXISTS epi_1335626_nsp2_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp2_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp3 protein
-- 1335626 can be replaced with the virus taxon id, while nsp3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp3_protein;

DROP INDEX IF EXISTS epi_1335626_nsp3_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp3_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp4 protein
-- 1335626 can be replaced with the virus taxon id, while nsp4_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp4_protein;

DROP INDEX IF EXISTS epi_1335626_nsp4_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp4_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp5 protein
-- 1335626 can be replaced with the virus taxon id, while nsp5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp5_protein;

DROP INDEX IF EXISTS epi_1335626_nsp5_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp5_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp6 protein
-- 1335626 can be replaced with the virus taxon id, while nsp6_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp6_protein;

DROP INDEX IF EXISTS epi_1335626_nsp6_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp6_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp7 protein
-- 1335626 can be replaced with the virus taxon id, while nsp7_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp7_protein;

DROP INDEX IF EXISTS epi_1335626_nsp7_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp7_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp8 protein
-- 1335626 can be replaced with the virus taxon id, while nsp8_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp8_protein;

DROP INDEX IF EXISTS epi_1335626_nsp8_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp8_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp9 protein
-- 1335626 can be replaced with the virus taxon id, while nsp9_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp9_protein;

DROP INDEX IF EXISTS epi_1335626_nsp9_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp9_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp10 protein
-- 1335626 can be replaced with the virus taxon id, while nsp10_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp10_protein;

DROP INDEX IF EXISTS epi_1335626_nsp10_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp10_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nsp11 protein
-- 1335626 can be replaced with the virus taxon id, while nsp11_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nsp11_protein;

DROP INDEX IF EXISTS epi_1335626_nsp11_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nsp11_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT spike protein
-- 1335626 can be replaced with the virus taxon id, while spike_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_spike_protein;

DROP INDEX IF EXISTS epi_1335626_spike_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_spike_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_spike_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_spike_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_spike_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_spike_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_spike_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_spike_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_spike_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_spike_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_spike_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_spike_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_spike_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_spike_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_spike_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_spike_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_spike_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_spike_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_spike_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_spike_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_spike_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_spike_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_spike_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT NS3 protein
-- 1335626 can be replaced with the virus taxon id, while ns3_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_ns3_protein;

DROP INDEX IF EXISTS epi_1335626_ns3_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_ns3_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT NS4A protein
-- 1335626 can be replaced with the virus taxon id, while ns4a_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_ns4a_protein;

DROP INDEX IF EXISTS epi_1335626_ns4a_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_ns4a_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT NS4B protein
-- 1335626 can be replaced with the virus taxon id, while ns4b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_ns4b_protein;

DROP INDEX IF EXISTS epi_1335626_ns4b_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_ns4b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT NS5 protein
-- 1335626 can be replaced with the virus taxon id, while ns5_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_ns5_protein;

DROP INDEX IF EXISTS epi_1335626_ns5_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_ns5_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT envelope protein
-- 1335626 can be replaced with the virus taxon id, while envelope_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_envelope_protein;

DROP INDEX IF EXISTS epi_1335626_envelope_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_envelope_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT membrane protein
-- 1335626 can be replaced with the virus taxon id, while membrane_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_membrane_protein;

DROP INDEX IF EXISTS epi_1335626_membrane_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_membrane_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT nucleocapsid protein
-- 1335626 can be replaced with the virus taxon id, while nucleocapsid_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_nucleocapsid_protein;

DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_nucleocapsid_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR mers and PROT ORF8b protein
-- 1335626 can be replaced with the virus taxon id, while orf8b_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_1335626_orf8b_protein;

DROP INDEX IF EXISTS epi_1335626_orf8b_protein__cell_type;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__epi_an_start;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__host_tax_id;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__host_tax_name;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__iedb_id;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__is_linear;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__mhc_allele;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__product_lower;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__taxon_id;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_host_product;
DROP INDEX IF EXISTS epi_1335626_orf8b_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT nucleoprotein
-- 186538 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_nucleoprotein;

DROP INDEX IF EXISTS epi_186538_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT polymerase complex protein
-- 186538 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT matrix protein
-- 186538 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_matrix_protein;

DROP INDEX IF EXISTS epi_186538_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_186538_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_186538_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_186538_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_186538_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_186538_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT spike glycoprotein
-- 186538 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_spike_glycoprotein;

DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT small secreted glycoprotein
-- 186538 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT second secreted glycoprotein
-- 186538 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT minor nucleoprotein
-- 186538 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT membrane-associated protein
-- 186538 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_membrane_associated_protein;

DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186538_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR zaire_ebolavirus and PROT RNA-dependent RNA polymerase
-- 186538 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186538_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_186538_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT nucleoprotein
-- 186540 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_nucleoprotein;

DROP INDEX IF EXISTS epi_186540_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT polymerase complex protein
-- 186540 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT matrix protein
-- 186540 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_matrix_protein;

DROP INDEX IF EXISTS epi_186540_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_186540_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_186540_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_186540_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_186540_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_186540_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT spike glycoprotein
-- 186540 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_spike_glycoprotein;

DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT small secreted glycoprotein
-- 186540 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT second secreted glycoprotein
-- 186540 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT minor nucleoprotein
-- 186540 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT membrane-associated protein
-- 186540 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_membrane_associated_protein;

DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186540_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR sudan_ebolavirus and PROT RNA-dependent RNA polymerase
-- 186540 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186540_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_186540_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT nucleoprotein
-- 186539 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_nucleoprotein;

DROP INDEX IF EXISTS epi_186539_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT polymerase complex protein
-- 186539 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT matrix protein
-- 186539 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_matrix_protein;

DROP INDEX IF EXISTS epi_186539_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_186539_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_186539_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_186539_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_186539_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_186539_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT spike glycoprotein
-- 186539 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_spike_glycoprotein;

DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT small secreted glycoprotein
-- 186539 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT second secreted glycoprotein
-- 186539 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT minor nucleoprotein
-- 186539 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT membrane-associated protein
-- 186539 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_membrane_associated_protein;

DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186539_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR reston_ebolavirus and PROT RNA-dependent RNA polymerase
-- 186539 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186539_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_186539_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT nucleoprotein
-- 565995 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_nucleoprotein;

DROP INDEX IF EXISTS epi_565995_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT polymerase complex protein
-- 565995 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT matrix protein
-- 565995 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_matrix_protein;

DROP INDEX IF EXISTS epi_565995_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_565995_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_565995_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_565995_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_565995_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_565995_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT spike glycoprotein
-- 565995 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_spike_glycoprotein;

DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT small secreted glycoprotein
-- 565995 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT second secreted glycoprotein
-- 565995 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT minor nucleoprotein
-- 565995 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT membrane-associated protein
-- 565995 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_membrane_associated_protein;

DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_565995_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bundibugyo_ebolavirus and PROT RNA-dependent RNA polymerase
-- 565995 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_565995_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_565995_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT nucleoprotein
-- 2010960 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_nucleoprotein;

DROP INDEX IF EXISTS epi_2010960_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT polymerase complex protein
-- 2010960 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT matrix protein
-- 2010960 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_matrix_protein;

DROP INDEX IF EXISTS epi_2010960_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT spike glycoprotein
-- 2010960 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_spike_glycoprotein;

DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT small secreted glycoprotein
-- 2010960 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT second secreted glycoprotein
-- 2010960 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT minor nucleoprotein
-- 2010960 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT membrane-associated protein
-- 2010960 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_membrane_associated_protein;

DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR bombali_ebolavirus and PROT RNA-dependent RNA polymerase
-- 2010960 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_2010960_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_2010960_rna_dependent_rna_polymerase__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT nucleoprotein
-- 186541 can be replaced with the virus taxon id, while nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_nucleoprotein;

DROP INDEX IF EXISTS epi_186541_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT polymerase complex protein
-- 186541 can be replaced with the virus taxon id, while polymerase_complex_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_polymerase_complex_protein;

DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__cell_type;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__iedb_id;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__is_linear;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__product_lower;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__taxon_id;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_polymerase_complex_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT matrix protein
-- 186541 can be replaced with the virus taxon id, while matrix_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_matrix_protein;

DROP INDEX IF EXISTS epi_186541_matrix_protein__cell_type;
DROP INDEX IF EXISTS epi_186541_matrix_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_matrix_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_matrix_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_matrix_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_matrix_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_matrix_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_matrix_protein__iedb_id;
DROP INDEX IF EXISTS epi_186541_matrix_protein__is_linear;
DROP INDEX IF EXISTS epi_186541_matrix_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_matrix_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_matrix_protein__product_lower;
DROP INDEX IF EXISTS epi_186541_matrix_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_matrix_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_matrix_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_matrix_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_matrix_protein__taxon_id;
DROP INDEX IF EXISTS epi_186541_matrix_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_matrix_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_matrix_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_matrix_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_matrix_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_matrix_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT spike glycoprotein
-- 186541 can be replaced with the virus taxon id, while spike_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_spike_glycoprotein;

DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_spike_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT small secreted glycoprotein
-- 186541 can be replaced with the virus taxon id, while small_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_small_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_small_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT second secreted glycoprotein
-- 186541 can be replaced with the virus taxon id, while second_secreted_glycoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_second_secreted_glycoprotein;

DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__cell_type;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__is_linear;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__product_lower;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_second_secreted_glycoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT minor nucleoprotein
-- 186541 can be replaced with the virus taxon id, while minor_nucleoprotein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_minor_nucleoprotein;

DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__cell_type;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__iedb_id;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__is_linear;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__product_lower;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__taxon_id;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_minor_nucleoprotein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT membrane-associated protein
-- 186541 can be replaced with the virus taxon id, while membrane_associated_protein can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_membrane_associated_protein;

DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__cell_type;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__epi_an_start;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__host_tax_id;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__host_tax_name;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__iedb_id;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__is_linear;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__mhc_allele;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__product_lower;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__taxon_id;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_host_product;
DROP INDEX IF EXISTS epi_186541_membrane_associated_protein__vir_host_resp_freq;


-- DROP TABLES 'N INDEXES OF VIR tai_forest_ebolavirus and PROT RNA-dependent RNA polymerase
-- 186541 can be replaced with the virus taxon id, while rna_dependent_rna_polymerase can be long 28 chars max to comply
-- with postgres limit on DB object names (max 63 chars allowed) on views, tables, constraints and indexes.
DROP TABLE IF EXISTS public.epitope_186541_rna_dependent_rna_polymerase;

DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__cell_type;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__epi_an_start;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__epi_an_nstop;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__epi_frag_an_start;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__epi_frag_an_stop;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__host_tax_id;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__host_tax_name;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__iedb_id;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__is_linear;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__mhc_allele;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__mhc_class_lower;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__product_lower;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__response_freq_pos;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__seq_aa_alt;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__seq_aa_orig;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__start_aa_orig;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__taxon_id;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__taxon_name_lower;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__variant_aa_length;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__variant_aa_type;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_n_host_tax_id;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_host_cell_type;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_host_epi_start;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__virus_host_epi_stop;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__virus_host_is_linear;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_host_mhc_allele;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_host_product;
DROP INDEX IF EXISTS epi_186541_rna_dependent_rna_polymerase__vir_host_resp_freq;


