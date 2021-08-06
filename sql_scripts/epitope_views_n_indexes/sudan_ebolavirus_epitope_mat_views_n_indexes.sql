-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN nucleoprotein
CREATE MATERIALIZED VIEW public.epitope_186540_nucleoprote
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
             AND ann.product = 'nucleoprotein'
             AND epi.protein_name = 'nucleoprotein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_nucleoprote
    OWNER TO geco;

CREATE INDEX epi_186540_nucleoprote__cell_type__idx
    ON public.epitope_186540_nucleoprote USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__epi_annotation_start__idx
    ON public.epitope_186540_nucleoprote USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__epi_annotation_stop__idx
    ON public.epitope_186540_nucleoprote USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__epi_frag_annotation_start__i
    ON public.epitope_186540_nucleoprote USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__epi_frag_annotation_stop__id
    ON public.epitope_186540_nucleoprote USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__host_taxon_id__idx
    ON public.epitope_186540_nucleoprote USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__host_taxon_name_lower__idx
    ON public.epitope_186540_nucleoprote USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__iedb_epitope_id__idx
    ON public.epitope_186540_nucleoprote USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__is_linear__idx
    ON public.epitope_186540_nucleoprote USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__mhc_allele__idx
    ON public.epitope_186540_nucleoprote USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__mhc_class_lower__idx
    ON public.epitope_186540_nucleoprote USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__product_lower__idx
    ON public.epitope_186540_nucleoprote USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__response_frequency_pos__idx
    ON public.epitope_186540_nucleoprote USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__sequence_aa_alternative__idx
    ON public.epitope_186540_nucleoprote USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__sequence_aa_original__idx
    ON public.epitope_186540_nucleoprote USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__start_aa_original__idx
    ON public.epitope_186540_nucleoprote USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__taxon_id__idx
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__taxon_name_lower__idx
    ON public.epitope_186540_nucleoprote USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__variant_aa_length__idx
    ON public.epitope_186540_nucleoprote USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__variant_aa_type__idx
    ON public.epitope_186540_nucleoprote USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_cell_type__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_epi_start__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_epi_stop__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_is_linear__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_mhc_allele__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_product__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_nucleoprote__virus_host_resp_freq__i
    ON public.epitope_186540_nucleoprote USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN polymerase complex protein
CREATE MATERIALIZED VIEW public.epitope_186540_polymerase_
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
             AND ann.product = 'polymerase complex protein'
             AND epi.protein_name = 'polymerase complex protein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_polymerase_
    OWNER TO geco;

CREATE INDEX epi_186540_polymerase___cell_type__idx
    ON public.epitope_186540_polymerase_ USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___epi_annotation_start__idx
    ON public.epitope_186540_polymerase_ USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___epi_annotation_stop__idx
    ON public.epitope_186540_polymerase_ USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___epi_frag_annotation_start__i
    ON public.epitope_186540_polymerase_ USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___epi_frag_annotation_stop__id
    ON public.epitope_186540_polymerase_ USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___host_taxon_id__idx
    ON public.epitope_186540_polymerase_ USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___host_taxon_name_lower__idx
    ON public.epitope_186540_polymerase_ USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___iedb_epitope_id__idx
    ON public.epitope_186540_polymerase_ USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___is_linear__idx
    ON public.epitope_186540_polymerase_ USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___mhc_allele__idx
    ON public.epitope_186540_polymerase_ USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___mhc_class_lower__idx
    ON public.epitope_186540_polymerase_ USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___product_lower__idx
    ON public.epitope_186540_polymerase_ USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___response_frequency_pos__idx
    ON public.epitope_186540_polymerase_ USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___sequence_aa_alternative__idx
    ON public.epitope_186540_polymerase_ USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___sequence_aa_original__idx
    ON public.epitope_186540_polymerase_ USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___start_aa_original__idx
    ON public.epitope_186540_polymerase_ USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___taxon_id__idx
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___taxon_name_lower__idx
    ON public.epitope_186540_polymerase_ USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___variant_aa_length__idx
    ON public.epitope_186540_polymerase_ USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___variant_aa_type__idx
    ON public.epitope_186540_polymerase_ USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_cell_type__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_epi_start__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_epi_stop__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_is_linear__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_mhc_allele__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_product__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_polymerase___virus_host_resp_freq__i
    ON public.epitope_186540_polymerase_ USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN matrix protein
CREATE MATERIALIZED VIEW public.epitope_186540_matrix_prot
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
             AND ann.product = 'matrix protein'
             AND epi.protein_name = 'matrix protein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_matrix_prot
    OWNER TO geco;

CREATE INDEX epi_186540_matrix_prot__cell_type__idx
    ON public.epitope_186540_matrix_prot USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__epi_annotation_start__idx
    ON public.epitope_186540_matrix_prot USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__epi_annotation_stop__idx
    ON public.epitope_186540_matrix_prot USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__epi_frag_annotation_start__i
    ON public.epitope_186540_matrix_prot USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__epi_frag_annotation_stop__id
    ON public.epitope_186540_matrix_prot USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__host_taxon_id__idx
    ON public.epitope_186540_matrix_prot USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__host_taxon_name_lower__idx
    ON public.epitope_186540_matrix_prot USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__iedb_epitope_id__idx
    ON public.epitope_186540_matrix_prot USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__is_linear__idx
    ON public.epitope_186540_matrix_prot USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__mhc_allele__idx
    ON public.epitope_186540_matrix_prot USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__mhc_class_lower__idx
    ON public.epitope_186540_matrix_prot USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__product_lower__idx
    ON public.epitope_186540_matrix_prot USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__response_frequency_pos__idx
    ON public.epitope_186540_matrix_prot USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__sequence_aa_alternative__idx
    ON public.epitope_186540_matrix_prot USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__sequence_aa_original__idx
    ON public.epitope_186540_matrix_prot USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__start_aa_original__idx
    ON public.epitope_186540_matrix_prot USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__taxon_id__idx
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__taxon_name_lower__idx
    ON public.epitope_186540_matrix_prot USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__variant_aa_length__idx
    ON public.epitope_186540_matrix_prot USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__variant_aa_type__idx
    ON public.epitope_186540_matrix_prot USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_cell_type__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_epi_start__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_epi_stop__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_is_linear__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_mhc_allele__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_product__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_matrix_prot__virus_host_resp_freq__i
    ON public.epitope_186540_matrix_prot USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN spike glycoprotein
CREATE MATERIALIZED VIEW public.epitope_186540_spike_glyco
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
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_spike_glyco
    OWNER TO geco;

CREATE INDEX epi_186540_spike_glyco__cell_type__idx
    ON public.epitope_186540_spike_glyco USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__epi_annotation_start__idx
    ON public.epitope_186540_spike_glyco USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__epi_annotation_stop__idx
    ON public.epitope_186540_spike_glyco USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__epi_frag_annotation_start__i
    ON public.epitope_186540_spike_glyco USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__epi_frag_annotation_stop__id
    ON public.epitope_186540_spike_glyco USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__host_taxon_id__idx
    ON public.epitope_186540_spike_glyco USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__host_taxon_name_lower__idx
    ON public.epitope_186540_spike_glyco USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__iedb_epitope_id__idx
    ON public.epitope_186540_spike_glyco USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__is_linear__idx
    ON public.epitope_186540_spike_glyco USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__mhc_allele__idx
    ON public.epitope_186540_spike_glyco USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__mhc_class_lower__idx
    ON public.epitope_186540_spike_glyco USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__product_lower__idx
    ON public.epitope_186540_spike_glyco USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__response_frequency_pos__idx
    ON public.epitope_186540_spike_glyco USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__sequence_aa_alternative__idx
    ON public.epitope_186540_spike_glyco USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__sequence_aa_original__idx
    ON public.epitope_186540_spike_glyco USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__start_aa_original__idx
    ON public.epitope_186540_spike_glyco USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__taxon_id__idx
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__taxon_name_lower__idx
    ON public.epitope_186540_spike_glyco USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__variant_aa_length__idx
    ON public.epitope_186540_spike_glyco USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__variant_aa_type__idx
    ON public.epitope_186540_spike_glyco USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_cell_type__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_epi_start__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_epi_stop__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_is_linear__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_mhc_allele__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_product__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_spike_glyco__virus_host_resp_freq__i
    ON public.epitope_186540_spike_glyco USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN small secreted glycoprotein
CREATE MATERIALIZED VIEW public.epitope_186540_small_secre
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
             AND ann.product = 'small secreted glycoprotein'
             AND epi.protein_name = 'small secreted glycoprotein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_small_secre
    OWNER TO geco;

CREATE INDEX epi_186540_small_secre__cell_type__idx
    ON public.epitope_186540_small_secre USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__epi_annotation_start__idx
    ON public.epitope_186540_small_secre USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__epi_annotation_stop__idx
    ON public.epitope_186540_small_secre USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__epi_frag_annotation_start__i
    ON public.epitope_186540_small_secre USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__epi_frag_annotation_stop__id
    ON public.epitope_186540_small_secre USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__host_taxon_id__idx
    ON public.epitope_186540_small_secre USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__host_taxon_name_lower__idx
    ON public.epitope_186540_small_secre USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__iedb_epitope_id__idx
    ON public.epitope_186540_small_secre USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__is_linear__idx
    ON public.epitope_186540_small_secre USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__mhc_allele__idx
    ON public.epitope_186540_small_secre USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__mhc_class_lower__idx
    ON public.epitope_186540_small_secre USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__product_lower__idx
    ON public.epitope_186540_small_secre USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__response_frequency_pos__idx
    ON public.epitope_186540_small_secre USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__sequence_aa_alternative__idx
    ON public.epitope_186540_small_secre USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__sequence_aa_original__idx
    ON public.epitope_186540_small_secre USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__start_aa_original__idx
    ON public.epitope_186540_small_secre USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__taxon_id__idx
    ON public.epitope_186540_small_secre USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__taxon_name_lower__idx
    ON public.epitope_186540_small_secre USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__variant_aa_length__idx
    ON public.epitope_186540_small_secre USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__variant_aa_type__idx
    ON public.epitope_186540_small_secre USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_cell_type__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_epi_start__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_epi_stop__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_is_linear__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_mhc_allele__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_product__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_small_secre__virus_host_resp_freq__i
    ON public.epitope_186540_small_secre USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN second secreted glycoprotein
CREATE MATERIALIZED VIEW public.epitope_186540_second_secr
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
             AND ann.product = 'second secreted glycoprotein'
             AND epi.protein_name = 'second secreted glycoprotein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_second_secr
    OWNER TO geco;

CREATE INDEX epi_186540_second_secr__cell_type__idx
    ON public.epitope_186540_second_secr USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__epi_annotation_start__idx
    ON public.epitope_186540_second_secr USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__epi_annotation_stop__idx
    ON public.epitope_186540_second_secr USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__epi_frag_annotation_start__i
    ON public.epitope_186540_second_secr USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__epi_frag_annotation_stop__id
    ON public.epitope_186540_second_secr USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__host_taxon_id__idx
    ON public.epitope_186540_second_secr USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__host_taxon_name_lower__idx
    ON public.epitope_186540_second_secr USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__iedb_epitope_id__idx
    ON public.epitope_186540_second_secr USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__is_linear__idx
    ON public.epitope_186540_second_secr USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__mhc_allele__idx
    ON public.epitope_186540_second_secr USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__mhc_class_lower__idx
    ON public.epitope_186540_second_secr USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__product_lower__idx
    ON public.epitope_186540_second_secr USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__response_frequency_pos__idx
    ON public.epitope_186540_second_secr USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__sequence_aa_alternative__idx
    ON public.epitope_186540_second_secr USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__sequence_aa_original__idx
    ON public.epitope_186540_second_secr USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__start_aa_original__idx
    ON public.epitope_186540_second_secr USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__taxon_id__idx
    ON public.epitope_186540_second_secr USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__taxon_name_lower__idx
    ON public.epitope_186540_second_secr USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__variant_aa_length__idx
    ON public.epitope_186540_second_secr USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__variant_aa_type__idx
    ON public.epitope_186540_second_secr USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_cell_type__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_epi_start__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_epi_stop__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_is_linear__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_mhc_allele__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_product__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_second_secr__virus_host_resp_freq__i
    ON public.epitope_186540_second_secr USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN minor nucleoprotein
CREATE MATERIALIZED VIEW public.epitope_186540_minor_nucle
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
             AND ann.product = 'minor nucleoprotein'
             AND epi.protein_name = 'minor nucleoprotein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_minor_nucle
    OWNER TO geco;

CREATE INDEX epi_186540_minor_nucle__cell_type__idx
    ON public.epitope_186540_minor_nucle USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__epi_annotation_start__idx
    ON public.epitope_186540_minor_nucle USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__epi_annotation_stop__idx
    ON public.epitope_186540_minor_nucle USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__epi_frag_annotation_start__i
    ON public.epitope_186540_minor_nucle USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__epi_frag_annotation_stop__id
    ON public.epitope_186540_minor_nucle USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__host_taxon_id__idx
    ON public.epitope_186540_minor_nucle USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__host_taxon_name_lower__idx
    ON public.epitope_186540_minor_nucle USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__iedb_epitope_id__idx
    ON public.epitope_186540_minor_nucle USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__is_linear__idx
    ON public.epitope_186540_minor_nucle USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__mhc_allele__idx
    ON public.epitope_186540_minor_nucle USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__mhc_class_lower__idx
    ON public.epitope_186540_minor_nucle USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__product_lower__idx
    ON public.epitope_186540_minor_nucle USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__response_frequency_pos__idx
    ON public.epitope_186540_minor_nucle USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__sequence_aa_alternative__idx
    ON public.epitope_186540_minor_nucle USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__sequence_aa_original__idx
    ON public.epitope_186540_minor_nucle USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__start_aa_original__idx
    ON public.epitope_186540_minor_nucle USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__taxon_id__idx
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__taxon_name_lower__idx
    ON public.epitope_186540_minor_nucle USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__variant_aa_length__idx
    ON public.epitope_186540_minor_nucle USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__variant_aa_type__idx
    ON public.epitope_186540_minor_nucle USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_cell_type__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_epi_start__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_epi_stop__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_is_linear__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_mhc_allele__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_product__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_minor_nucle__virus_host_resp_freq__i
    ON public.epitope_186540_minor_nucle USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN membrane-associated protein
CREATE MATERIALIZED VIEW public.epitope_186540_membrane_as
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
             AND ann.product = 'membrane-associated protein'
             AND epi.protein_name = 'membrane-associated protein'
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_membrane_as
    OWNER TO geco;

CREATE INDEX epi_186540_membrane_as__cell_type__idx
    ON public.epitope_186540_membrane_as USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__epi_annotation_start__idx
    ON public.epitope_186540_membrane_as USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__epi_annotation_stop__idx
    ON public.epitope_186540_membrane_as USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__epi_frag_annotation_start__i
    ON public.epitope_186540_membrane_as USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__epi_frag_annotation_stop__id
    ON public.epitope_186540_membrane_as USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__host_taxon_id__idx
    ON public.epitope_186540_membrane_as USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__host_taxon_name_lower__idx
    ON public.epitope_186540_membrane_as USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__iedb_epitope_id__idx
    ON public.epitope_186540_membrane_as USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__is_linear__idx
    ON public.epitope_186540_membrane_as USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__mhc_allele__idx
    ON public.epitope_186540_membrane_as USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__mhc_class_lower__idx
    ON public.epitope_186540_membrane_as USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__product_lower__idx
    ON public.epitope_186540_membrane_as USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__response_frequency_pos__idx
    ON public.epitope_186540_membrane_as USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__sequence_aa_alternative__idx
    ON public.epitope_186540_membrane_as USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__sequence_aa_original__idx
    ON public.epitope_186540_membrane_as USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__start_aa_original__idx
    ON public.epitope_186540_membrane_as USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__taxon_id__idx
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__taxon_name_lower__idx
    ON public.epitope_186540_membrane_as USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__variant_aa_length__idx
    ON public.epitope_186540_membrane_as USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__variant_aa_type__idx
    ON public.epitope_186540_membrane_as USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_cell_type__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_epi_start__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_epi_stop__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_is_linear__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_mhc_allele__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_product__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_membrane_as__virus_host_resp_freq__i
    ON public.epitope_186540_membrane_as USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


-- MATERIALIZED VIEW AND INDEXES FOR VIRUS 186540 AND PROTEIN RNA-dependent RNA polymerase
CREATE MATERIALIZED VIEW public.epitope_186540_rna_depende
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
             AND vir.taxon_id = 186540)
  ORDER BY epi.iedb_epitope_id
WITH DATA;

ALTER TABLE public.epitope_186540_rna_depende
    OWNER TO geco;

CREATE INDEX epi_186540_rna_depende__cell_type__idx
    ON public.epitope_186540_rna_depende USING btree
    ((cell_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__epi_annotation_start__idx
    ON public.epitope_186540_rna_depende USING btree
    (epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__epi_annotation_stop__idx
    ON public.epitope_186540_rna_depende USING btree
    (epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__epi_frag_annotation_start__i
    ON public.epitope_186540_rna_depende USING btree
    (epi_frag_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__epi_frag_annotation_stop__id
    ON public.epitope_186540_rna_depende USING btree
    (epi_frag_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__host_taxon_id__idx
    ON public.epitope_186540_rna_depende USING btree
    (host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__host_taxon_name_lower__idx
    ON public.epitope_186540_rna_depende USING btree
    (lower(host_taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__iedb_epitope_id__idx
    ON public.epitope_186540_rna_depende USING btree
    (iedb_epitope_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__is_linear__idx
    ON public.epitope_186540_rna_depende USING btree
    (is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__mhc_allele__idx
    ON public.epitope_186540_rna_depende USING btree
    ((mhc_allele::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__mhc_class_lower__idx
    ON public.epitope_186540_rna_depende USING btree
    (lower(mhc_class::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__product_lower__idx
    ON public.epitope_186540_rna_depende USING btree
    (lower(product::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__response_frequency_pos__idx
    ON public.epitope_186540_rna_depende USING btree
    (response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__sequence_aa_alternative__idx
    ON public.epitope_186540_rna_depende USING btree
    ((sequence_aa_alternative::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__sequence_aa_original__idx
    ON public.epitope_186540_rna_depende USING btree
    ((sequence_aa_original::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__start_aa_original__idx
    ON public.epitope_186540_rna_depende USING btree
    (start_aa_original)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__taxon_id__idx
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__taxon_name_lower__idx
    ON public.epitope_186540_rna_depende USING btree
    (lower(taxon_name::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__variant_aa_length__idx
    ON public.epitope_186540_rna_depende USING btree
    (variant_aa_length)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__variant_aa_type__idx
    ON public.epitope_186540_rna_depende USING btree
    ((variant_aa_type::text) COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_taxon_and_host_taxon_id__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_cell_type__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, cell_type COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_epi_start__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_start)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_epi_stop__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, epi_annotation_stop)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_is_linear__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, is_linear)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_mhc_allele__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, mhc_allele COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_product__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, product COLLATE pg_catalog."default")
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;
CREATE INDEX epi_186540_rna_depende__virus_host_resp_freq__i
    ON public.epitope_186540_rna_depende USING btree
    (taxon_id, host_taxon_id, response_frequency_pos)
    WITH (FILLFACTOR=100)
    TABLESPACE default_ts;


