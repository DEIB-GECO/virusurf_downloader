CREATE OR REPLACE VIEW public.annotation_cds
AS SELECT annotation.annotation_id,
    annotation.sequence_id,
    annotation.start,
    annotation.stop,
    annotation.gene_name,
    annotation.product,
    annotation.external_reference,
    annotation.aminoacid_sequence
   FROM annotation
  WHERE annotation.feature_type::text = 'CDS'::text;


CREATE OR REPLACE VIEW public.annotation_view
AS SELECT annotation.sequence_id,
    annotation.product AS annotation_view_product,
    annotation.aminoacid_sequence AS annotation_view_aminoacid_sequence,
    annotation.annotation_nucleotide_sequence AS annotation_view_nucleotide_sequence
   FROM annotation
  WHERE annotation.product IS NOT NULL AND
        (annotation.aminoacid_sequence IS NOT NULL OR annotation.annotation_nucleotide_sequence IS NOT NULL);


CREATE OR REPLACE VIEW public.host_sample_view
AS SELECT host_sample.*,
    host_specie.*
   FROM host_sample
     JOIN host_specie USING (host_id);


CREATE OR REPLACE VIEW public.nucleotide_variant_limited
AS SELECT nucleotide_variant.*
   FROM nucleotide_variant
  WHERE nucleotide_variant.variant_length <= 20;


CREATE MATERIALIZED VIEW public.nucleotide_variant_annotation
TABLESPACE default_ts
AS SELECT nucleotide_variant.nucleotide_variant_id,
    annotation.feature_type AS n_feature_type,
    annotation.gene_name AS n_gene_name,
    annotation.product AS n_product
   FROM annotation JOIN nucleotide_variant
         ON nucleotide_variant.start_original >= annotation.start
                AND nucleotide_variant.start_original <= annotation.stop
                AND nucleotide_variant.sequence_id = annotation.sequence_id
WITH DATA;


CREATE MATERIALIZED VIEW nucleotide_variant_annotated
WITH (FILLFACTOR = 100)
AS
 SELECT DISTINCT
  nc.nucleotide_variant_id,
  -- ann.annotation_id, -- --> possibly have annotation_id, too. but for now, we an skip that.

  nc.sequence_id,
  nc.variant_type,
  nc.start_original,
  nc.sequence_original,
  nc.sequence_alternative,
  nc.variant_length,

    ann.feature_type AS n_feature_type,
    ann.gene_name AS n_gene_name,
    ann.product AS n_product

   FROM nucleotide_variant as nc
     LEFT JOIN annotation as ann  -- --> LEFT?
   ON nc.start_original >= ann.start
  AND nc.start_original <= ann.stop
  AND nc.sequence_id = ann.sequence_id
  WHERE variant_length <= 20  -----------------> ?remove LIMIT 20
  -- and nc.sequence_id < 10 and ann.sequence_id < 10 -- FOR TESTING
  ;
