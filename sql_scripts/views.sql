-- CREATE OR REPLACE VIEW public.annotation_cds
-- AS SELECT annotation.annotation_id,
--     annotation.sequence_id,
--     annotation.start,
--     annotation.stop,
--     annotation.gene_name,
--     annotation.product,
--     annotation.external_reference,
--     annotation.aminoacid_sequence
--    FROM annotation
--   WHERE annotation.feature_type::text = 'CDS'::text;
--
--
-- CREATE OR REPLACE VIEW public.annotation_view
-- AS SELECT annotation.sequence_id,
--     annotation.product AS annotation_view_product,
--     annotation.aminoacid_sequence AS annotation_view_aminoacid_sequence,
--     annotation.annotation_nucleotide_sequence AS annotation_view_nucleotide_sequence
--    FROM annotation
--   WHERE annotation.product IS NOT NULL AND
--         (annotation.aminoacid_sequence IS NOT NULL OR annotation.annotation_nucleotide_sequence IS NOT NULL);


CREATE OR REPLACE VIEW public.host_sample_view
AS SELECT host_sample.*,
    host_specie.host_taxon_id, host_specie.host_taxon_name
   FROM host_sample
     JOIN host_specie USING (host_id);


CREATE OR REPLACE VIEW public.nucleotide_variant_limited
AS SELECT nucleotide_variant.*
   FROM nucleotide_variant
  WHERE nucleotide_variant.variant_length <= 20;


DROP MATERIALIZED VIEW IF EXISTS nucleotide_variant_annotated;
CREATE MATERIALIZED VIEW public.nucleotide_variant_annotated
WITH (
    FILLFACTOR = 100
)
TABLESPACE default_ts
AS
SELECT nc.nucleotide_variant_id,
    nc.sequence_id,
    nc.variant_type,
    nc.start_original,
    CASE nc.variant_type  when 'SUB' then  nc.sequence_original    ELSE null END as sequence_original,
    CASE nc.variant_type  when 'SUB' then  nc.sequence_alternative ELSE null END as sequence_alternative,
    nc.variant_length,
    ann.gene_name AS n_gene_name
   FROM nucleotide_variant nc
     LEFT JOIN annotation ann ON nc.start_original >= ann.start AND nc.start_original <= ann.stop AND nc.sequence_id = ann.sequence_id AND ann.feature_type = 'gene'
WITH DATA;
