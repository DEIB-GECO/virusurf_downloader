DROP MATERIALIZED VIEW IF EXISTS nucleotide_variant_annotation;
DROP MATERIALIZED VIEW IF EXISTS nucleotide_variant_annotated;

DROP VIEW IF EXISTS annotation_cds;
DROP VIEW IF EXISTS annotation_view;
DROP VIEW IF EXISTS host_sample_view;
DROP VIEW IF EXISTS nucleotide_variant_limited;

TRUNCATE db_meta, aminoacid_variant, annotation, variant_impact, nucleotide_variant, sequencing_project, host_specie,
    host_sample, experiment_type, virus, sequence RESTART IDENTITY;