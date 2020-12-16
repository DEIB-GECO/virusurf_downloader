------------------- INDEXES
--SEQUENCE
DROP INDEX IF EXISTS seq__experiment_id;
DROP INDEX IF EXISTS seq__host_id;
DROP INDEX IF EXISTS seq__seq_proj_id;
DROP INDEX IF EXISTS seq__virus_id;
DROP INDEX IF EXISTS sequence__is_reference__idx;
DROP INDEX IF EXISTS seq__accession_id;
DROP INDEX IF EXISTS seq__alternative_accession_id;

--NUCLEOTIDE VARIANT
DROP INDEX IF EXISTS nuc_var__length;
DROP INDEX IF EXISTS nuc_var__seq_id;
DROP INDEX IF EXISTS nuc_var__start_alt;
DROP INDEX IF EXISTS nuc_var__start_orig;

--VARIANT IMPACT
DROP INDEX IF EXISTS impact__var_id;

--ANNOTATION
DROP INDEX IF EXISTS ann__seq_id;
DROP INDEX IF EXISTS ann__start;
DROP INDEX IF EXISTS ann__stop;

--AMINO ACID VARIANT
DROP INDEX IF EXISTS aa__ann_id;
DROP INDEX IF EXISTS aa__start_original;
DROP INDEX IF EXISTS aa__var_type_lower;
DROP INDEX IF EXISTS aa__var_type_normal;


--NUCLEOTIDE VARIANT ANNOTATION
DROP INDEX IF EXISTS nuc_var_ann__var_id;


--NUCLEOTIDE VARIANT ANNOTATED
DROP INDEX IF EXISTS nucleotide_variant_annotated__sequence_id__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__nucleotide_variant_id__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__variant_type_lower__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__start_original__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__sequence_original_lower__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__sequence_alternative_lower__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__n_feature_type_lower__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__n_gene_name_lower__idx;
DROP INDEX IF EXISTS nucleotide_variant_annotated__n_product__idx;


------------------  VIEWS
DROP MATERIALIZED VIEW IF EXISTS nucleotide_variant_annotation;
DROP MATERIALIZED VIEW IF EXISTS nucleotide_variant_annotated;

DROP VIEW IF EXISTS annotation_cds;
DROP VIEW IF EXISTS annotation_view;
DROP VIEW IF EXISTS host_sample_view;
DROP VIEW IF EXISTS nucleotide_variant_limited;


-------------------- CONSTRAINTS
-- EPITOPE
ALTER TABLE epitope DROP CONSTRAINT IF EXISTS epitope_host_id_fkey;
ALTER TABLE epitope DROP CONSTRAINT IF EXISTS epitope_virus_id_fkey;

--EPITOPE FRAGMENT
ALTER TABLE epitope_fragment DROP CONSTRAINT IF EXISTS epitope_fragment_epitope_id_fkey;

--HOST SAMPLE
ALTER TABLE host_sample DROP CONSTRAINT IF EXISTS host_sample_host_id_fkey;

--SEQUENCE
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_accession_id_key;
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_alternative_accession_id_key;
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_experiment_type_id_fkey;
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_host_sample_id_fkey;
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_sequencing_project_id_fkey;
ALTER TABLE sequence DROP CONSTRAINT IF EXISTS sequence_virus_id_fkey;

--ANNOTATION
ALTER TABLE annotation DROP CONSTRAINT IF EXISTS annotation_sequence_id_fkey;

--NUCLEOTIDE VARIANT
ALTER TABLE nucleotide_variant DROP CONSTRAINT IF EXISTS nucleotide_variant_sequence_id_fkey;

--VARIANT IMPACT
ALTER TABLE variant_impact DROP CONSTRAINT IF EXISTS variant_impact_nucleotide_variant_id_fkey;

--AMINO ACID VARIANT
ALTER TABLE aminoacid_variant DROP CONSTRAINT IF EXISTS aminoacid_variant_annotation_id_fkey;

-- DB META
ALTER TABLE db_meta DROP CONSTRAINT IF EXISTS db_meta_pkey;