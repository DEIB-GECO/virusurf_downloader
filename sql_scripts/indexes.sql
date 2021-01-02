--SEQUENCE
CREATE INDEX seq__experiment_id ON sequence USING btree (experiment_type_id) TABLESPACE default_ts;
CREATE INDEX seq__host_id ON sequence USING btree (host_sample_id) TABLESPACE default_ts;
CREATE INDEX seq__seq_proj_id ON sequence USING btree (sequencing_project_id) TABLESPACE default_ts;
CREATE INDEX seq__virus_id ON sequence USING btree (virus_id) TABLESPACE default_ts;
CREATE INDEX sequence__is_reference__idx ON sequence USING btree (is_reference) TABLESPACE default_ts;
CREATE UNIQUE INDEX seq__accession_id ON sequence USING btree (lower((accession_id)::text)) TABLESPACE default_ts;
CREATE UNIQUE INDEX seq__alternative_accession_id ON sequence USING btree (lower((alternative_accession_id)::text)) TABLESPACE default_ts;

--NUCLEOTIDE VARIANT
CREATE INDEX nuc_var__length ON nucleotide_variant USING btree (variant_length) TABLESPACE default_ts;
CREATE INDEX nuc_var__seq_id ON nucleotide_variant USING btree (sequence_id) TABLESPACE default_ts;
CREATE INDEX nuc_var__start_alt ON nucleotide_variant USING btree (start_alternative) TABLESPACE default_ts;
CREATE INDEX nuc_var__start_orig ON nucleotide_variant USING btree (start_original) TABLESPACE default_ts;

--VARIANT IMPACT
CREATE INDEX impact__var_id ON variant_impact USING btree (nucleotide_variant_id) TABLESPACE default_ts;

--ANNOTATION
CREATE INDEX ann__seq_id ON annotation USING btree (sequence_id) TABLESPACE default_ts;
CREATE INDEX ann__start ON annotation USING btree (start) TABLESPACE default_ts;
CREATE INDEX ann__stop ON annotation USING btree (stop) TABLESPACE default_ts;

-- ANNOTATION SEQUENCE
CREATE UNIQUE INDEX ann_seq__seq_id__product ON annotation_sequence(sequence_id, product) TABLESPACE default_ts;

--AMINO ACID VARIANT
CREATE INDEX aa__ann_id ON aminoacid_variant USING btree (annotation_id) TABLESPACE default_ts;
CREATE INDEX aa__start_original ON aminoacid_variant USING btree (start_aa_original) TABLESPACE default_ts;
CREATE INDEX aa__var_type_lower ON aminoacid_variant USING btree (lower((variant_aa_type)::text)) TABLESPACE default_ts;
CREATE INDEX aa__var_type_normal ON aminoacid_variant USING btree (variant_aa_type) TABLESPACE default_ts;

--NUCLEOTIDE VARIANT ANNOTATED
-- JOIN KEY
CREATE INDEX nucleotide_variant_annotated__sequence_id__idx ON nucleotide_variant_annotated(sequence_id) WITH (FILLFACTOR = 100) TABLESPACE default_ts;

CLUSTER VERBOSE nucleotide_variant_annotated USING nucleotide_variant_annotated__sequence_id__idx;

CREATE INDEX nucleotide_variant_annotated__nucleotide_variant_id__idx ON public.nucleotide_variant_annotated USING btree (nucleotide_variant_id) WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__sequence_alternative_lower__idx ON public.nucleotide_variant_annotated USING btree (lower(sequence_alternative::text) COLLATE pg_catalog."default") WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__sequence_id__idx ON public.nucleotide_variant_annotated USING btree (sequence_id) WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__sequence_original_lower__idx ON public.nucleotide_variant_annotated USING btree (lower(sequence_original::text) COLLATE pg_catalog."default") WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__start_original__idx ON public.nucleotide_variant_annotated USING btree (start_original) WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__variant_type_lower__idx ON public.nucleotide_variant_annotated USING btree (lower(variant_type::text) COLLATE pg_catalog."default") WITH (FILLFACTOR=100) TABLESPACE default_ts;
CREATE INDEX nucleotide_variant_annotated__n_gene_name_lower__idx ON public.nucleotide_variant_annotated USING btree (lower(n_gene_name::text) COLLATE pg_catalog."default") WITH (FILLFACTOR=100) TABLESPACE default_ts;