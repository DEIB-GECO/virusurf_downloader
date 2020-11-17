-- THIS CONSTRAINT HAS BEEN REPLACED BY AN INDEX ON LOWER(SEQUENCE.ACCESSION_ID)
-- Constraint: sequence_accession_id_key
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_accession_id_key;
-- ALTER TABLE public.sequence
--     ADD CONSTRAINT sequence_accession_id_key UNIQUE (accession_id)
--     USING INDEX TABLESPACE default_ts;


-- THIS CONSTRAINT HAS BEEN REPLACED BY AN INDEX ON LOWER(SEQUENCE.ALTERNATIVE_ACCESSION_ID)
-- Constraint: sequence_alternative_accession_id_key
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_alternative_accession_id_key;
-- ALTER TABLE public.sequence
--     ADD CONSTRAINT sequence_alternative_accession_id_key UNIQUE (alternative_accession_id)
--     USING INDEX TABLESPACE default_ts;



-- Constraint: aminoacid_variant_annotation_id_fkey
-- ALTER TABLE public.aminoacid_variant DROP CONSTRAINT aminoacid_variant_annotation_id_fkey;
ALTER TABLE public.aminoacid_variant
    ADD CONSTRAINT aminoacid_variant_annotation_id_fkey FOREIGN KEY (annotation_id)
    REFERENCES public.annotation (annotation_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: annotation_sequence_id_fkey
-- ALTER TABLE public.annotation DROP CONSTRAINT annotation_sequence_id_fkey;
ALTER TABLE public.annotation
    ADD CONSTRAINT annotation_sequence_id_fkey FOREIGN KEY (sequence_id)
    REFERENCES public.sequence (sequence_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: variant_impact_nucleotide_variant_id_fkey
-- ALTER TABLE public.variant_impact DROP CONSTRAINT variant_impact_nucleotide_variant_id_fkey;
ALTER TABLE public.variant_impact
    ADD CONSTRAINT variant_impact_nucleotide_variant_id_fkey FOREIGN KEY (nucleotide_variant_id)
    REFERENCES public.nucleotide_variant (nucleotide_variant_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;


-- Constraint: nucleotide_variant_sequence_id_fkey
-- ALTER TABLE public.nucleotide_variant DROP CONSTRAINT nucleotide_variant_sequence_id_fkey;
ALTER TABLE public.nucleotide_variant
    ADD CONSTRAINT nucleotide_variant_sequence_id_fkey FOREIGN KEY (sequence_id)
    REFERENCES public.sequence (sequence_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;


-- Constraint: sequence_experiment_type_id_fkey
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_experiment_type_id_fkey;
ALTER TABLE public.sequence
    ADD CONSTRAINT sequence_experiment_type_id_fkey FOREIGN KEY (experiment_type_id)
    REFERENCES public.experiment_type (experiment_type_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: sequence_host_sample_id_fkey
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_host_sample_id_fkey;
ALTER TABLE public.sequence
    ADD CONSTRAINT sequence_host_sample_id_fkey FOREIGN KEY (host_sample_id)
    REFERENCES public.host_sample (host_sample_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: sequence_sequencing_project_id_fkey
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_sequencing_project_id_fkey;
ALTER TABLE public.sequence
    ADD CONSTRAINT sequence_sequencing_project_id_fkey FOREIGN KEY (sequencing_project_id)
    REFERENCES public.sequencing_project (sequencing_project_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: sequence_virus_id_fkey
-- ALTER TABLE public.sequence DROP CONSTRAINT sequence_virus_id_fkey;
ALTER TABLE public.sequence
    ADD CONSTRAINT sequence_virus_id_fkey FOREIGN KEY (virus_id)
    REFERENCES public.virus (virus_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: host_sample_host_id_fkey
-- ALTER TABLE public.host_sample DROP CONSTRAINT host_sample_host_id_fkey;
ALTER TABLE public.host_sample
    ADD CONSTRAINT host_sample_host_id_fkey FOREIGN KEY (host_id)
    REFERENCES public.host_specie (host_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: epitope_host_id_fkey
-- ALTER TABLE public.epitope DROP CONSTRAINT epitope_host_id_fkey;
ALTER TABLE public.epitope
    ADD CONSTRAINT epitope_host_id_fkey FOREIGN KEY (host_id)
    REFERENCES public.host_specie (host_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: epitope_virus_id_fkey
-- ALTER TABLE public.epitope DROP CONSTRAINT epitope_virus_id_fkey;
ALTER TABLE public.epitope
    ADD CONSTRAINT epitope_virus_id_fkey FOREIGN KEY (virus_id)
    REFERENCES public.virus (virus_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;



-- Constraint: epitope_fragment_epitope_id_fkey
-- ALTER TABLE public.epitope_fragment DROP CONSTRAINT epitope_fragment_epitope_id_fkey;
ALTER TABLE public.epitope_fragment
    ADD CONSTRAINT epitope_fragment_epitope_id_fkey FOREIGN KEY (epitope_id)
    REFERENCES public.epitope (epitope_id) MATCH SIMPLE
    ON UPDATE NO ACTION
    ON DELETE NO ACTION;

