--EPITOPE
ALTER TABLE epitope ADD CONSTRAINT epitope_host_id_fkey FOREIGN KEY (host_id) REFERENCES host_specie(host_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;
ALTER TABLE epitope ADD CONSTRAINT epitope_virus_id_fkey FOREIGN KEY (virus_id) REFERENCES virus(virus_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--EPITOPE FRAGMENT 
ALTER TABLE epitope_fragment ADD CONSTRAINT epitope_fragment_epitope_id_fkey FOREIGN KEY (epitope_id) REFERENCES epitope(epitope_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--HOST SAMPLE
ALTER TABLE host_sample ADD CONSTRAINT host_sample_host_id_fkey FOREIGN KEY (host_id) REFERENCES host_specie(host_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--SEQUENCE
-- ALTER TABLE sequence ADD CONSTRAINT sequence_accession_id_key UNIQUE (accession_id); -- constraint replaced by unique index on lower(...)
-- ALTER TABLE sequence ADD CONSTRAINT sequence_alternative_accession_id_key UNIQUE (alternative_accession_id); -- constraint replaced by unique index on lower(...)
ALTER TABLE sequence ADD CONSTRAINT sequence_experiment_type_id_fkey FOREIGN KEY (experiment_type_id) REFERENCES experiment_type(experiment_type_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;
ALTER TABLE sequence ADD CONSTRAINT sequence_host_sample_id_fkey FOREIGN KEY (host_sample_id) REFERENCES host_sample(host_sample_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;
ALTER TABLE sequence ADD CONSTRAINT sequence_sequencing_project_id_fkey FOREIGN KEY (sequencing_project_id) REFERENCES sequencing_project(sequencing_project_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;
ALTER TABLE sequence ADD CONSTRAINT sequence_virus_id_fkey FOREIGN KEY (virus_id) REFERENCES virus(virus_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--ANNOTATION
ALTER TABLE annotation ADD CONSTRAINT annotation_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--NUCLEOTIDE VARIANT
ALTER TABLE nucleotide_variant ADD CONSTRAINT nucleotide_variant_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--VARIANT IMPACT
ALTER TABLE variant_impact ADD CONSTRAINT variant_impact_nucleotide_variant_id_fkey FOREIGN KEY (nucleotide_variant_id) REFERENCES nucleotide_variant(nucleotide_variant_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;

--AMINO ACID VARIANT
ALTER TABLE aminoacid_variant ADD CONSTRAINT aminoacid_variant_annotation_id_fkey FOREIGN KEY (annotation_id) REFERENCES annotation(annotation_id) MATCH SIMPLE
ON UPDATE NO ACTION
ON DELETE NO ACTION;
