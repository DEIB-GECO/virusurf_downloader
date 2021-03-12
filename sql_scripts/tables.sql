CREATE TABLE public.db_meta (
	virus_id int4 NOT NULL,
	"source" varchar NOT NULL,
	date_of_import date NULL,
	PRIMARY KEY (virus_id, "source")
);

CREATE TABLE public.experiment_type (
	experiment_type_id serial primary key,
	sequencing_technology varchar NULL,
	assembly_method varchar NULL,
	coverage int4 NULL
);

CREATE TABLE public.sequencing_project (
	sequencing_project_id serial primary key ,
	sequencing_lab varchar NULL,
	submission_date date NULL,
	database_source varchar NULL,
	bioproject_id varchar NULL
);

CREATE TABLE public.virus (
	virus_id serial primary key,
	taxon_id int4 NULL,
	taxon_name varchar NULL,
	"family" varchar NULL,
	sub_family varchar NULL,
	genus varchar NULL,
	species varchar NULL,
	equivalent_list varchar NULL,
	molecule_type varchar NULL,
	is_single_stranded bool NULL,
	is_positive_stranded bool NULL
);

CREATE TABLE public.epitope (
	epitope_id serial primary key,
	epitope_iri varchar NULL,
	iedb_epitope_id int4 NULL,
	virus_id int4 NOT NULL,
	host_id int4 NOT NULL,
	protein_ncbi_id varchar NULL,
	cell_type varchar NULL,
	mhc_class varchar NULL,
	response_frequency_pos float4 NULL,
	epitope_sequence varchar NULL,
	epi_annotation_start int4 NULL,
	epi_annotation_stop int4 NULL,
	external_link varchar NULL,
	prediction_process varchar NULL,
	is_linear bool NULL,
	source_host_name varchar NULL,
	source_host_iri varchar NULL,
	mhc_allele varchar NULL,
	assay_type varchar NULL
);

CREATE TABLE public.epitope_fragment (
	epi_fragment_id serial primary key,
	epitope_id int4 NULL,
	epi_fragment_sequence varchar NULL,
	epi_frag_annotation_start int4 NULL,
	epi_frag_annotation_stop int4 NULL
);

CREATE TABLE public.host_specie (
	host_id serial primary key,
	host_taxon_id int4 NULL,
	host_taxon_name varchar NULL
);

CREATE TABLE public.host_sample (
	host_sample_id serial primary key,
	host_id int4 NULL,
	collection_date varchar NULL,
	isolation_source varchar NULL,
	originating_lab varchar NULL,
	country varchar NULL,
	region varchar NULL,
	geo_group varchar NULL,
	age int4 NULL,
	gender varchar NULL
);

CREATE TABLE public."sequence" (
	sequence_id serial primary key,
	experiment_type_id int4 NOT NULL,
	virus_id int4 NOT NULL,
	host_sample_id int4 NOT NULL,
	sequencing_project_id int4 NOT NULL,
	accession_id varchar NOT NULL,
	alternative_accession_id varchar NULL,
	strain_name varchar NULL,
	is_reference bool NOT NULL,
	is_complete bool NULL,
	strand varchar NULL,
	length int4 NULL,
	gc_percentage float8 NULL,
	n_percentage float8 NULL,
	lineage varchar NULL,
	clade varchar NULL,
	gisaid_only bool NOT NULL
);

CREATE TABLE public.nucleotide_sequence (
    sequence_id int4 primary key,
    nucleotide_sequence varchar NULL
);

CREATE TABLE public.nucleotide_variant (
	nucleotide_variant_id serial primary key,
	sequence_id int4 NOT NULL,
	sequence_original varchar NOT NULL,
	sequence_alternative varchar NOT NULL,
	start_original int4 NULL,
	start_alternative int4 NULL,
	variant_length int4 NOT NULL,
	variant_type varchar NOT NULL
);

CREATE TABLE public.variant_impact (
	variant_impact_id serial primary key,
	nucleotide_variant_id int4 NOT NULL,
	effect varchar NULL,
	putative_impact varchar NULL,
	impact_gene_name varchar NULL
);

CREATE TABLE public.annotation (
	annotation_id serial primary key,
	sequence_id int4 NOT NULL,
	feature_type varchar NOT NULL,
	"start" int4 NULL,
	stop int4 NULL,
	gene_name varchar NULL,
	product varchar NULL,
	external_reference varchar NULL
);

CREATE TABLE public.annotation_sequence (
    annotation_id int4 primary key,
	sequence_id int4 NOT NULL,
	product varchar NULL,
	aminoacid_sequence varchar NULL,
	annotation_nucleotide_sequence varchar NULL
);

CREATE TABLE public.aminoacid_variant (
	aminoacid_variant_id serial primary key,
	annotation_id int4 NOT NULL,
	sequence_aa_original varchar NOT NULL,
	sequence_aa_alternative varchar NOT NULL,
	start_aa_original int4 NULL,
	variant_aa_length int4 NOT NULL,
	variant_aa_type varchar NOT NULL
);

CREATE TABLE public.overlap (
    sequence_id int4 NOT NULL,
    accession_id varchar NOT NULL,
    overlapping_accession_id varchar NOT NULL,
    overlapping_source varchar NOT NULL
);

CREATE TABLE public.pipeline_event (
    event_id    serial primary key,
    event_name    varchar NOT NULL,
    event_date    varchar,
    added_items   integer,
    removed_items integer,
    changed_items integer
);