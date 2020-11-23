-- public.db_meta definition

-- Drop table

-- DROP TABLE db_meta;

CREATE TABLE db_meta (
	virus_id serial primary key ,
	date_of_import varchar NULL
);


-- public.experiment_type definition

-- Drop table

-- DROP TABLE experiment_type;

CREATE TABLE experiment_type (
	experiment_type_id serial primary key,
	sequencing_technology varchar NULL,
	assembly_method varchar NULL,
	coverage int4 NULL
);


-- public.sequencing_project definition

-- Drop table

-- DROP TABLE sequencing_project;

CREATE TABLE sequencing_project (
	sequencing_project_id serial primary key ,
	sequencing_lab varchar NULL,
	submission_date date NULL,
	database_source varchar NULL,
	bioproject_id varchar NULL
);


-- public.virus definition

-- Drop table

-- DROP TABLE virus;

CREATE TABLE virus (
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


-- public.epitope definition

-- Drop table

-- DROP TABLE epitope;

CREATE TABLE epitope (
	epitope_id serial primary key,
	virus_id int4 NOT NULL,
	host_id int4 NOT NULL,
	protein_ncbi_id varchar NULL,
	cell_type varchar NULL,
	mhc_class varchar NULL,
	response_frequency float4 NULL,
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


-- public.epitope_fragment definition

-- Drop table

-- DROP TABLE epitope_fragment;

CREATE TABLE epitope_fragment (
	epi_fragment_id serial primary key,
	epitope_id int4 NULL,
	epi_fragment_sequence varchar NULL,
	epi_frag_annotation_start int4 NULL,
	epi_frag_annotation_stop int4 NULL
);


-- public.host_specie definition

-- Drop table

-- DROP TABLE host_specie;

CREATE TABLE host_specie (
	host_id serial primary key,
	host_taxon_id int4 NULL,
	host_taxon_name varchar NULL
);


-- public.host_sample definition

-- Drop table

-- DROP TABLE host_sample;

CREATE TABLE host_sample (
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


-- public."sequence" definition

-- Drop table

-- DROP TABLE "sequence";

CREATE TABLE "sequence" (
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
	nucleotide_sequence varchar NULL,
	strand varchar NULL,
	length int4 NULL,
	gc_percentage float8 NULL,
	n_percentage float8 NULL,
	lineage varchar NULL,
	clade varchar NULL,
	gisaid_only bool NOT NULL
);


-- public.nucleotide_variant definition

-- Drop table

-- DROP TABLE nucleotide_variant;

CREATE TABLE nucleotide_variant (
	nucleotide_variant_id serial primary key,
	sequence_id int4 NOT NULL,
	sequence_original varchar NOT NULL,
	sequence_alternative varchar NOT NULL,
	start_original int4 NULL,
	start_alternative int4 NULL,
	variant_length int4 NOT NULL,
	variant_type varchar NOT NULL
);


-- public.variant_impact definition

-- Drop table

-- DROP TABLE variant_impact;

CREATE TABLE variant_impact (
	variant_impact_id bigserial primary key,
	nucleotide_variant_id int4 NOT NULL,
	effect varchar NULL,
	putative_impact varchar NULL,
	impact_gene_name varchar NULL
);


-- public.annotation definition

-- Drop table

-- DROP TABLE annotation;

CREATE TABLE annotation (
	annotation_id serial primary key,
	sequence_id int4 NOT NULL,
	feature_type varchar NOT NULL,
	"start" int4 NULL,
	stop int4 NULL,
	gene_name varchar NULL,
	product varchar NULL,
	external_reference varchar NULL,
	aminoacid_sequence varchar NULL,
	annotation_nucleotide_sequence varchar NULL
);


-- public.aminoacid_variant definition

-- Drop table

-- DROP TABLE aminoacid_variant;

CREATE TABLE aminoacid_variant (
	aminoacid_variant_id serial primary key,
	annotation_id int4 NOT NULL,
	sequence_aa_original varchar NOT NULL,
	sequence_aa_alternative varchar NOT NULL,
	start_aa_original int4 NULL,
	variant_aa_length int4 NOT NULL,
	variant_aa_type varchar NOT NULL
);