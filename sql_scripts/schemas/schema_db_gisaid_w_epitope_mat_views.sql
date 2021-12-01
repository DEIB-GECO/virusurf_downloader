--
-- PostgreSQL database dump
--

-- Dumped from database version 10.18 (Ubuntu 10.18-0ubuntu0.18.04.1)
-- Dumped by pg_dump version 12.5 (Ubuntu 12.5-1.pgdg16.04+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: btree_gin; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS btree_gin WITH SCHEMA public;


--
-- Name: EXTENSION btree_gin; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION btree_gin IS 'support for indexing common datatypes in GIN';


--
-- Name: btree_gist; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS btree_gist WITH SCHEMA public;


--
-- Name: EXTENSION btree_gist; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION btree_gist IS 'support for indexing common datatypes in GiST';


--
-- Name: pg_trgm; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pg_trgm WITH SCHEMA public;


--
-- Name: EXTENSION pg_trgm; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION pg_trgm IS 'text similarity measurement and index searching based on trigrams';


SET default_tablespace = '';

--
-- Name: aminoacid_variant; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.aminoacid_variant (
    aminoacid_variant_id integer NOT NULL,
    annotation_id integer NOT NULL,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    start_aa_original integer,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.aminoacid_variant OWNER TO geco;

--
-- Name: aminoacid_variant_aminoacid_variant_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.aminoacid_variant_aminoacid_variant_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.aminoacid_variant_aminoacid_variant_id_seq OWNER TO geco;

--
-- Name: aminoacid_variant_aminoacid_variant_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.aminoacid_variant_aminoacid_variant_id_seq OWNED BY public.aminoacid_variant.aminoacid_variant_id;


--
-- Name: annotation; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.annotation (
    annotation_id integer NOT NULL,
    sequence_id integer NOT NULL,
    feature_type character varying NOT NULL,
    start integer,
    stop integer,
    gene_name character varying,
    product character varying,
    external_reference character varying
);


ALTER TABLE public.annotation OWNER TO geco;

--
-- Name: annotation_annotation_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.annotation_annotation_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.annotation_annotation_id_seq OWNER TO geco;

--
-- Name: annotation_annotation_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.annotation_annotation_id_seq OWNED BY public.annotation.annotation_id;


--
-- Name: annotation_sequence; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.annotation_sequence (
    annotation_id integer NOT NULL,
    sequence_id integer NOT NULL,
    product character varying,
    aminoacid_sequence character varying,
    annotation_nucleotide_sequence character varying
);


ALTER TABLE public.annotation_sequence OWNER TO geco;

--
-- Name: db_meta; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.db_meta (
    virus_id integer NOT NULL,
    source character varying NOT NULL,
    date_of_import date
);


ALTER TABLE public.db_meta OWNER TO geco;

--
-- Name: epitope; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope (
    epitope_id integer NOT NULL,
    epitope_iri character varying,
    iedb_epitope_id integer,
    virus_id integer NOT NULL,
    host_id integer NOT NULL,
    protein_name character varying,
    cell_type character varying,
    mhc_class character varying,
    response_frequency_pos real,
    epitope_sequence character varying,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    external_link character varying,
    prediction_process character varying,
    is_linear boolean,
    source_host_name character varying,
    source_host_iri character varying,
    mhc_allele character varying,
    assay_type character varying
);


ALTER TABLE public.epitope OWNER TO geco;

--
-- Name: epitope_fragment; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_fragment (
    epi_fragment_id integer NOT NULL,
    epitope_id integer,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer
);


ALTER TABLE public.epitope_fragment OWNER TO geco;

--
-- Name: host_sample; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.host_sample (
    host_sample_id integer NOT NULL,
    host_id integer,
    collection_date character varying,
    isolation_source character varying,
    originating_lab character varying,
    region character varying,
    country character varying,
    geo_group character varying,
    age integer,
    gender character varying,
    province character varying,
    coll_date_precision smallint
);


ALTER TABLE public.host_sample OWNER TO geco;

--
-- Name: host_specie; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.host_specie (
    host_id integer NOT NULL,
    host_taxon_id integer,
    host_taxon_name character varying
);


ALTER TABLE public.host_specie OWNER TO geco;

--
-- Name: sequence; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.sequence (
    sequence_id integer NOT NULL,
    experiment_type_id integer NOT NULL,
    virus_id integer NOT NULL,
    host_sample_id integer NOT NULL,
    sequencing_project_id integer NOT NULL,
    accession_id character varying NOT NULL,
    alternative_accession_id character varying,
    strain_name character varying,
    is_reference boolean NOT NULL,
    is_complete boolean,
    strand character varying,
    length integer,
    gc_percentage double precision,
    n_percentage double precision,
    lineage character varying,
    clade character varying,
    gisaid_only boolean NOT NULL,
    nucleotide_sequence character varying
);


ALTER TABLE public.sequence OWNER TO geco;

--
-- Name: virus; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.virus (
    virus_id integer NOT NULL,
    taxon_id integer,
    taxon_name character varying,
    family character varying,
    sub_family character varying,
    genus character varying,
    species character varying,
    equivalent_list character varying,
    molecule_type character varying,
    is_single_stranded boolean,
    is_positive_stranded boolean
);


ALTER TABLE public.virus OWNER TO geco;

--
-- Name: epitope_2697049_e_envelope_; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_e_envelope_ AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'E (envelope protein)'::text) AND ((epi.protein_name)::text = 'E (envelope protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_e_envelope_ OWNER TO geco;

--
-- Name: epitope_2697049_m_membrane_; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_m_membrane_ AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'M (membrane glycoprotein)'::text) AND ((epi.protein_name)::text = 'M (membrane glycoprotein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_m_membrane_ OWNER TO geco;

--
-- Name: epitope_2697049_n_nucleocap; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_n_nucleocap AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'N (nucleocapsid phosphoprotein)'::text) AND ((epi.protein_name)::text = 'N (nucleocapsid phosphoprotein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_n_nucleocap OWNER TO geco;

--
-- Name: epitope_2697049_ns3_orf3a_p; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_ns3_orf3a_p AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NS3 (ORF3a protein)'::text) AND ((epi.protein_name)::text = 'NS3 (ORF3a protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_ns3_orf3a_p OWNER TO geco;

--
-- Name: epitope_2697049_ns6_orf6_pr; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_ns6_orf6_pr AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NS6 (ORF6 protein)'::text) AND ((epi.protein_name)::text = 'NS6 (ORF6 protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_ns6_orf6_pr OWNER TO geco;

--
-- Name: epitope_2697049_ns7a_orf7a_; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_ns7a_orf7a_ AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NS7a (ORF7a protein)'::text) AND ((epi.protein_name)::text = 'NS7a (ORF7a protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_ns7a_orf7a_ OWNER TO geco;

--
-- Name: epitope_2697049_ns7b_orf7b; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_ns7b_orf7b AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NS7b (ORF7b)'::text) AND ((epi.protein_name)::text = 'NS7b (ORF7b)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_ns7b_orf7b OWNER TO geco;

--
-- Name: epitope_2697049_ns8_orf8_pr; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_ns8_orf8_pr AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NS8 (ORF8 protein)'::text) AND ((epi.protein_name)::text = 'NS8 (ORF8 protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_ns8_orf8_pr OWNER TO geco;

--
-- Name: epitope_2697049_nsp10; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp10 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP10'::text) AND ((epi.protein_name)::text = 'NSP10'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp10 OWNER TO geco;

--
-- Name: epitope_2697049_nsp11; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp11 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP11'::text) AND ((epi.protein_name)::text = 'NSP11'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp11 OWNER TO geco;

--
-- Name: epitope_2697049_nsp12_rna_d; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp12_rna_d AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP12 (RNA-dependent RNA polymerase)'::text) AND ((epi.protein_name)::text = 'NSP12 (RNA-dependent RNA polymerase)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp12_rna_d OWNER TO geco;

--
-- Name: epitope_2697049_nsp13_helic; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp13_helic AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP13 (helicase)'::text) AND ((epi.protein_name)::text = 'NSP13 (helicase)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp13_helic OWNER TO geco;

--
-- Name: epitope_2697049_nsp14_3_to_; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp14_3_to_ AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP14 (3''-to-5'' exonuclease)'::text) AND ((epi.protein_name)::text = 'NSP14 (3''-to-5'' exonuclease)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp14_3_to_ OWNER TO geco;

--
-- Name: epitope_2697049_nsp15_endor; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp15_endor AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP15 (endoRNAse)'::text) AND ((epi.protein_name)::text = 'NSP15 (endoRNAse)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp15_endor OWNER TO geco;

--
-- Name: epitope_2697049_nsp16_2_o_r; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp16_2_o_r AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP16 (2''-O-ribose methyltransferase)'::text) AND ((epi.protein_name)::text = 'NSP16 (2''-O-ribose methyltransferase)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp16_2_o_r OWNER TO geco;

--
-- Name: epitope_2697049_nsp1_leader; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp1_leader AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP1 (leader protein)'::text) AND ((epi.protein_name)::text = 'NSP1 (leader protein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp1_leader OWNER TO geco;

--
-- Name: epitope_2697049_nsp2; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp2 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP2'::text) AND ((epi.protein_name)::text = 'NSP2'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp2 OWNER TO geco;

--
-- Name: epitope_2697049_nsp3; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp3 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP3'::text) AND ((epi.protein_name)::text = 'NSP3'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp3 OWNER TO geco;

--
-- Name: epitope_2697049_nsp4; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp4 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP4'::text) AND ((epi.protein_name)::text = 'NSP4'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp4 OWNER TO geco;

--
-- Name: epitope_2697049_nsp5_3c_lik; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp5_3c_lik AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP5 (3C-like proteinase)'::text) AND ((epi.protein_name)::text = 'NSP5 (3C-like proteinase)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp5_3c_lik OWNER TO geco;

--
-- Name: epitope_2697049_nsp6; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp6 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP6'::text) AND ((epi.protein_name)::text = 'NSP6'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp6 OWNER TO geco;

--
-- Name: epitope_2697049_nsp7; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp7 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP7'::text) AND ((epi.protein_name)::text = 'NSP7'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp7 OWNER TO geco;

--
-- Name: epitope_2697049_nsp8; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp8 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP8'::text) AND ((epi.protein_name)::text = 'NSP8'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp8 OWNER TO geco;

--
-- Name: epitope_2697049_nsp9; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_nsp9 AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'NSP9'::text) AND ((epi.protein_name)::text = 'NSP9'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_nsp9 OWNER TO geco;

--
-- Name: epitope_2697049_orf10_prote; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_orf10_prote AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'ORF10 protein'::text) AND ((epi.protein_name)::text = 'ORF10 protein'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_orf10_prote OWNER TO geco;

--
-- Name: epitope_2697049_orf1a_polyp; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_orf1a_polyp AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'ORF1a polyprotein'::text) AND ((epi.protein_name)::text = 'ORF1a polyprotein'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_orf1a_polyp OWNER TO geco;

--
-- Name: epitope_2697049_orf1ab_poly; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_orf1ab_poly AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'ORF1ab polyprotein'::text) AND ((epi.protein_name)::text = 'ORF1ab polyprotein'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_orf1ab_poly OWNER TO geco;

--
-- Name: epitope_2697049_spike_surfa; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.epitope_2697049_spike_surfa AS
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
   FROM (((((((public.epitope epi
     JOIN public.epitope_fragment epif ON ((epi.epitope_id = epif.epitope_id)))
     JOIN public.virus vir ON ((epi.virus_id = vir.virus_id)))
     JOIN public.host_specie hspec ON ((epi.host_id = hspec.host_id)))
     JOIN public.host_sample hsamp ON ((hspec.host_id = hsamp.host_id)))
     JOIN public.sequence seq ON (((hsamp.host_sample_id = seq.host_sample_id) AND (vir.virus_id = seq.virus_id))))
     JOIN public.annotation ann ON ((seq.sequence_id = ann.sequence_id)))
     JOIN public.aminoacid_variant amin ON ((ann.annotation_id = amin.annotation_id)))
  WHERE (((epi.protein_name)::text = (ann.product)::text) AND (amin.start_aa_original <= epif.epi_frag_annotation_stop) AND (amin.start_aa_original >= epif.epi_frag_annotation_start) AND ((ann.product)::text = 'Spike (surface glycoprotein)'::text) AND ((epi.protein_name)::text = 'Spike (surface glycoprotein)'::text) AND (vir.taxon_id = 2697049))
  ORDER BY epi.iedb_epitope_id
  WITH NO DATA;


ALTER TABLE public.epitope_2697049_spike_surfa OWNER TO geco;

--
-- Name: epitope_epitope_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.epitope_epitope_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.epitope_epitope_id_seq OWNER TO geco;

--
-- Name: epitope_epitope_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.epitope_epitope_id_seq OWNED BY public.epitope.epitope_id;


--
-- Name: epitope_fragment_epi_fragment_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.epitope_fragment_epi_fragment_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.epitope_fragment_epi_fragment_id_seq OWNER TO geco;

--
-- Name: epitope_fragment_epi_fragment_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.epitope_fragment_epi_fragment_id_seq OWNED BY public.epitope_fragment.epi_fragment_id;


--
-- Name: experiment_type; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.experiment_type (
    experiment_type_id integer NOT NULL,
    sequencing_technology character varying,
    assembly_method character varying,
    coverage integer
);


ALTER TABLE public.experiment_type OWNER TO geco;

--
-- Name: experiment_type_experiment_type_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.experiment_type_experiment_type_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.experiment_type_experiment_type_id_seq OWNER TO geco;

--
-- Name: experiment_type_experiment_type_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.experiment_type_experiment_type_id_seq OWNED BY public.experiment_type.experiment_type_id;


--
-- Name: host_sample_host_sample_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.host_sample_host_sample_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.host_sample_host_sample_id_seq OWNER TO geco;

--
-- Name: host_sample_host_sample_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.host_sample_host_sample_id_seq OWNED BY public.host_sample.host_sample_id;


--
-- Name: host_sample_view; Type: VIEW; Schema: public; Owner: geco
--

CREATE VIEW public.host_sample_view AS
 SELECT host_sample.host_sample_id,
    host_sample.host_id,
    host_sample.collection_date,
    host_sample.isolation_source,
    host_sample.originating_lab,
    host_sample.region,
    host_sample.country,
    host_sample.geo_group,
    host_sample.age,
    host_sample.gender,
    host_sample.province,
    host_sample.coll_date_precision,
    host_specie.host_taxon_id,
    host_specie.host_taxon_name
   FROM (public.host_sample
     JOIN public.host_specie USING (host_id));


ALTER TABLE public.host_sample_view OWNER TO geco;

--
-- Name: host_specie_host_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.host_specie_host_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.host_specie_host_id_seq OWNER TO geco;

--
-- Name: host_specie_host_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.host_specie_host_id_seq OWNED BY public.host_specie.host_id;


--
-- Name: location; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.location AS
 SELECT DISTINCT host_sample.province AS loc
   FROM public.host_sample
  WHERE ((host_sample.province IS NOT NULL) AND ((host_sample.province)::text <> ''::text))
UNION
 SELECT DISTINCT host_sample.region AS loc
   FROM public.host_sample
  WHERE ((host_sample.region IS NOT NULL) AND ((host_sample.region)::text <> ''::text))
UNION
 SELECT DISTINCT host_sample.country AS loc
   FROM public.host_sample
  WHERE ((host_sample.country IS NOT NULL) AND ((host_sample.country)::text <> ''::text))
UNION
 SELECT DISTINCT host_sample.geo_group AS loc
   FROM public.host_sample
  WHERE ((host_sample.geo_group IS NOT NULL) AND ((host_sample.geo_group)::text <> ''::text))
  WITH NO DATA;


ALTER TABLE public.location OWNER TO geco;

--
-- Name: nucleotide_sequence; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.nucleotide_sequence (
    sequence_id integer NOT NULL,
    nucleotide_sequence character varying
);


ALTER TABLE public.nucleotide_sequence OWNER TO geco;

--
-- Name: nucleotide_variant; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.nucleotide_variant (
    nucleotide_variant_id integer NOT NULL,
    sequence_id integer NOT NULL,
    sequence_original character varying NOT NULL,
    sequence_alternative character varying NOT NULL,
    start_original integer,
    start_alternative integer,
    variant_length integer NOT NULL,
    variant_type character varying NOT NULL
);


ALTER TABLE public.nucleotide_variant OWNER TO geco;

--
-- Name: nucleotide_variant_annotated; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.nucleotide_variant_annotated
WITH (fillfactor='100') AS
 SELECT nc.nucleotide_variant_id,
    nc.sequence_id,
    nc.variant_type,
    nc.start_original,
        CASE nc.variant_type
            WHEN 'SUB'::text THEN nc.sequence_original
            ELSE NULL::character varying
        END AS sequence_original,
        CASE nc.variant_type
            WHEN 'SUB'::text THEN nc.sequence_alternative
            ELSE NULL::character varying
        END AS sequence_alternative,
    nc.variant_length,
    ann.gene_name AS n_gene_name
   FROM (public.nucleotide_variant nc
     LEFT JOIN public.annotation ann ON (((nc.start_original >= ann.start) AND (nc.start_original <= ann.stop) AND (nc.sequence_id = ann.sequence_id) AND ((ann.feature_type)::text = 'gene'::text))))
  WITH NO DATA;


ALTER TABLE public.nucleotide_variant_annotated OWNER TO geco;

--
-- Name: nucleotide_variant_limited; Type: VIEW; Schema: public; Owner: geco
--

CREATE VIEW public.nucleotide_variant_limited AS
 SELECT nucleotide_variant.nucleotide_variant_id,
    nucleotide_variant.sequence_id,
    nucleotide_variant.sequence_original,
    nucleotide_variant.sequence_alternative,
    nucleotide_variant.start_original,
    nucleotide_variant.start_alternative,
    nucleotide_variant.variant_length,
    nucleotide_variant.variant_type
   FROM public.nucleotide_variant
  WHERE (nucleotide_variant.variant_length <= 20);


ALTER TABLE public.nucleotide_variant_limited OWNER TO geco;

--
-- Name: nucleotide_variant_nucleotide_variant_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.nucleotide_variant_nucleotide_variant_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.nucleotide_variant_nucleotide_variant_id_seq OWNER TO geco;

--
-- Name: nucleotide_variant_nucleotide_variant_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.nucleotide_variant_nucleotide_variant_id_seq OWNED BY public.nucleotide_variant.nucleotide_variant_id;


--
-- Name: overlap; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.overlap (
    sequence_id integer NOT NULL,
    accession_id character varying NOT NULL,
    overlapping_accession_id character varying NOT NULL,
    overlapping_source character varying NOT NULL
);


ALTER TABLE public.overlap OWNER TO geco;

--
-- Name: pipeline_event; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.pipeline_event (
    event_id integer NOT NULL,
    event_name character varying NOT NULL,
    event_date character varying,
    added_items integer,
    removed_items integer,
    changed_items integer
);


ALTER TABLE public.pipeline_event OWNER TO geco;

--
-- Name: pipeline_event_event_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.pipeline_event_event_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.pipeline_event_event_id_seq OWNER TO geco;

--
-- Name: pipeline_event_event_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.pipeline_event_event_id_seq OWNED BY public.pipeline_event.event_id;


--
-- Name: sequence_sequence_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.sequence_sequence_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.sequence_sequence_id_seq OWNER TO geco;

--
-- Name: sequence_sequence_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.sequence_sequence_id_seq OWNED BY public.sequence.sequence_id;


--
-- Name: sequencing_project; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.sequencing_project (
    sequencing_project_id integer NOT NULL,
    sequencing_lab character varying,
    submission_date date,
    database_source character varying,
    bioproject_id character varying
);


ALTER TABLE public.sequencing_project OWNER TO geco;

--
-- Name: sequencing_project_sequencing_project_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.sequencing_project_sequencing_project_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.sequencing_project_sequencing_project_id_seq OWNER TO geco;

--
-- Name: sequencing_project_sequencing_project_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.sequencing_project_sequencing_project_id_seq OWNED BY public.sequencing_project.sequencing_project_id;


--
-- Name: variant_impact; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.variant_impact (
    variant_impact_id integer NOT NULL,
    nucleotide_variant_id integer NOT NULL,
    effect character varying,
    putative_impact character varying,
    impact_gene_name character varying
);


ALTER TABLE public.variant_impact OWNER TO geco;

--
-- Name: variant_impact_variant_impact_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.variant_impact_variant_impact_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.variant_impact_variant_impact_id_seq OWNER TO geco;

--
-- Name: variant_impact_variant_impact_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.variant_impact_variant_impact_id_seq OWNED BY public.variant_impact.variant_impact_id;


--
-- Name: variant_observation; Type: MATERIALIZED VIEW; Schema: public; Owner: geco
--

CREATE MATERIALIZED VIEW public.variant_observation AS
 SELECT a.change,
    a.obs_location,
    a.collect_period,
    a.period_type,
    a.mutated_sequences,
    a.tot_sequences
   FROM ( WITH host_sample_periods AS (
                 SELECT host_sample.host_sample_id,
                    host_sample.host_id,
                    host_sample.collection_date,
                    host_sample.isolation_source,
                    host_sample.originating_lab,
                    host_sample.region,
                    host_sample.country,
                    host_sample.geo_group,
                    host_sample.age,
                    host_sample.gender,
                    host_sample.province,
                    host_sample.coll_date_precision,
                    "left"((host_sample.collection_date)::text, 7) AS month
                   FROM public.host_sample
                  WHERE (host_sample.coll_date_precision >= 1)
                ), change_location AS (
                 SELECT concat(split_part((a_1.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location
                   FROM ((((public.sequence s
                     JOIN public.host_sample h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a_1 USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a_1.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc
                ), time_period AS (
                 SELECT DISTINCT host_sample_periods.month AS collect_period
                   FROM host_sample_periods
                ), mutated AS (
                 SELECT concat(split_part((a_1.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location,
                    h.month AS collect_period,
                    count(*) AS mutated_sequences
                   FROM ((((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a_1 USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a_1.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc, h.month
                ), totals AS (
                 SELECT l.loc AS obs_location,
                    h.month AS collect_period,
                    count(*) AS tot_sequences
                   FROM ((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                  GROUP BY l.loc, h.month
                )
         SELECT change_location.change,
            change_location.obs_location,
            time_period.collect_period,
            'month'::text AS period_type,
            COALESCE(mutated.mutated_sequences, (0)::bigint) AS mutated_sequences,
            COALESCE(totals.tot_sequences, (0)::bigint) AS tot_sequences
           FROM (((change_location
             CROSS JOIN time_period)
             LEFT JOIN mutated ON (((change_location.change = mutated.change) AND ((change_location.obs_location)::text = (mutated.obs_location)::text) AND (time_period.collect_period = mutated.collect_period))))
             LEFT JOIN totals ON ((((change_location.obs_location)::text = (totals.obs_location)::text) AND (time_period.collect_period = totals.collect_period))))) a
UNION
 SELECT b.change,
    b.obs_location,
    b.collect_period,
    b.period_type,
    b.mutated_sequences,
    b.tot_sequences
   FROM ( WITH host_sample_periods AS (
                 SELECT host_sample.host_sample_id,
                    host_sample.host_id,
                    host_sample.collection_date,
                    host_sample.isolation_source,
                    host_sample.originating_lab,
                    host_sample.region,
                    host_sample.country,
                    host_sample.geo_group,
                    host_sample.age,
                    host_sample.gender,
                    host_sample.province,
                    host_sample.coll_date_precision,
                        CASE
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '01'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '15'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '01|15')
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '16'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '31'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '16|31')
                            ELSE NULL::text
                        END AS biweek
                   FROM public.host_sample
                  WHERE (host_sample.coll_date_precision >= 2)
                ), change_location AS (
                 SELECT concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location
                   FROM ((((public.sequence s
                     JOIN public.host_sample h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc
                ), time_period AS (
                 SELECT DISTINCT host_sample_periods.biweek AS collect_period
                   FROM host_sample_periods
                ), mutated AS (
                 SELECT concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location,
                    h.biweek AS collect_period,
                    count(*) AS mutated_sequences
                   FROM ((((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc, h.biweek
                ), totals AS (
                 SELECT l.loc AS obs_location,
                    h.biweek AS collect_period,
                    count(*) AS tot_sequences
                   FROM ((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                  GROUP BY l.loc, h.biweek
                )
         SELECT change_location.change,
            change_location.obs_location,
            time_period.collect_period,
            'biweek'::text AS period_type,
            COALESCE(mutated.mutated_sequences, (0)::bigint) AS mutated_sequences,
            COALESCE(totals.tot_sequences, (0)::bigint) AS tot_sequences
           FROM (((change_location
             CROSS JOIN time_period)
             LEFT JOIN mutated ON (((change_location.change = mutated.change) AND ((change_location.obs_location)::text = (mutated.obs_location)::text) AND (time_period.collect_period = mutated.collect_period))))
             LEFT JOIN totals ON ((((change_location.obs_location)::text = (totals.obs_location)::text) AND (time_period.collect_period = totals.collect_period))))) b
UNION
 SELECT c.change,
    c.obs_location,
    c.collect_period,
    c.period_type,
    c.mutated_sequences,
    c.tot_sequences
   FROM ( WITH host_sample_periods AS (
                 SELECT host_sample.host_sample_id,
                    host_sample.host_id,
                    host_sample.collection_date,
                    host_sample.isolation_source,
                    host_sample.originating_lab,
                    host_sample.region,
                    host_sample.country,
                    host_sample.geo_group,
                    host_sample.age,
                    host_sample.gender,
                    host_sample.province,
                    host_sample.coll_date_precision,
                        CASE
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '01'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '07'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '01|07')
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '08'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '15'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '08|15')
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '16'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '23'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '16|23')
                            WHEN (("right"((host_sample.collection_date)::text, 2) >= '24'::text) AND ("right"((host_sample.collection_date)::text, 2) <= '31'::text)) THEN concat("left"((host_sample.collection_date)::text, 8), '24|31')
                            ELSE NULL::text
                        END AS week
                   FROM public.host_sample
                  WHERE (host_sample.coll_date_precision >= 2)
                ), change_location AS (
                 SELECT concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location
                   FROM ((((public.sequence s
                     JOIN public.host_sample h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc
                ), time_period AS (
                 SELECT DISTINCT host_sample_periods.week AS collect_period
                   FROM host_sample_periods
                ), mutated AS (
                 SELECT concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative) AS change,
                    l.loc AS obs_location,
                    h.week AS collect_period,
                    count(*) AS mutated_sequences
                   FROM ((((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                     JOIN public.annotation a USING (sequence_id))
                     JOIN public.aminoacid_variant v USING (annotation_id))
                  GROUP BY (concat(split_part((a.product)::text, ' '::text, 1), '_', v.sequence_aa_original, v.start_aa_original, v.sequence_aa_alternative)), l.loc, h.week
                ), totals AS (
                 SELECT l.loc AS obs_location,
                    h.week AS collect_period,
                    count(*) AS tot_sequences
                   FROM ((public.sequence s
                     JOIN host_sample_periods h USING (host_sample_id))
                     JOIN public.location l ON ((((h.region)::text = (l.loc)::text) OR ((h.country)::text = (l.loc)::text) OR ((h.geo_group)::text = (l.loc)::text))))
                  GROUP BY l.loc, h.week
                )
         SELECT change_location.change,
            change_location.obs_location,
            time_period.collect_period,
            'week'::text AS period_type,
            COALESCE(mutated.mutated_sequences, (0)::bigint) AS mutated_sequences,
            COALESCE(totals.tot_sequences, (0)::bigint) AS tot_sequences
           FROM (((change_location
             CROSS JOIN time_period)
             LEFT JOIN mutated ON (((change_location.change = mutated.change) AND ((change_location.obs_location)::text = (mutated.obs_location)::text) AND (time_period.collect_period = mutated.collect_period))))
             LEFT JOIN totals ON ((((change_location.obs_location)::text = (totals.obs_location)::text) AND (time_period.collect_period = totals.collect_period))))) c
  WITH NO DATA;


ALTER TABLE public.variant_observation OWNER TO geco;

--
-- Name: virus_virus_id_seq; Type: SEQUENCE; Schema: public; Owner: geco
--

CREATE SEQUENCE public.virus_virus_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.virus_virus_id_seq OWNER TO geco;

--
-- Name: virus_virus_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: geco
--

ALTER SEQUENCE public.virus_virus_id_seq OWNED BY public.virus.virus_id;


--
-- Name: aminoacid_variant aminoacid_variant_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.aminoacid_variant ALTER COLUMN aminoacid_variant_id SET DEFAULT nextval('public.aminoacid_variant_aminoacid_variant_id_seq'::regclass);


--
-- Name: annotation annotation_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation ALTER COLUMN annotation_id SET DEFAULT nextval('public.annotation_annotation_id_seq'::regclass);


--
-- Name: epitope epitope_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope ALTER COLUMN epitope_id SET DEFAULT nextval('public.epitope_epitope_id_seq'::regclass);


--
-- Name: epitope_fragment epi_fragment_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope_fragment ALTER COLUMN epi_fragment_id SET DEFAULT nextval('public.epitope_fragment_epi_fragment_id_seq'::regclass);


--
-- Name: experiment_type experiment_type_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.experiment_type ALTER COLUMN experiment_type_id SET DEFAULT nextval('public.experiment_type_experiment_type_id_seq'::regclass);


--
-- Name: host_sample host_sample_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.host_sample ALTER COLUMN host_sample_id SET DEFAULT nextval('public.host_sample_host_sample_id_seq'::regclass);


--
-- Name: host_specie host_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.host_specie ALTER COLUMN host_id SET DEFAULT nextval('public.host_specie_host_id_seq'::regclass);


--
-- Name: nucleotide_variant nucleotide_variant_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.nucleotide_variant ALTER COLUMN nucleotide_variant_id SET DEFAULT nextval('public.nucleotide_variant_nucleotide_variant_id_seq'::regclass);


--
-- Name: pipeline_event event_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.pipeline_event ALTER COLUMN event_id SET DEFAULT nextval('public.pipeline_event_event_id_seq'::regclass);


--
-- Name: sequence sequence_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence ALTER COLUMN sequence_id SET DEFAULT nextval('public.sequence_sequence_id_seq'::regclass);


--
-- Name: sequencing_project sequencing_project_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequencing_project ALTER COLUMN sequencing_project_id SET DEFAULT nextval('public.sequencing_project_sequencing_project_id_seq'::regclass);


--
-- Name: variant_impact variant_impact_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.variant_impact ALTER COLUMN variant_impact_id SET DEFAULT nextval('public.variant_impact_variant_impact_id_seq'::regclass);


--
-- Name: virus virus_id; Type: DEFAULT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.virus ALTER COLUMN virus_id SET DEFAULT nextval('public.virus_virus_id_seq'::regclass);


--
-- Name: aminoacid_variant aminoacid_variant_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.aminoacid_variant
    ADD CONSTRAINT aminoacid_variant_pkey PRIMARY KEY (aminoacid_variant_id);


--
-- Name: annotation annotation_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation
    ADD CONSTRAINT annotation_pkey PRIMARY KEY (annotation_id);


--
-- Name: annotation_sequence annotation_sequence_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation_sequence
    ADD CONSTRAINT annotation_sequence_pkey PRIMARY KEY (annotation_id);


--
-- Name: db_meta db_meta_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.db_meta
    ADD CONSTRAINT db_meta_pkey PRIMARY KEY (virus_id, source);


--
-- Name: epitope_fragment epitope_fragment_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope_fragment
    ADD CONSTRAINT epitope_fragment_pkey PRIMARY KEY (epi_fragment_id);


--
-- Name: epitope epitope_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope
    ADD CONSTRAINT epitope_pkey PRIMARY KEY (epitope_id);


--
-- Name: experiment_type experiment_type_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.experiment_type
    ADD CONSTRAINT experiment_type_pkey PRIMARY KEY (experiment_type_id);


--
-- Name: host_sample host_sample_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.host_sample
    ADD CONSTRAINT host_sample_pkey PRIMARY KEY (host_sample_id);


--
-- Name: host_specie host_specie_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.host_specie
    ADD CONSTRAINT host_specie_pkey PRIMARY KEY (host_id);


--
-- Name: nucleotide_sequence nucleotide_sequence_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.nucleotide_sequence
    ADD CONSTRAINT nucleotide_sequence_pkey PRIMARY KEY (sequence_id);


--
-- Name: nucleotide_variant nucleotide_variant_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.nucleotide_variant
    ADD CONSTRAINT nucleotide_variant_pkey PRIMARY KEY (nucleotide_variant_id);


--
-- Name: overlap overlap_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.overlap
    ADD CONSTRAINT overlap_pkey PRIMARY KEY (sequence_id, overlapping_accession_id);


--
-- Name: pipeline_event pipeline_event_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.pipeline_event
    ADD CONSTRAINT pipeline_event_pkey PRIMARY KEY (event_id);


--
-- Name: sequence sequence_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence
    ADD CONSTRAINT sequence_pkey PRIMARY KEY (sequence_id);


--
-- Name: sequencing_project sequencing_project_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequencing_project
    ADD CONSTRAINT sequencing_project_pkey PRIMARY KEY (sequencing_project_id);


--
-- Name: variant_impact variant_impact_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.variant_impact
    ADD CONSTRAINT variant_impact_pkey PRIMARY KEY (variant_impact_id);


--
-- Name: virus virus_pkey; Type: CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.virus
    ADD CONSTRAINT virus_pkey PRIMARY KEY (virus_id);


--
-- Name: aa__ann_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX aa__ann_id ON public.aminoacid_variant USING btree (annotation_id);


--
-- Name: aa__start_original; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX aa__start_original ON public.aminoacid_variant USING btree (start_aa_original);


--
-- Name: aa__var_type_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX aa__var_type_lower ON public.aminoacid_variant USING btree (lower((variant_aa_type)::text));


--
-- Name: aa__var_type_normal; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX aa__var_type_normal ON public.aminoacid_variant USING btree (variant_aa_type);


--
-- Name: ann__seq_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX ann__seq_id ON public.annotation USING btree (sequence_id);


--
-- Name: ann__start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX ann__start ON public.annotation USING btree (start);


--
-- Name: ann__stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX ann__stop ON public.annotation USING btree (stop);


--
-- Name: ann_seq__seq_id__product; Type: INDEX; Schema: public; Owner: geco
--

CREATE UNIQUE INDEX ann_seq__seq_id__product ON public.annotation_sequence USING btree (sequence_id, product);


--
-- Name: epi_2697049_e_envelope___cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___cell_type__idx ON public.epitope_2697049_e_envelope_ USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___epi_annotation_start__idx ON public.epitope_2697049_e_envelope_ USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___epi_annotation_stop__idx ON public.epitope_2697049_e_envelope_ USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___epi_frag_annotation_start__i ON public.epitope_2697049_e_envelope_ USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___epi_frag_annotation_stop__id ON public.epitope_2697049_e_envelope_ USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___host_taxon_id__idx ON public.epitope_2697049_e_envelope_ USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___host_taxon_name_lower__idx ON public.epitope_2697049_e_envelope_ USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___iedb_epitope_id__idx ON public.epitope_2697049_e_envelope_ USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___is_linear__idx ON public.epitope_2697049_e_envelope_ USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___mhc_allele__idx ON public.epitope_2697049_e_envelope_ USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___mhc_class_lower__idx ON public.epitope_2697049_e_envelope_ USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___product_lower__idx ON public.epitope_2697049_e_envelope_ USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___response_frequency_pos__idx ON public.epitope_2697049_e_envelope_ USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___sequence_aa_alternative__idx ON public.epitope_2697049_e_envelope_ USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___sequence_aa_original__idx ON public.epitope_2697049_e_envelope_ USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___start_aa_original__idx ON public.epitope_2697049_e_envelope_ USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___taxon_id__idx ON public.epitope_2697049_e_envelope_ USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___taxon_name_lower__idx ON public.epitope_2697049_e_envelope_ USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___variant_aa_length__idx ON public.epitope_2697049_e_envelope_ USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___variant_aa_type__idx ON public.epitope_2697049_e_envelope_ USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_cell_type__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_epi_start__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_epi_stop__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_is_linear__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_mhc_allele__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_product__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_host_resp_freq__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope___virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope___virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_e_envelope_ USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___cell_type__idx ON public.epitope_2697049_m_membrane_ USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___epi_annotation_start__idx ON public.epitope_2697049_m_membrane_ USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___epi_annotation_stop__idx ON public.epitope_2697049_m_membrane_ USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___epi_frag_annotation_start__i ON public.epitope_2697049_m_membrane_ USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___epi_frag_annotation_stop__id ON public.epitope_2697049_m_membrane_ USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___host_taxon_id__idx ON public.epitope_2697049_m_membrane_ USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___host_taxon_name_lower__idx ON public.epitope_2697049_m_membrane_ USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___iedb_epitope_id__idx ON public.epitope_2697049_m_membrane_ USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___is_linear__idx ON public.epitope_2697049_m_membrane_ USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___mhc_allele__idx ON public.epitope_2697049_m_membrane_ USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___mhc_class_lower__idx ON public.epitope_2697049_m_membrane_ USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___product_lower__idx ON public.epitope_2697049_m_membrane_ USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___response_frequency_pos__idx ON public.epitope_2697049_m_membrane_ USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___sequence_aa_alternative__idx ON public.epitope_2697049_m_membrane_ USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___sequence_aa_original__idx ON public.epitope_2697049_m_membrane_ USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___start_aa_original__idx ON public.epitope_2697049_m_membrane_ USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___taxon_id__idx ON public.epitope_2697049_m_membrane_ USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___taxon_name_lower__idx ON public.epitope_2697049_m_membrane_ USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___variant_aa_length__idx ON public.epitope_2697049_m_membrane_ USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___variant_aa_type__idx ON public.epitope_2697049_m_membrane_ USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_cell_type__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_epi_start__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_epi_stop__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_is_linear__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_mhc_allele__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_product__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_host_resp_freq__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane___virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane___virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_m_membrane_ USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__cell_type__idx ON public.epitope_2697049_n_nucleocap USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__epi_annotation_start__idx ON public.epitope_2697049_n_nucleocap USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__epi_annotation_stop__idx ON public.epitope_2697049_n_nucleocap USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__epi_frag_annotation_start__i ON public.epitope_2697049_n_nucleocap USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__epi_frag_annotation_stop__id ON public.epitope_2697049_n_nucleocap USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__host_taxon_id__idx ON public.epitope_2697049_n_nucleocap USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__host_taxon_name_lower__idx ON public.epitope_2697049_n_nucleocap USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__iedb_epitope_id__idx ON public.epitope_2697049_n_nucleocap USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__is_linear__idx ON public.epitope_2697049_n_nucleocap USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__mhc_allele__idx ON public.epitope_2697049_n_nucleocap USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__mhc_class_lower__idx ON public.epitope_2697049_n_nucleocap USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__product_lower__idx ON public.epitope_2697049_n_nucleocap USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__response_frequency_pos__idx ON public.epitope_2697049_n_nucleocap USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__sequence_aa_alternative__idx ON public.epitope_2697049_n_nucleocap USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__sequence_aa_original__idx ON public.epitope_2697049_n_nucleocap USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__start_aa_original__idx ON public.epitope_2697049_n_nucleocap USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__taxon_id__idx ON public.epitope_2697049_n_nucleocap USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__taxon_name_lower__idx ON public.epitope_2697049_n_nucleocap USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__variant_aa_length__idx ON public.epitope_2697049_n_nucleocap USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__variant_aa_type__idx ON public.epitope_2697049_n_nucleocap USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_cell_type__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_epi_start__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_epi_stop__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_is_linear__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_mhc_allele__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_product__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_host_resp_freq__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocap__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocap__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_n_nucleocap USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__cell_type__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__epi_annotation_start__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__epi_annotation_stop__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__epi_frag_annotation_start__i ON public.epitope_2697049_ns3_orf3a_p USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__epi_frag_annotation_stop__id ON public.epitope_2697049_ns3_orf3a_p USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__host_taxon_id__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__host_taxon_name_lower__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__iedb_epitope_id__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__is_linear__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__mhc_allele__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__mhc_class_lower__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__product_lower__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__response_frequency_pos__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__sequence_aa_alternative__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__sequence_aa_original__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__start_aa_original__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__taxon_id__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__taxon_name_lower__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__variant_aa_length__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__variant_aa_type__idx ON public.epitope_2697049_ns3_orf3a_p USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_cell_type__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_epi_start__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_epi_stop__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_is_linear__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_mhc_allele__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_product__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_host_resp_freq__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_p__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_p__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_ns3_orf3a_p USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__cell_type__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__epi_annotation_start__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__epi_annotation_stop__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__epi_frag_annotation_start__i ON public.epitope_2697049_ns6_orf6_pr USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__epi_frag_annotation_stop__id ON public.epitope_2697049_ns6_orf6_pr USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__host_taxon_id__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__host_taxon_name_lower__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__iedb_epitope_id__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__is_linear__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__mhc_allele__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__mhc_class_lower__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__product_lower__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__response_frequency_pos__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__sequence_aa_alternative__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__sequence_aa_original__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__start_aa_original__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__taxon_id__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__taxon_name_lower__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__variant_aa_length__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__variant_aa_type__idx ON public.epitope_2697049_ns6_orf6_pr USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_cell_type__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_epi_start__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_epi_stop__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_is_linear__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_mhc_allele__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_product__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_host_resp_freq__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_pr__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_pr__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_ns6_orf6_pr USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___cell_type__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___epi_annotation_start__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___epi_annotation_stop__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___epi_frag_annotation_start__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___epi_frag_annotation_stop__id ON public.epitope_2697049_ns7a_orf7a_ USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___host_taxon_id__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___host_taxon_name_lower__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___iedb_epitope_id__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___is_linear__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___mhc_allele__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___mhc_class_lower__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___product_lower__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___response_frequency_pos__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___sequence_aa_alternative__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___sequence_aa_original__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___start_aa_original__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___taxon_id__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___taxon_name_lower__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___variant_aa_length__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___variant_aa_type__idx ON public.epitope_2697049_ns7a_orf7a_ USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_cell_type__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_epi_start__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_epi_stop__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_is_linear__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_mhc_allele__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_product__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_host_resp_freq__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a___virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a___virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_ns7a_orf7a_ USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__cell_type__idx ON public.epitope_2697049_ns7b_orf7b USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_annotation_start__idx ON public.epitope_2697049_ns7b_orf7b USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_annotation_stop__idx ON public.epitope_2697049_ns7b_orf7b USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_annotation_start__i ON public.epitope_2697049_ns7b_orf7b USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_annotation_stop__id ON public.epitope_2697049_ns7b_orf7b USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__host_taxon_id__idx ON public.epitope_2697049_ns7b_orf7b USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__host_taxon_name_lower__idx ON public.epitope_2697049_ns7b_orf7b USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__iedb_epitope_id__idx ON public.epitope_2697049_ns7b_orf7b USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__is_linear__idx ON public.epitope_2697049_ns7b_orf7b USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__mhc_allele__idx ON public.epitope_2697049_ns7b_orf7b USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__mhc_class_lower__idx ON public.epitope_2697049_ns7b_orf7b USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__product_lower__idx ON public.epitope_2697049_ns7b_orf7b USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__response_frequency_pos__idx ON public.epitope_2697049_ns7b_orf7b USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__sequence_aa_alternative__idx ON public.epitope_2697049_ns7b_orf7b USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__sequence_aa_original__idx ON public.epitope_2697049_ns7b_orf7b USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__start_aa_original__idx ON public.epitope_2697049_ns7b_orf7b USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__taxon_id__idx ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__taxon_name_lower__idx ON public.epitope_2697049_ns7b_orf7b USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_length__idx ON public.epitope_2697049_ns7b_orf7b USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_type__idx ON public.epitope_2697049_ns7b_orf7b USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_cell_type__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_start__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_stop__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_is_linear__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_mhc_allele__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_product__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_resp_freq__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__cell_type__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__epi_annotation_start__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__epi_annotation_stop__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__epi_frag_annotation_start__i ON public.epitope_2697049_ns8_orf8_pr USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__epi_frag_annotation_stop__id ON public.epitope_2697049_ns8_orf8_pr USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__host_taxon_id__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__host_taxon_name_lower__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__iedb_epitope_id__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__is_linear__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__mhc_allele__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__mhc_class_lower__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__product_lower__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__response_frequency_pos__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__sequence_aa_alternative__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__sequence_aa_original__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__start_aa_original__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__taxon_id__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__taxon_name_lower__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__variant_aa_length__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__variant_aa_type__idx ON public.epitope_2697049_ns8_orf8_pr USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_cell_type__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_epi_start__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_epi_stop__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_is_linear__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_mhc_allele__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_product__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_host_resp_freq__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_pr__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_pr__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_ns8_orf8_pr USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__cell_type__idx ON public.epitope_2697049_nsp10 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_annotation_start__idx ON public.epitope_2697049_nsp10 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_annotation_stop__idx ON public.epitope_2697049_nsp10 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_frag_annotation_start__i ON public.epitope_2697049_nsp10 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp10 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__host_taxon_id__idx ON public.epitope_2697049_nsp10 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__host_taxon_name_lower__idx ON public.epitope_2697049_nsp10 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__iedb_epitope_id__idx ON public.epitope_2697049_nsp10 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__is_linear__idx ON public.epitope_2697049_nsp10 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__mhc_allele__idx ON public.epitope_2697049_nsp10 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__mhc_class_lower__idx ON public.epitope_2697049_nsp10 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__product_lower__idx ON public.epitope_2697049_nsp10 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__response_frequency_pos__idx ON public.epitope_2697049_nsp10 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__sequence_aa_alternative__idx ON public.epitope_2697049_nsp10 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__sequence_aa_original__idx ON public.epitope_2697049_nsp10 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__start_aa_original__idx ON public.epitope_2697049_nsp10 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__taxon_id__idx ON public.epitope_2697049_nsp10 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__taxon_name_lower__idx ON public.epitope_2697049_nsp10 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__variant_aa_length__idx ON public.epitope_2697049_nsp10 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__variant_aa_type__idx ON public.epitope_2697049_nsp10 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_cell_type__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_epi_start__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_epi_stop__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_is_linear__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_mhc_allele__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_product__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_resp_freq__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__cell_type__idx ON public.epitope_2697049_nsp11 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_annotation_start__idx ON public.epitope_2697049_nsp11 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_annotation_stop__idx ON public.epitope_2697049_nsp11 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_frag_annotation_start__i ON public.epitope_2697049_nsp11 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp11 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__host_taxon_id__idx ON public.epitope_2697049_nsp11 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__host_taxon_name_lower__idx ON public.epitope_2697049_nsp11 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__iedb_epitope_id__idx ON public.epitope_2697049_nsp11 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__is_linear__idx ON public.epitope_2697049_nsp11 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__mhc_allele__idx ON public.epitope_2697049_nsp11 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__mhc_class_lower__idx ON public.epitope_2697049_nsp11 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__product_lower__idx ON public.epitope_2697049_nsp11 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__response_frequency_pos__idx ON public.epitope_2697049_nsp11 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__sequence_aa_alternative__idx ON public.epitope_2697049_nsp11 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__sequence_aa_original__idx ON public.epitope_2697049_nsp11 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__start_aa_original__idx ON public.epitope_2697049_nsp11 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__taxon_id__idx ON public.epitope_2697049_nsp11 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__taxon_name_lower__idx ON public.epitope_2697049_nsp11 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__variant_aa_length__idx ON public.epitope_2697049_nsp11 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__variant_aa_type__idx ON public.epitope_2697049_nsp11 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_cell_type__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_epi_start__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_epi_stop__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_is_linear__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_mhc_allele__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_product__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_resp_freq__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__cell_type__idx ON public.epitope_2697049_nsp12_rna_d USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__epi_annotation_start__idx ON public.epitope_2697049_nsp12_rna_d USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__epi_annotation_stop__idx ON public.epitope_2697049_nsp12_rna_d USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__epi_frag_annotation_start__i ON public.epitope_2697049_nsp12_rna_d USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp12_rna_d USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__host_taxon_id__idx ON public.epitope_2697049_nsp12_rna_d USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__host_taxon_name_lower__idx ON public.epitope_2697049_nsp12_rna_d USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__iedb_epitope_id__idx ON public.epitope_2697049_nsp12_rna_d USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__is_linear__idx ON public.epitope_2697049_nsp12_rna_d USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__mhc_allele__idx ON public.epitope_2697049_nsp12_rna_d USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__mhc_class_lower__idx ON public.epitope_2697049_nsp12_rna_d USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__product_lower__idx ON public.epitope_2697049_nsp12_rna_d USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__response_frequency_pos__idx ON public.epitope_2697049_nsp12_rna_d USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__sequence_aa_alternative__idx ON public.epitope_2697049_nsp12_rna_d USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__sequence_aa_original__idx ON public.epitope_2697049_nsp12_rna_d USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__start_aa_original__idx ON public.epitope_2697049_nsp12_rna_d USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__taxon_id__idx ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__taxon_name_lower__idx ON public.epitope_2697049_nsp12_rna_d USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__variant_aa_length__idx ON public.epitope_2697049_nsp12_rna_d USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__variant_aa_type__idx ON public.epitope_2697049_nsp12_rna_d USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_cell_type__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_epi_start__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_epi_stop__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_is_linear__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_mhc_allele__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_product__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_host_resp_freq__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_d__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_d__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp12_rna_d USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__cell_type__idx ON public.epitope_2697049_nsp13_helic USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__epi_annotation_start__idx ON public.epitope_2697049_nsp13_helic USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__epi_annotation_stop__idx ON public.epitope_2697049_nsp13_helic USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__epi_frag_annotation_start__i ON public.epitope_2697049_nsp13_helic USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp13_helic USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__host_taxon_id__idx ON public.epitope_2697049_nsp13_helic USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__host_taxon_name_lower__idx ON public.epitope_2697049_nsp13_helic USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__iedb_epitope_id__idx ON public.epitope_2697049_nsp13_helic USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__is_linear__idx ON public.epitope_2697049_nsp13_helic USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__mhc_allele__idx ON public.epitope_2697049_nsp13_helic USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__mhc_class_lower__idx ON public.epitope_2697049_nsp13_helic USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__product_lower__idx ON public.epitope_2697049_nsp13_helic USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__response_frequency_pos__idx ON public.epitope_2697049_nsp13_helic USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__sequence_aa_alternative__idx ON public.epitope_2697049_nsp13_helic USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__sequence_aa_original__idx ON public.epitope_2697049_nsp13_helic USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__start_aa_original__idx ON public.epitope_2697049_nsp13_helic USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__taxon_id__idx ON public.epitope_2697049_nsp13_helic USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__taxon_name_lower__idx ON public.epitope_2697049_nsp13_helic USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__variant_aa_length__idx ON public.epitope_2697049_nsp13_helic USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__variant_aa_type__idx ON public.epitope_2697049_nsp13_helic USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_cell_type__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_epi_start__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_epi_stop__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_is_linear__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_mhc_allele__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_product__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_host_resp_freq__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helic__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helic__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp13_helic USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___cell_type__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___epi_annotation_start__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___epi_annotation_stop__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___epi_frag_annotation_start__i ON public.epitope_2697049_nsp14_3_to_ USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___epi_frag_annotation_stop__id ON public.epitope_2697049_nsp14_3_to_ USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___host_taxon_id__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___host_taxon_name_lower__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___iedb_epitope_id__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___is_linear__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___mhc_allele__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___mhc_class_lower__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___product_lower__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___response_frequency_pos__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___sequence_aa_alternative__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___sequence_aa_original__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___start_aa_original__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___taxon_id__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___taxon_name_lower__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___variant_aa_length__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___variant_aa_type__idx ON public.epitope_2697049_nsp14_3_to_ USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_cell_type__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_epi_start__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_epi_stop__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_is_linear__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_mhc_allele__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_product__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_host_resp_freq__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to___virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to___virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp14_3_to_ USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__cell_type__idx ON public.epitope_2697049_nsp15_endor USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__epi_annotation_start__idx ON public.epitope_2697049_nsp15_endor USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__epi_annotation_stop__idx ON public.epitope_2697049_nsp15_endor USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__epi_frag_annotation_start__i ON public.epitope_2697049_nsp15_endor USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp15_endor USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__host_taxon_id__idx ON public.epitope_2697049_nsp15_endor USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__host_taxon_name_lower__idx ON public.epitope_2697049_nsp15_endor USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__iedb_epitope_id__idx ON public.epitope_2697049_nsp15_endor USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__is_linear__idx ON public.epitope_2697049_nsp15_endor USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__mhc_allele__idx ON public.epitope_2697049_nsp15_endor USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__mhc_class_lower__idx ON public.epitope_2697049_nsp15_endor USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__product_lower__idx ON public.epitope_2697049_nsp15_endor USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__response_frequency_pos__idx ON public.epitope_2697049_nsp15_endor USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__sequence_aa_alternative__idx ON public.epitope_2697049_nsp15_endor USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__sequence_aa_original__idx ON public.epitope_2697049_nsp15_endor USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__start_aa_original__idx ON public.epitope_2697049_nsp15_endor USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__taxon_id__idx ON public.epitope_2697049_nsp15_endor USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__taxon_name_lower__idx ON public.epitope_2697049_nsp15_endor USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__variant_aa_length__idx ON public.epitope_2697049_nsp15_endor USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__variant_aa_type__idx ON public.epitope_2697049_nsp15_endor USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_cell_type__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_epi_start__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_epi_stop__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_is_linear__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_mhc_allele__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_product__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_host_resp_freq__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endor__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endor__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp15_endor USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__cell_type__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__epi_annotation_start__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__epi_annotation_stop__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__epi_frag_annotation_start__i ON public.epitope_2697049_nsp16_2_o_r USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp16_2_o_r USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__host_taxon_id__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__host_taxon_name_lower__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__iedb_epitope_id__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__is_linear__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__mhc_allele__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__mhc_class_lower__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__product_lower__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__response_frequency_pos__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__sequence_aa_alternative__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__sequence_aa_original__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__start_aa_original__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__taxon_id__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__taxon_name_lower__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__variant_aa_length__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__variant_aa_type__idx ON public.epitope_2697049_nsp16_2_o_r USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_cell_type__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_epi_start__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_epi_stop__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_is_linear__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_mhc_allele__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_product__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_host_resp_freq__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_r__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_r__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp16_2_o_r USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__cell_type__idx ON public.epitope_2697049_nsp1_leader USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__epi_annotation_start__idx ON public.epitope_2697049_nsp1_leader USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__epi_annotation_stop__idx ON public.epitope_2697049_nsp1_leader USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__epi_frag_annotation_start__i ON public.epitope_2697049_nsp1_leader USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp1_leader USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__host_taxon_id__idx ON public.epitope_2697049_nsp1_leader USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__host_taxon_name_lower__idx ON public.epitope_2697049_nsp1_leader USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__iedb_epitope_id__idx ON public.epitope_2697049_nsp1_leader USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__is_linear__idx ON public.epitope_2697049_nsp1_leader USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__mhc_allele__idx ON public.epitope_2697049_nsp1_leader USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__mhc_class_lower__idx ON public.epitope_2697049_nsp1_leader USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__product_lower__idx ON public.epitope_2697049_nsp1_leader USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__response_frequency_pos__idx ON public.epitope_2697049_nsp1_leader USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__sequence_aa_alternative__idx ON public.epitope_2697049_nsp1_leader USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__sequence_aa_original__idx ON public.epitope_2697049_nsp1_leader USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__start_aa_original__idx ON public.epitope_2697049_nsp1_leader USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__taxon_id__idx ON public.epitope_2697049_nsp1_leader USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__taxon_name_lower__idx ON public.epitope_2697049_nsp1_leader USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__variant_aa_length__idx ON public.epitope_2697049_nsp1_leader USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__variant_aa_type__idx ON public.epitope_2697049_nsp1_leader USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_cell_type__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_epi_start__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_epi_stop__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_is_linear__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_mhc_allele__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_product__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_host_resp_freq__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp1_leader USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__cell_type__idx ON public.epitope_2697049_nsp2 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_annotation_start__idx ON public.epitope_2697049_nsp2 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_annotation_stop__idx ON public.epitope_2697049_nsp2 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_frag_annotation_start__i ON public.epitope_2697049_nsp2 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp2 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__host_taxon_id__idx ON public.epitope_2697049_nsp2 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__host_taxon_name_lower__idx ON public.epitope_2697049_nsp2 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__iedb_epitope_id__idx ON public.epitope_2697049_nsp2 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__is_linear__idx ON public.epitope_2697049_nsp2 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__mhc_allele__idx ON public.epitope_2697049_nsp2 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__mhc_class_lower__idx ON public.epitope_2697049_nsp2 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__product_lower__idx ON public.epitope_2697049_nsp2 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__response_frequency_pos__idx ON public.epitope_2697049_nsp2 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__sequence_aa_alternative__idx ON public.epitope_2697049_nsp2 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__sequence_aa_original__idx ON public.epitope_2697049_nsp2 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__start_aa_original__idx ON public.epitope_2697049_nsp2 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__taxon_id__idx ON public.epitope_2697049_nsp2 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__taxon_name_lower__idx ON public.epitope_2697049_nsp2 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__variant_aa_length__idx ON public.epitope_2697049_nsp2 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__variant_aa_type__idx ON public.epitope_2697049_nsp2 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_cell_type__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_epi_start__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_epi_stop__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_is_linear__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_mhc_allele__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_product__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_resp_freq__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__cell_type__idx ON public.epitope_2697049_nsp3 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_annotation_start__idx ON public.epitope_2697049_nsp3 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_annotation_stop__idx ON public.epitope_2697049_nsp3 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_frag_annotation_start__i ON public.epitope_2697049_nsp3 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp3 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__host_taxon_id__idx ON public.epitope_2697049_nsp3 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__host_taxon_name_lower__idx ON public.epitope_2697049_nsp3 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__iedb_epitope_id__idx ON public.epitope_2697049_nsp3 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__is_linear__idx ON public.epitope_2697049_nsp3 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__mhc_allele__idx ON public.epitope_2697049_nsp3 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__mhc_class_lower__idx ON public.epitope_2697049_nsp3 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__product_lower__idx ON public.epitope_2697049_nsp3 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__response_frequency_pos__idx ON public.epitope_2697049_nsp3 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__sequence_aa_alternative__idx ON public.epitope_2697049_nsp3 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__sequence_aa_original__idx ON public.epitope_2697049_nsp3 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__start_aa_original__idx ON public.epitope_2697049_nsp3 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__taxon_id__idx ON public.epitope_2697049_nsp3 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__taxon_name_lower__idx ON public.epitope_2697049_nsp3 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__variant_aa_length__idx ON public.epitope_2697049_nsp3 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__variant_aa_type__idx ON public.epitope_2697049_nsp3 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_cell_type__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_epi_start__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_epi_stop__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_is_linear__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_mhc_allele__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_product__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_resp_freq__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__cell_type__idx ON public.epitope_2697049_nsp4 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_annotation_start__idx ON public.epitope_2697049_nsp4 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_annotation_stop__idx ON public.epitope_2697049_nsp4 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_frag_annotation_start__i ON public.epitope_2697049_nsp4 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp4 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__host_taxon_id__idx ON public.epitope_2697049_nsp4 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__host_taxon_name_lower__idx ON public.epitope_2697049_nsp4 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__iedb_epitope_id__idx ON public.epitope_2697049_nsp4 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__is_linear__idx ON public.epitope_2697049_nsp4 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__mhc_allele__idx ON public.epitope_2697049_nsp4 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__mhc_class_lower__idx ON public.epitope_2697049_nsp4 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__product_lower__idx ON public.epitope_2697049_nsp4 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__response_frequency_pos__idx ON public.epitope_2697049_nsp4 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__sequence_aa_alternative__idx ON public.epitope_2697049_nsp4 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__sequence_aa_original__idx ON public.epitope_2697049_nsp4 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__start_aa_original__idx ON public.epitope_2697049_nsp4 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__taxon_id__idx ON public.epitope_2697049_nsp4 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__taxon_name_lower__idx ON public.epitope_2697049_nsp4 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__variant_aa_length__idx ON public.epitope_2697049_nsp4 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__variant_aa_type__idx ON public.epitope_2697049_nsp4 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_cell_type__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_epi_start__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_epi_stop__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_is_linear__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_mhc_allele__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_product__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_resp_freq__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__cell_type__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__epi_annotation_start__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__epi_annotation_stop__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__epi_frag_annotation_start__i ON public.epitope_2697049_nsp5_3c_lik USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp5_3c_lik USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__host_taxon_id__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__host_taxon_name_lower__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__iedb_epitope_id__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__is_linear__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__mhc_allele__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__mhc_class_lower__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__product_lower__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__response_frequency_pos__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__sequence_aa_alternative__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__sequence_aa_original__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__start_aa_original__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__taxon_id__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__taxon_name_lower__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__variant_aa_length__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__variant_aa_type__idx ON public.epitope_2697049_nsp5_3c_lik USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_cell_type__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_epi_start__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_epi_stop__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_is_linear__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_mhc_allele__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_product__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_host_resp_freq__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_lik__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_lik__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp5_3c_lik USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__cell_type__idx ON public.epitope_2697049_nsp6 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_annotation_start__idx ON public.epitope_2697049_nsp6 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_annotation_stop__idx ON public.epitope_2697049_nsp6 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_frag_annotation_start__i ON public.epitope_2697049_nsp6 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp6 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__host_taxon_id__idx ON public.epitope_2697049_nsp6 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__host_taxon_name_lower__idx ON public.epitope_2697049_nsp6 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__iedb_epitope_id__idx ON public.epitope_2697049_nsp6 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__is_linear__idx ON public.epitope_2697049_nsp6 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__mhc_allele__idx ON public.epitope_2697049_nsp6 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__mhc_class_lower__idx ON public.epitope_2697049_nsp6 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__product_lower__idx ON public.epitope_2697049_nsp6 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__response_frequency_pos__idx ON public.epitope_2697049_nsp6 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__sequence_aa_alternative__idx ON public.epitope_2697049_nsp6 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__sequence_aa_original__idx ON public.epitope_2697049_nsp6 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__start_aa_original__idx ON public.epitope_2697049_nsp6 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__taxon_id__idx ON public.epitope_2697049_nsp6 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__taxon_name_lower__idx ON public.epitope_2697049_nsp6 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__variant_aa_length__idx ON public.epitope_2697049_nsp6 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__variant_aa_type__idx ON public.epitope_2697049_nsp6 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_cell_type__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_epi_start__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_epi_stop__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_is_linear__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_mhc_allele__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_product__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_resp_freq__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__cell_type__idx ON public.epitope_2697049_nsp7 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_annotation_start__idx ON public.epitope_2697049_nsp7 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_annotation_stop__idx ON public.epitope_2697049_nsp7 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_frag_annotation_start__i ON public.epitope_2697049_nsp7 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp7 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__host_taxon_id__idx ON public.epitope_2697049_nsp7 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__host_taxon_name_lower__idx ON public.epitope_2697049_nsp7 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__iedb_epitope_id__idx ON public.epitope_2697049_nsp7 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__is_linear__idx ON public.epitope_2697049_nsp7 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__mhc_allele__idx ON public.epitope_2697049_nsp7 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__mhc_class_lower__idx ON public.epitope_2697049_nsp7 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__product_lower__idx ON public.epitope_2697049_nsp7 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__response_frequency_pos__idx ON public.epitope_2697049_nsp7 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__sequence_aa_alternative__idx ON public.epitope_2697049_nsp7 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__sequence_aa_original__idx ON public.epitope_2697049_nsp7 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__start_aa_original__idx ON public.epitope_2697049_nsp7 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__taxon_id__idx ON public.epitope_2697049_nsp7 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__taxon_name_lower__idx ON public.epitope_2697049_nsp7 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__variant_aa_length__idx ON public.epitope_2697049_nsp7 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__variant_aa_type__idx ON public.epitope_2697049_nsp7 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_cell_type__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_epi_start__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_epi_stop__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_is_linear__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_mhc_allele__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_product__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_resp_freq__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__cell_type__idx ON public.epitope_2697049_nsp8 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_annotation_start__idx ON public.epitope_2697049_nsp8 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_annotation_stop__idx ON public.epitope_2697049_nsp8 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_frag_annotation_start__i ON public.epitope_2697049_nsp8 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp8 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__host_taxon_id__idx ON public.epitope_2697049_nsp8 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__host_taxon_name_lower__idx ON public.epitope_2697049_nsp8 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__iedb_epitope_id__idx ON public.epitope_2697049_nsp8 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__is_linear__idx ON public.epitope_2697049_nsp8 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__mhc_allele__idx ON public.epitope_2697049_nsp8 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__mhc_class_lower__idx ON public.epitope_2697049_nsp8 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__product_lower__idx ON public.epitope_2697049_nsp8 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__response_frequency_pos__idx ON public.epitope_2697049_nsp8 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__sequence_aa_alternative__idx ON public.epitope_2697049_nsp8 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__sequence_aa_original__idx ON public.epitope_2697049_nsp8 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__start_aa_original__idx ON public.epitope_2697049_nsp8 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__taxon_id__idx ON public.epitope_2697049_nsp8 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__taxon_name_lower__idx ON public.epitope_2697049_nsp8 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__variant_aa_length__idx ON public.epitope_2697049_nsp8 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__variant_aa_type__idx ON public.epitope_2697049_nsp8 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_cell_type__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_epi_start__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_epi_stop__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_is_linear__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_mhc_allele__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_product__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_resp_freq__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__cell_type__idx ON public.epitope_2697049_nsp9 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_annotation_start__idx ON public.epitope_2697049_nsp9 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_annotation_stop__idx ON public.epitope_2697049_nsp9 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_frag_annotation_start__i ON public.epitope_2697049_nsp9 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_frag_annotation_stop__id ON public.epitope_2697049_nsp9 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__host_taxon_id__idx ON public.epitope_2697049_nsp9 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__host_taxon_name_lower__idx ON public.epitope_2697049_nsp9 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__iedb_epitope_id__idx ON public.epitope_2697049_nsp9 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__is_linear__idx ON public.epitope_2697049_nsp9 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__mhc_allele__idx ON public.epitope_2697049_nsp9 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__mhc_class_lower__idx ON public.epitope_2697049_nsp9 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__product_lower__idx ON public.epitope_2697049_nsp9 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__response_frequency_pos__idx ON public.epitope_2697049_nsp9 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__sequence_aa_alternative__idx ON public.epitope_2697049_nsp9 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__sequence_aa_original__idx ON public.epitope_2697049_nsp9 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__start_aa_original__idx ON public.epitope_2697049_nsp9 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__taxon_id__idx ON public.epitope_2697049_nsp9 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__taxon_name_lower__idx ON public.epitope_2697049_nsp9 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__variant_aa_length__idx ON public.epitope_2697049_nsp9 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__variant_aa_type__idx ON public.epitope_2697049_nsp9 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_cell_type__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_epi_start__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_epi_stop__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_is_linear__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_mhc_allele__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_product__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_resp_freq__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__cell_type__idx ON public.epitope_2697049_orf10_prote USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__epi_annotation_start__idx ON public.epitope_2697049_orf10_prote USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__epi_annotation_stop__idx ON public.epitope_2697049_orf10_prote USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__epi_frag_annotation_start__i ON public.epitope_2697049_orf10_prote USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__epi_frag_annotation_stop__id ON public.epitope_2697049_orf10_prote USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__host_taxon_id__idx ON public.epitope_2697049_orf10_prote USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__host_taxon_name_lower__idx ON public.epitope_2697049_orf10_prote USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__iedb_epitope_id__idx ON public.epitope_2697049_orf10_prote USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__is_linear__idx ON public.epitope_2697049_orf10_prote USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__mhc_allele__idx ON public.epitope_2697049_orf10_prote USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__mhc_class_lower__idx ON public.epitope_2697049_orf10_prote USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__product_lower__idx ON public.epitope_2697049_orf10_prote USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__response_frequency_pos__idx ON public.epitope_2697049_orf10_prote USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__sequence_aa_alternative__idx ON public.epitope_2697049_orf10_prote USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__sequence_aa_original__idx ON public.epitope_2697049_orf10_prote USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__start_aa_original__idx ON public.epitope_2697049_orf10_prote USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__taxon_id__idx ON public.epitope_2697049_orf10_prote USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__taxon_name_lower__idx ON public.epitope_2697049_orf10_prote USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__variant_aa_length__idx ON public.epitope_2697049_orf10_prote USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__variant_aa_type__idx ON public.epitope_2697049_orf10_prote USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_cell_type__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_epi_start__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_epi_stop__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_is_linear__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_mhc_allele__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_product__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_host_resp_freq__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_prote__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_prote__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_orf10_prote USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__cell_type__idx ON public.epitope_2697049_orf1a_polyp USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__epi_annotation_start__idx ON public.epitope_2697049_orf1a_polyp USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__epi_annotation_stop__idx ON public.epitope_2697049_orf1a_polyp USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__epi_frag_annotation_start__i ON public.epitope_2697049_orf1a_polyp USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__epi_frag_annotation_stop__id ON public.epitope_2697049_orf1a_polyp USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__host_taxon_id__idx ON public.epitope_2697049_orf1a_polyp USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__host_taxon_name_lower__idx ON public.epitope_2697049_orf1a_polyp USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__iedb_epitope_id__idx ON public.epitope_2697049_orf1a_polyp USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__is_linear__idx ON public.epitope_2697049_orf1a_polyp USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__mhc_allele__idx ON public.epitope_2697049_orf1a_polyp USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__mhc_class_lower__idx ON public.epitope_2697049_orf1a_polyp USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__product_lower__idx ON public.epitope_2697049_orf1a_polyp USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__response_frequency_pos__idx ON public.epitope_2697049_orf1a_polyp USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__sequence_aa_alternative__idx ON public.epitope_2697049_orf1a_polyp USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__sequence_aa_original__idx ON public.epitope_2697049_orf1a_polyp USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__start_aa_original__idx ON public.epitope_2697049_orf1a_polyp USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__taxon_id__idx ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__taxon_name_lower__idx ON public.epitope_2697049_orf1a_polyp USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__variant_aa_length__idx ON public.epitope_2697049_orf1a_polyp USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__variant_aa_type__idx ON public.epitope_2697049_orf1a_polyp USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_cell_type__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_epi_start__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_epi_stop__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_is_linear__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_mhc_allele__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_product__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_host_resp_freq__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyp__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyp__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_orf1a_polyp USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__cell_type__idx ON public.epitope_2697049_orf1ab_poly USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__epi_annotation_start__idx ON public.epitope_2697049_orf1ab_poly USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__epi_annotation_stop__idx ON public.epitope_2697049_orf1ab_poly USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__epi_frag_annotation_start__i ON public.epitope_2697049_orf1ab_poly USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__epi_frag_annotation_stop__id ON public.epitope_2697049_orf1ab_poly USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__host_taxon_id__idx ON public.epitope_2697049_orf1ab_poly USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__host_taxon_name_lower__idx ON public.epitope_2697049_orf1ab_poly USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__iedb_epitope_id__idx ON public.epitope_2697049_orf1ab_poly USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__is_linear__idx ON public.epitope_2697049_orf1ab_poly USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__mhc_allele__idx ON public.epitope_2697049_orf1ab_poly USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__mhc_class_lower__idx ON public.epitope_2697049_orf1ab_poly USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__product_lower__idx ON public.epitope_2697049_orf1ab_poly USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__response_frequency_pos__idx ON public.epitope_2697049_orf1ab_poly USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__sequence_aa_alternative__idx ON public.epitope_2697049_orf1ab_poly USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__sequence_aa_original__idx ON public.epitope_2697049_orf1ab_poly USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__start_aa_original__idx ON public.epitope_2697049_orf1ab_poly USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__taxon_id__idx ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__taxon_name_lower__idx ON public.epitope_2697049_orf1ab_poly USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__variant_aa_length__idx ON public.epitope_2697049_orf1ab_poly USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__variant_aa_type__idx ON public.epitope_2697049_orf1ab_poly USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_cell_type__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_epi_start__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_epi_stop__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_is_linear__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_mhc_allele__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_product__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_host_resp_freq__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_poly__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_poly__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_orf1ab_poly USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__cell_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__cell_type__idx ON public.epitope_2697049_spike_surfa USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__epi_annotation_start__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__epi_annotation_start__idx ON public.epitope_2697049_spike_surfa USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__epi_annotation_stop__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__epi_annotation_stop__idx ON public.epitope_2697049_spike_surfa USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__epi_frag_annotation_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__epi_frag_annotation_start__i ON public.epitope_2697049_spike_surfa USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__epi_frag_annotation_stop__id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__epi_frag_annotation_stop__id ON public.epitope_2697049_spike_surfa USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__host_taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__host_taxon_id__idx ON public.epitope_2697049_spike_surfa USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__host_taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__host_taxon_name_lower__idx ON public.epitope_2697049_spike_surfa USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__iedb_epitope_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__iedb_epitope_id__idx ON public.epitope_2697049_spike_surfa USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__is_linear__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__is_linear__idx ON public.epitope_2697049_spike_surfa USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__mhc_allele__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__mhc_allele__idx ON public.epitope_2697049_spike_surfa USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__mhc_class_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__mhc_class_lower__idx ON public.epitope_2697049_spike_surfa USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__product_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__product_lower__idx ON public.epitope_2697049_spike_surfa USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__response_frequency_pos__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__response_frequency_pos__idx ON public.epitope_2697049_spike_surfa USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__sequence_aa_alternative__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__sequence_aa_alternative__idx ON public.epitope_2697049_spike_surfa USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__sequence_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__sequence_aa_original__idx ON public.epitope_2697049_spike_surfa USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__start_aa_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__start_aa_original__idx ON public.epitope_2697049_spike_surfa USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__taxon_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__taxon_id__idx ON public.epitope_2697049_spike_surfa USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__taxon_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__taxon_name_lower__idx ON public.epitope_2697049_spike_surfa USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__variant_aa_length__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__variant_aa_length__idx ON public.epitope_2697049_spike_surfa USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__variant_aa_type__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__variant_aa_type__idx ON public.epitope_2697049_spike_surfa USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_cell_type__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_cell_type__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_epi_start__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_epi_start__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_epi_stop__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_epi_stop__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_is_linear__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_is_linear__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_mhc_allele__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_mhc_allele__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_product__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_product__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_host_resp_freq__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_host_resp_freq__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surfa__virus_taxon_and_host_taxon_id__i; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surfa__virus_taxon_and_host_taxon_id__i ON public.epitope_2697049_spike_surfa USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: idx_loc; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX idx_loc ON public.location USING gin (loc);


--
-- Name: impact__var_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX impact__var_id ON public.variant_impact USING btree (nucleotide_variant_id);


--
-- Name: nuc_var__length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nuc_var__length ON public.nucleotide_variant USING btree (variant_length);


--
-- Name: nuc_var__seq_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nuc_var__seq_id ON public.nucleotide_variant USING btree (sequence_id);


--
-- Name: nuc_var__start_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nuc_var__start_alt ON public.nucleotide_variant USING btree (start_alternative);


--
-- Name: nuc_var__start_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nuc_var__start_orig ON public.nucleotide_variant USING btree (start_original);


--
-- Name: nucleotide_variant_annotated__n_gene_name_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__n_gene_name_lower__idx ON public.nucleotide_variant_annotated USING btree (lower((n_gene_name)::text)) WITH (fillfactor='100');


--
-- Name: nucleotide_variant_annotated__nucleotide_variant_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__nucleotide_variant_id__idx ON public.nucleotide_variant_annotated USING btree (nucleotide_variant_id) WITH (fillfactor='100');


--
-- Name: nucleotide_variant_annotated__sequence_alternative_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__sequence_alternative_lower__idx ON public.nucleotide_variant_annotated USING btree (lower((sequence_alternative)::text)) WITH (fillfactor='100');


--
-- Name: nucleotide_variant_annotated__sequence_id__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__sequence_id__idx ON public.nucleotide_variant_annotated USING btree (sequence_id) WITH (fillfactor='100');

ALTER TABLE public.nucleotide_variant_annotated CLUSTER ON nucleotide_variant_annotated__sequence_id__idx;


--
-- Name: nucleotide_variant_annotated__sequence_original_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__sequence_original_lower__idx ON public.nucleotide_variant_annotated USING btree (lower((sequence_original)::text)) WITH (fillfactor='100');


--
-- Name: nucleotide_variant_annotated__start_original__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__start_original__idx ON public.nucleotide_variant_annotated USING btree (start_original) WITH (fillfactor='100');


--
-- Name: nucleotide_variant_annotated__variant_type_lower__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX nucleotide_variant_annotated__variant_type_lower__idx ON public.nucleotide_variant_annotated USING btree (lower((variant_type)::text)) WITH (fillfactor='100');


--
-- Name: seq__accession_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE UNIQUE INDEX seq__accession_id ON public.sequence USING btree (lower((accession_id)::text));


--
-- Name: seq__alternative_accession_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE UNIQUE INDEX seq__alternative_accession_id ON public.sequence USING btree (lower((alternative_accession_id)::text));


--
-- Name: seq__experiment_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX seq__experiment_id ON public.sequence USING btree (experiment_type_id);


--
-- Name: seq__host_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX seq__host_id ON public.sequence USING btree (host_sample_id);


--
-- Name: seq__seq_proj_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX seq__seq_proj_id ON public.sequence USING btree (sequencing_project_id);


--
-- Name: seq__virus_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX seq__virus_id ON public.sequence USING btree (virus_id);


--
-- Name: sequence__is_reference__idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX sequence__is_reference__idx ON public.sequence USING btree (is_reference);


--
-- Name: strain_name_gist_idx; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX strain_name_gist_idx ON public.sequence USING gist (lower((strain_name)::text) public.gist_trgm_ops);


--
-- Name: aminoacid_variant aminoacid_variant_annotation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.aminoacid_variant
    ADD CONSTRAINT aminoacid_variant_annotation_id_fkey FOREIGN KEY (annotation_id) REFERENCES public.annotation(annotation_id);


--
-- Name: annotation_sequence annotation_sequence_annotation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation_sequence
    ADD CONSTRAINT annotation_sequence_annotation_id_fkey FOREIGN KEY (annotation_id) REFERENCES public.annotation(annotation_id);


--
-- Name: annotation annotation_sequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation
    ADD CONSTRAINT annotation_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES public.sequence(sequence_id);


--
-- Name: annotation_sequence annotation_sequence_sequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.annotation_sequence
    ADD CONSTRAINT annotation_sequence_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES public.sequence(sequence_id);


--
-- Name: epitope_fragment epitope_fragment_epitope_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope_fragment
    ADD CONSTRAINT epitope_fragment_epitope_id_fkey FOREIGN KEY (epitope_id) REFERENCES public.epitope(epitope_id);


--
-- Name: epitope epitope_host_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope
    ADD CONSTRAINT epitope_host_id_fkey FOREIGN KEY (host_id) REFERENCES public.host_specie(host_id);


--
-- Name: epitope epitope_virus_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.epitope
    ADD CONSTRAINT epitope_virus_id_fkey FOREIGN KEY (virus_id) REFERENCES public.virus(virus_id);


--
-- Name: host_sample host_sample_host_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.host_sample
    ADD CONSTRAINT host_sample_host_id_fkey FOREIGN KEY (host_id) REFERENCES public.host_specie(host_id);


--
-- Name: nucleotide_sequence nucleotide_sequence_sequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.nucleotide_sequence
    ADD CONSTRAINT nucleotide_sequence_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES public.sequence(sequence_id);


--
-- Name: nucleotide_variant nucleotide_variant_sequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.nucleotide_variant
    ADD CONSTRAINT nucleotide_variant_sequence_id_fkey FOREIGN KEY (sequence_id) REFERENCES public.sequence(sequence_id);


--
-- Name: sequence sequence_experiment_type_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence
    ADD CONSTRAINT sequence_experiment_type_id_fkey FOREIGN KEY (experiment_type_id) REFERENCES public.experiment_type(experiment_type_id);


--
-- Name: sequence sequence_host_sample_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence
    ADD CONSTRAINT sequence_host_sample_id_fkey FOREIGN KEY (host_sample_id) REFERENCES public.host_sample(host_sample_id);


--
-- Name: sequence sequence_sequencing_project_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence
    ADD CONSTRAINT sequence_sequencing_project_id_fkey FOREIGN KEY (sequencing_project_id) REFERENCES public.sequencing_project(sequencing_project_id);


--
-- Name: sequence sequence_virus_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.sequence
    ADD CONSTRAINT sequence_virus_id_fkey FOREIGN KEY (virus_id) REFERENCES public.virus(virus_id);


--
-- Name: variant_impact variant_impact_nucleotide_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: geco
--

ALTER TABLE ONLY public.variant_impact
    ADD CONSTRAINT variant_impact_nucleotide_variant_id_fkey FOREIGN KEY (nucleotide_variant_id) REFERENCES public.nucleotide_variant(nucleotide_variant_id);


--
-- PostgreSQL database dump complete
--

