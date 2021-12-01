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
-- Name: epitope_2697049_e_envelope_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_e_envelope_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_e_envelope_protein OWNER TO geco;

--
-- Name: epitope_2697049_m_membrane_glycoprotein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_m_membrane_glycoprotein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_m_membrane_glycoprotein OWNER TO geco;

--
-- Name: epitope_2697049_n_nucleocapsid_phosphoprotei; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_n_nucleocapsid_phosphoprotei (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_n_nucleocapsid_phosphoprotei OWNER TO geco;

--
-- Name: epitope_2697049_ns3_orf3a_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_ns3_orf3a_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_ns3_orf3a_protein OWNER TO geco;

--
-- Name: epitope_2697049_ns6_orf6_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_ns6_orf6_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_ns6_orf6_protein OWNER TO geco;

--
-- Name: epitope_2697049_ns7a_orf7a_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_ns7a_orf7a_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_ns7a_orf7a_protein OWNER TO geco;

--
-- Name: epitope_2697049_ns7b_orf7b; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_ns7b_orf7b (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_ns7b_orf7b OWNER TO geco;

--
-- Name: epitope_2697049_ns8_orf8_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_ns8_orf8_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_ns8_orf8_protein OWNER TO geco;

--
-- Name: epitope_2697049_nsp10; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp10 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp10 OWNER TO geco;

--
-- Name: epitope_2697049_nsp11; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp11 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp11 OWNER TO geco;

--
-- Name: epitope_2697049_nsp12_rna_dependent_rna_poly; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp12_rna_dependent_rna_poly (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp12_rna_dependent_rna_poly OWNER TO geco;

--
-- Name: epitope_2697049_nsp13_helicase; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp13_helicase (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp13_helicase OWNER TO geco;

--
-- Name: epitope_2697049_nsp14_3_to_5_exonuclease; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp14_3_to_5_exonuclease (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp14_3_to_5_exonuclease OWNER TO geco;

--
-- Name: epitope_2697049_nsp15_endornase; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp15_endornase (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp15_endornase OWNER TO geco;

--
-- Name: epitope_2697049_nsp16_2_o_ribose_methyltrans; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp16_2_o_ribose_methyltrans (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp16_2_o_ribose_methyltrans OWNER TO geco;

--
-- Name: epitope_2697049_nsp1_leader_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp1_leader_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp1_leader_protein OWNER TO geco;

--
-- Name: epitope_2697049_nsp2; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp2 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp2 OWNER TO geco;

--
-- Name: epitope_2697049_nsp3; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp3 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp3 OWNER TO geco;

--
-- Name: epitope_2697049_nsp4; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp4 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp4 OWNER TO geco;

--
-- Name: epitope_2697049_nsp5_3c_like_proteinase; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp5_3c_like_proteinase (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp5_3c_like_proteinase OWNER TO geco;

--
-- Name: epitope_2697049_nsp6; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp6 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp6 OWNER TO geco;

--
-- Name: epitope_2697049_nsp7; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp7 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp7 OWNER TO geco;

--
-- Name: epitope_2697049_nsp8; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp8 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp8 OWNER TO geco;

--
-- Name: epitope_2697049_nsp9; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_nsp9 (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_nsp9 OWNER TO geco;

--
-- Name: epitope_2697049_orf10_protein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_orf10_protein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_orf10_protein OWNER TO geco;

--
-- Name: epitope_2697049_orf1a_polyprotein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_orf1a_polyprotein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_orf1a_polyprotein OWNER TO geco;

--
-- Name: epitope_2697049_orf1ab_polyprotein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_orf1ab_polyprotein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_orf1ab_polyprotein OWNER TO geco;

--
-- Name: epitope_2697049_spike_surface_glycoprotein; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.epitope_2697049_spike_surface_glycoprotein (
    iedb_epitope_id integer,
    epitope_iri character varying,
    cell_type character varying,
    mhc_class character varying,
    mhc_allele character varying,
    response_frequency_pos real,
    epi_annotation_start integer,
    epi_annotation_stop integer,
    is_linear boolean,
    assay_type character varying,
    epi_fragment_sequence character varying,
    epi_frag_annotation_start integer,
    epi_frag_annotation_stop integer,
    taxon_id integer,
    taxon_name character varying,
    host_taxon_id integer,
    host_taxon_name character varying,
    sequence_id integer,
    product character varying,
    aminoacid_variant_id integer,
    start_aa_original integer,
    sequence_aa_original character varying NOT NULL,
    sequence_aa_alternative character varying NOT NULL,
    variant_aa_length integer NOT NULL,
    variant_aa_type character varying NOT NULL
);


ALTER TABLE public.epitope_2697049_spike_surface_glycoprotein OWNER TO geco;

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
-- Name: host_specie; Type: TABLE; Schema: public; Owner: geco
--

CREATE TABLE public.host_specie (
    host_id integer NOT NULL,
    host_taxon_id integer,
    host_taxon_name character varying
);


ALTER TABLE public.host_specie OWNER TO geco;

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
-- Name: epi_2697049_e_envelope_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__cell_type ON public.epitope_2697049_e_envelope_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__epi_an_nstop ON public.epitope_2697049_e_envelope_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__epi_an_start ON public.epitope_2697049_e_envelope_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__epi_frag_an_start ON public.epitope_2697049_e_envelope_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__epi_frag_an_stop ON public.epitope_2697049_e_envelope_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__host_tax_id ON public.epitope_2697049_e_envelope_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__host_tax_name ON public.epitope_2697049_e_envelope_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__iedb_id ON public.epitope_2697049_e_envelope_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__is_linear ON public.epitope_2697049_e_envelope_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__mhc_allele ON public.epitope_2697049_e_envelope_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__mhc_class_lower ON public.epitope_2697049_e_envelope_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__product_lower ON public.epitope_2697049_e_envelope_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__response_freq_pos ON public.epitope_2697049_e_envelope_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__seq_aa_alt ON public.epitope_2697049_e_envelope_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__seq_aa_orig ON public.epitope_2697049_e_envelope_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__start_aa_orig ON public.epitope_2697049_e_envelope_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__taxon_id ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__taxon_name_lower ON public.epitope_2697049_e_envelope_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__variant_aa_length ON public.epitope_2697049_e_envelope_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__variant_aa_type ON public.epitope_2697049_e_envelope_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_host_cell_type ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_host_epi_start ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_host_mhc_allele ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_host_product ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_host_resp_freq ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__vir_n_host_tax_id ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__virus_host_epi_stop ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_e_envelope_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_e_envelope_protein__virus_host_is_linear ON public.epitope_2697049_e_envelope_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__cell_type ON public.epitope_2697049_m_membrane_glycoprotein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_an_nstop ON public.epitope_2697049_m_membrane_glycoprotein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_an_start ON public.epitope_2697049_m_membrane_glycoprotein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_frag_an_start ON public.epitope_2697049_m_membrane_glycoprotein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__epi_frag_an_stop ON public.epitope_2697049_m_membrane_glycoprotein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__host_tax_id ON public.epitope_2697049_m_membrane_glycoprotein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__host_tax_name ON public.epitope_2697049_m_membrane_glycoprotein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__iedb_id ON public.epitope_2697049_m_membrane_glycoprotein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__is_linear ON public.epitope_2697049_m_membrane_glycoprotein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__mhc_allele ON public.epitope_2697049_m_membrane_glycoprotein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__mhc_class_lower ON public.epitope_2697049_m_membrane_glycoprotein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__product_lower ON public.epitope_2697049_m_membrane_glycoprotein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__response_freq_pos ON public.epitope_2697049_m_membrane_glycoprotein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__seq_aa_alt ON public.epitope_2697049_m_membrane_glycoprotein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__seq_aa_orig ON public.epitope_2697049_m_membrane_glycoprotein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__start_aa_orig ON public.epitope_2697049_m_membrane_glycoprotein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__taxon_id ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__taxon_name_lower ON public.epitope_2697049_m_membrane_glycoprotein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__variant_aa_length ON public.epitope_2697049_m_membrane_glycoprotein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__variant_aa_type ON public.epitope_2697049_m_membrane_glycoprotein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_cell_type ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_epi_start ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_mhc_allele ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_product ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_host_resp_freq ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__vir_n_host_tax_id ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__virus_host_epi_stop ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_m_membrane_glycoprotein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_m_membrane_glycoprotein__virus_host_is_linear ON public.epitope_2697049_m_membrane_glycoprotein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__cell_type ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_nstop ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_an_start ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_start ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__epi_frag_an_stop ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_id ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__host_tax_name ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__iedb_id ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__is_linear ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__mhc_allele ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__mhc_class_lower ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__product_lower ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__response_freq_pos ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_alt ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__seq_aa_orig ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__start_aa_orig ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__taxon_id ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__taxon_name_lower ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_length ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__variant_aa_type ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_cell_type ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_epi_start ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_mhc_allele ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_product ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_host_resp_freq ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__vir_n_host_tax_id ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_epi_stop ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_n_nucleocapsid_phosphoprotei__virus_host_is_linear ON public.epitope_2697049_n_nucleocapsid_phosphoprotei USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__cell_type ON public.epitope_2697049_ns3_orf3a_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_an_nstop ON public.epitope_2697049_ns3_orf3a_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_an_start ON public.epitope_2697049_ns3_orf3a_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_frag_an_start ON public.epitope_2697049_ns3_orf3a_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__epi_frag_an_stop ON public.epitope_2697049_ns3_orf3a_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__host_tax_id ON public.epitope_2697049_ns3_orf3a_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__host_tax_name ON public.epitope_2697049_ns3_orf3a_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__iedb_id ON public.epitope_2697049_ns3_orf3a_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__is_linear ON public.epitope_2697049_ns3_orf3a_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__mhc_allele ON public.epitope_2697049_ns3_orf3a_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__mhc_class_lower ON public.epitope_2697049_ns3_orf3a_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__product_lower ON public.epitope_2697049_ns3_orf3a_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__response_freq_pos ON public.epitope_2697049_ns3_orf3a_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__seq_aa_alt ON public.epitope_2697049_ns3_orf3a_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__seq_aa_orig ON public.epitope_2697049_ns3_orf3a_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__start_aa_orig ON public.epitope_2697049_ns3_orf3a_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__taxon_id ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__taxon_name_lower ON public.epitope_2697049_ns3_orf3a_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__variant_aa_length ON public.epitope_2697049_ns3_orf3a_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__variant_aa_type ON public.epitope_2697049_ns3_orf3a_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_cell_type ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_epi_start ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_mhc_allele ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_product ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_host_resp_freq ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__vir_n_host_tax_id ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__virus_host_epi_stop ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns3_orf3a_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns3_orf3a_protein__virus_host_is_linear ON public.epitope_2697049_ns3_orf3a_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__cell_type ON public.epitope_2697049_ns6_orf6_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__epi_an_nstop ON public.epitope_2697049_ns6_orf6_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__epi_an_start ON public.epitope_2697049_ns6_orf6_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__epi_frag_an_start ON public.epitope_2697049_ns6_orf6_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__epi_frag_an_stop ON public.epitope_2697049_ns6_orf6_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__host_tax_id ON public.epitope_2697049_ns6_orf6_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__host_tax_name ON public.epitope_2697049_ns6_orf6_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__iedb_id ON public.epitope_2697049_ns6_orf6_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__is_linear ON public.epitope_2697049_ns6_orf6_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__mhc_allele ON public.epitope_2697049_ns6_orf6_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__mhc_class_lower ON public.epitope_2697049_ns6_orf6_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__product_lower ON public.epitope_2697049_ns6_orf6_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__response_freq_pos ON public.epitope_2697049_ns6_orf6_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__seq_aa_alt ON public.epitope_2697049_ns6_orf6_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__seq_aa_orig ON public.epitope_2697049_ns6_orf6_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__start_aa_orig ON public.epitope_2697049_ns6_orf6_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__taxon_id ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__taxon_name_lower ON public.epitope_2697049_ns6_orf6_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__variant_aa_length ON public.epitope_2697049_ns6_orf6_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__variant_aa_type ON public.epitope_2697049_ns6_orf6_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_cell_type ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_epi_start ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_mhc_allele ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_product ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_host_resp_freq ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__vir_n_host_tax_id ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__virus_host_epi_stop ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns6_orf6_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns6_orf6_protein__virus_host_is_linear ON public.epitope_2697049_ns6_orf6_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__cell_type ON public.epitope_2697049_ns7a_orf7a_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_an_nstop ON public.epitope_2697049_ns7a_orf7a_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_an_start ON public.epitope_2697049_ns7a_orf7a_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_frag_an_start ON public.epitope_2697049_ns7a_orf7a_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__epi_frag_an_stop ON public.epitope_2697049_ns7a_orf7a_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__host_tax_id ON public.epitope_2697049_ns7a_orf7a_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__host_tax_name ON public.epitope_2697049_ns7a_orf7a_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__iedb_id ON public.epitope_2697049_ns7a_orf7a_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__is_linear ON public.epitope_2697049_ns7a_orf7a_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__mhc_allele ON public.epitope_2697049_ns7a_orf7a_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__mhc_class_lower ON public.epitope_2697049_ns7a_orf7a_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__product_lower ON public.epitope_2697049_ns7a_orf7a_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__response_freq_pos ON public.epitope_2697049_ns7a_orf7a_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__seq_aa_alt ON public.epitope_2697049_ns7a_orf7a_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__seq_aa_orig ON public.epitope_2697049_ns7a_orf7a_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__start_aa_orig ON public.epitope_2697049_ns7a_orf7a_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__taxon_id ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__taxon_name_lower ON public.epitope_2697049_ns7a_orf7a_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__variant_aa_length ON public.epitope_2697049_ns7a_orf7a_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__variant_aa_type ON public.epitope_2697049_ns7a_orf7a_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_cell_type ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_epi_start ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_mhc_allele ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_product ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_host_resp_freq ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__vir_n_host_tax_id ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__virus_host_epi_stop ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7a_orf7a_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7a_orf7a_protein__virus_host_is_linear ON public.epitope_2697049_ns7a_orf7a_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__cell_type ON public.epitope_2697049_ns7b_orf7b USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_an_nstop ON public.epitope_2697049_ns7b_orf7b USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_an_start ON public.epitope_2697049_ns7b_orf7b USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_an_start ON public.epitope_2697049_ns7b_orf7b USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__epi_frag_an_stop ON public.epitope_2697049_ns7b_orf7b USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__host_tax_id ON public.epitope_2697049_ns7b_orf7b USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__host_tax_name ON public.epitope_2697049_ns7b_orf7b USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__iedb_id ON public.epitope_2697049_ns7b_orf7b USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__is_linear ON public.epitope_2697049_ns7b_orf7b USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__mhc_allele ON public.epitope_2697049_ns7b_orf7b USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__mhc_class_lower ON public.epitope_2697049_ns7b_orf7b USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__product_lower ON public.epitope_2697049_ns7b_orf7b USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__response_freq_pos ON public.epitope_2697049_ns7b_orf7b USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__seq_aa_alt ON public.epitope_2697049_ns7b_orf7b USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__seq_aa_orig ON public.epitope_2697049_ns7b_orf7b USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__start_aa_orig ON public.epitope_2697049_ns7b_orf7b USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__taxon_id ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__taxon_name_lower ON public.epitope_2697049_ns7b_orf7b USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_length ON public.epitope_2697049_ns7b_orf7b USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__variant_aa_type ON public.epitope_2697049_ns7b_orf7b USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_cell_type ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_epi_start ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_mhc_allele ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_product ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_host_resp_freq ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__vir_n_host_tax_id ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_epi_stop ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns7b_orf7b__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns7b_orf7b__virus_host_is_linear ON public.epitope_2697049_ns7b_orf7b USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__cell_type ON public.epitope_2697049_ns8_orf8_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__epi_an_nstop ON public.epitope_2697049_ns8_orf8_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__epi_an_start ON public.epitope_2697049_ns8_orf8_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__epi_frag_an_start ON public.epitope_2697049_ns8_orf8_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__epi_frag_an_stop ON public.epitope_2697049_ns8_orf8_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__host_tax_id ON public.epitope_2697049_ns8_orf8_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__host_tax_name ON public.epitope_2697049_ns8_orf8_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__iedb_id ON public.epitope_2697049_ns8_orf8_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__is_linear ON public.epitope_2697049_ns8_orf8_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__mhc_allele ON public.epitope_2697049_ns8_orf8_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__mhc_class_lower ON public.epitope_2697049_ns8_orf8_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__product_lower ON public.epitope_2697049_ns8_orf8_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__response_freq_pos ON public.epitope_2697049_ns8_orf8_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__seq_aa_alt ON public.epitope_2697049_ns8_orf8_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__seq_aa_orig ON public.epitope_2697049_ns8_orf8_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__start_aa_orig ON public.epitope_2697049_ns8_orf8_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__taxon_id ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__taxon_name_lower ON public.epitope_2697049_ns8_orf8_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__variant_aa_length ON public.epitope_2697049_ns8_orf8_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__variant_aa_type ON public.epitope_2697049_ns8_orf8_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_cell_type ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_epi_start ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_mhc_allele ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_product ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_host_resp_freq ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__vir_n_host_tax_id ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__virus_host_epi_stop ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_ns8_orf8_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_ns8_orf8_protein__virus_host_is_linear ON public.epitope_2697049_ns8_orf8_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__cell_type ON public.epitope_2697049_nsp10 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_an_nstop ON public.epitope_2697049_nsp10 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_an_start ON public.epitope_2697049_nsp10 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_frag_an_start ON public.epitope_2697049_nsp10 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__epi_frag_an_stop ON public.epitope_2697049_nsp10 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__host_tax_id ON public.epitope_2697049_nsp10 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__host_tax_name ON public.epitope_2697049_nsp10 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__iedb_id ON public.epitope_2697049_nsp10 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__is_linear ON public.epitope_2697049_nsp10 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__mhc_allele ON public.epitope_2697049_nsp10 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__mhc_class_lower ON public.epitope_2697049_nsp10 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__product_lower ON public.epitope_2697049_nsp10 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__response_freq_pos ON public.epitope_2697049_nsp10 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__seq_aa_alt ON public.epitope_2697049_nsp10 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__seq_aa_orig ON public.epitope_2697049_nsp10 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__start_aa_orig ON public.epitope_2697049_nsp10 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__taxon_id ON public.epitope_2697049_nsp10 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__taxon_name_lower ON public.epitope_2697049_nsp10 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__variant_aa_length ON public.epitope_2697049_nsp10 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__variant_aa_type ON public.epitope_2697049_nsp10 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_host_cell_type ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_host_epi_start ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_host_mhc_allele ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_host_product ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_host_resp_freq ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__vir_n_host_tax_id ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_epi_stop ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp10__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp10__virus_host_is_linear ON public.epitope_2697049_nsp10 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__cell_type ON public.epitope_2697049_nsp11 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_an_nstop ON public.epitope_2697049_nsp11 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_an_start ON public.epitope_2697049_nsp11 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_frag_an_start ON public.epitope_2697049_nsp11 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__epi_frag_an_stop ON public.epitope_2697049_nsp11 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__host_tax_id ON public.epitope_2697049_nsp11 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__host_tax_name ON public.epitope_2697049_nsp11 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__iedb_id ON public.epitope_2697049_nsp11 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__is_linear ON public.epitope_2697049_nsp11 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__mhc_allele ON public.epitope_2697049_nsp11 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__mhc_class_lower ON public.epitope_2697049_nsp11 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__product_lower ON public.epitope_2697049_nsp11 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__response_freq_pos ON public.epitope_2697049_nsp11 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__seq_aa_alt ON public.epitope_2697049_nsp11 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__seq_aa_orig ON public.epitope_2697049_nsp11 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__start_aa_orig ON public.epitope_2697049_nsp11 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__taxon_id ON public.epitope_2697049_nsp11 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__taxon_name_lower ON public.epitope_2697049_nsp11 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__variant_aa_length ON public.epitope_2697049_nsp11 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__variant_aa_type ON public.epitope_2697049_nsp11 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_host_cell_type ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_host_epi_start ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_host_mhc_allele ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_host_product ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_host_resp_freq ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__vir_n_host_tax_id ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_epi_stop ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp11__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp11__virus_host_is_linear ON public.epitope_2697049_nsp11 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__cell_type ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_nstop ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_an_start ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_start ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__epi_frag_an_stop ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_id ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__host_tax_name ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__iedb_id ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__is_linear ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__mhc_allele ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__mhc_class_lower ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__product_lower ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__response_freq_pos ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_alt ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__seq_aa_orig ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__start_aa_orig ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__taxon_id ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__taxon_name_lower ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_length ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__variant_aa_type ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_cell_type ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_epi_start ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_mhc_allele ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_product ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_host_resp_freq ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__vir_n_host_tax_id ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_epi_stop ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp12_rna_dependent_rna_poly__virus_host_is_linear ON public.epitope_2697049_nsp12_rna_dependent_rna_poly USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__cell_type ON public.epitope_2697049_nsp13_helicase USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__epi_an_nstop ON public.epitope_2697049_nsp13_helicase USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__epi_an_start ON public.epitope_2697049_nsp13_helicase USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__epi_frag_an_start ON public.epitope_2697049_nsp13_helicase USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__epi_frag_an_stop ON public.epitope_2697049_nsp13_helicase USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__host_tax_id ON public.epitope_2697049_nsp13_helicase USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__host_tax_name ON public.epitope_2697049_nsp13_helicase USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__iedb_id ON public.epitope_2697049_nsp13_helicase USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__is_linear ON public.epitope_2697049_nsp13_helicase USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__mhc_allele ON public.epitope_2697049_nsp13_helicase USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__mhc_class_lower ON public.epitope_2697049_nsp13_helicase USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__product_lower ON public.epitope_2697049_nsp13_helicase USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__response_freq_pos ON public.epitope_2697049_nsp13_helicase USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__seq_aa_alt ON public.epitope_2697049_nsp13_helicase USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__seq_aa_orig ON public.epitope_2697049_nsp13_helicase USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__start_aa_orig ON public.epitope_2697049_nsp13_helicase USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__taxon_id ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__taxon_name_lower ON public.epitope_2697049_nsp13_helicase USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__variant_aa_length ON public.epitope_2697049_nsp13_helicase USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__variant_aa_type ON public.epitope_2697049_nsp13_helicase USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_host_cell_type ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_host_epi_start ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_host_mhc_allele ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_host_product ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_host_resp_freq ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__vir_n_host_tax_id ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__virus_host_epi_stop ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp13_helicase__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp13_helicase__virus_host_is_linear ON public.epitope_2697049_nsp13_helicase USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__cell_type ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_an_nstop ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_an_start ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_start ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__epi_frag_an_stop ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__host_tax_id ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__host_tax_name ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__iedb_id ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__is_linear ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__mhc_allele ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__mhc_class_lower ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__product_lower ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__response_freq_pos ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_alt ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__seq_aa_orig ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__start_aa_orig ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__taxon_id ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__taxon_name_lower ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_length ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__variant_aa_type ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_cell_type ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_epi_start ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_mhc_allele ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_product ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_host_resp_freq ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__vir_n_host_tax_id ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__virus_host_epi_stop ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp14_3_to_5_exonuclease__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp14_3_to_5_exonuclease__virus_host_is_linear ON public.epitope_2697049_nsp14_3_to_5_exonuclease USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__cell_type ON public.epitope_2697049_nsp15_endornase USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__epi_an_nstop ON public.epitope_2697049_nsp15_endornase USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__epi_an_start ON public.epitope_2697049_nsp15_endornase USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__epi_frag_an_start ON public.epitope_2697049_nsp15_endornase USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__epi_frag_an_stop ON public.epitope_2697049_nsp15_endornase USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__host_tax_id ON public.epitope_2697049_nsp15_endornase USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__host_tax_name ON public.epitope_2697049_nsp15_endornase USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__iedb_id ON public.epitope_2697049_nsp15_endornase USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__is_linear ON public.epitope_2697049_nsp15_endornase USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__mhc_allele ON public.epitope_2697049_nsp15_endornase USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__mhc_class_lower ON public.epitope_2697049_nsp15_endornase USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__product_lower ON public.epitope_2697049_nsp15_endornase USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__response_freq_pos ON public.epitope_2697049_nsp15_endornase USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__seq_aa_alt ON public.epitope_2697049_nsp15_endornase USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__seq_aa_orig ON public.epitope_2697049_nsp15_endornase USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__start_aa_orig ON public.epitope_2697049_nsp15_endornase USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__taxon_id ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__taxon_name_lower ON public.epitope_2697049_nsp15_endornase USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__variant_aa_length ON public.epitope_2697049_nsp15_endornase USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__variant_aa_type ON public.epitope_2697049_nsp15_endornase USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_host_cell_type ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_host_epi_start ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_host_mhc_allele ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_host_product ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_host_resp_freq ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__vir_n_host_tax_id ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__virus_host_epi_stop ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp15_endornase__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp15_endornase__virus_host_is_linear ON public.epitope_2697049_nsp15_endornase USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__cell_type ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_nstop ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_an_start ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_start ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__epi_frag_an_stop ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_id ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__host_tax_name ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__iedb_id ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__is_linear ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_allele ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__mhc_class_lower ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__product_lower ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__response_freq_pos ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_alt ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__seq_aa_orig ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__start_aa_orig ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_id ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__taxon_name_lower ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_length ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__variant_aa_type ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_cell_type ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_epi_start ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_mhc_allele ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_product ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_host_resp_freq ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__vir_n_host_tax_id ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_epi_stop ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp16_2_o_ribose_methyltrans__virus_host_is_linear ON public.epitope_2697049_nsp16_2_o_ribose_methyltrans USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__cell_type ON public.epitope_2697049_nsp1_leader_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__epi_an_nstop ON public.epitope_2697049_nsp1_leader_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__epi_an_start ON public.epitope_2697049_nsp1_leader_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__epi_frag_an_start ON public.epitope_2697049_nsp1_leader_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__epi_frag_an_stop ON public.epitope_2697049_nsp1_leader_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__host_tax_id ON public.epitope_2697049_nsp1_leader_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__host_tax_name ON public.epitope_2697049_nsp1_leader_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__iedb_id ON public.epitope_2697049_nsp1_leader_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__is_linear ON public.epitope_2697049_nsp1_leader_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__mhc_allele ON public.epitope_2697049_nsp1_leader_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__mhc_class_lower ON public.epitope_2697049_nsp1_leader_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__product_lower ON public.epitope_2697049_nsp1_leader_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__response_freq_pos ON public.epitope_2697049_nsp1_leader_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__seq_aa_alt ON public.epitope_2697049_nsp1_leader_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__seq_aa_orig ON public.epitope_2697049_nsp1_leader_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__start_aa_orig ON public.epitope_2697049_nsp1_leader_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__taxon_id ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__taxon_name_lower ON public.epitope_2697049_nsp1_leader_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__variant_aa_length ON public.epitope_2697049_nsp1_leader_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__variant_aa_type ON public.epitope_2697049_nsp1_leader_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_cell_type ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_epi_start ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_mhc_allele ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_product ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_host_resp_freq ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__vir_n_host_tax_id ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__virus_host_epi_stop ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp1_leader_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp1_leader_protein__virus_host_is_linear ON public.epitope_2697049_nsp1_leader_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__cell_type ON public.epitope_2697049_nsp2 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_an_nstop ON public.epitope_2697049_nsp2 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_an_start ON public.epitope_2697049_nsp2 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_frag_an_start ON public.epitope_2697049_nsp2 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__epi_frag_an_stop ON public.epitope_2697049_nsp2 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__host_tax_id ON public.epitope_2697049_nsp2 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__host_tax_name ON public.epitope_2697049_nsp2 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__iedb_id ON public.epitope_2697049_nsp2 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__is_linear ON public.epitope_2697049_nsp2 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__mhc_allele ON public.epitope_2697049_nsp2 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__mhc_class_lower ON public.epitope_2697049_nsp2 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__product_lower ON public.epitope_2697049_nsp2 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__response_freq_pos ON public.epitope_2697049_nsp2 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__seq_aa_alt ON public.epitope_2697049_nsp2 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__seq_aa_orig ON public.epitope_2697049_nsp2 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__start_aa_orig ON public.epitope_2697049_nsp2 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__taxon_id ON public.epitope_2697049_nsp2 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__taxon_name_lower ON public.epitope_2697049_nsp2 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__variant_aa_length ON public.epitope_2697049_nsp2 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__variant_aa_type ON public.epitope_2697049_nsp2 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_host_cell_type ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_host_epi_start ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_host_mhc_allele ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_host_product ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_host_resp_freq ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__vir_n_host_tax_id ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_epi_stop ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp2__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp2__virus_host_is_linear ON public.epitope_2697049_nsp2 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__cell_type ON public.epitope_2697049_nsp3 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_an_nstop ON public.epitope_2697049_nsp3 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_an_start ON public.epitope_2697049_nsp3 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_frag_an_start ON public.epitope_2697049_nsp3 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__epi_frag_an_stop ON public.epitope_2697049_nsp3 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__host_tax_id ON public.epitope_2697049_nsp3 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__host_tax_name ON public.epitope_2697049_nsp3 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__iedb_id ON public.epitope_2697049_nsp3 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__is_linear ON public.epitope_2697049_nsp3 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__mhc_allele ON public.epitope_2697049_nsp3 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__mhc_class_lower ON public.epitope_2697049_nsp3 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__product_lower ON public.epitope_2697049_nsp3 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__response_freq_pos ON public.epitope_2697049_nsp3 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__seq_aa_alt ON public.epitope_2697049_nsp3 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__seq_aa_orig ON public.epitope_2697049_nsp3 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__start_aa_orig ON public.epitope_2697049_nsp3 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__taxon_id ON public.epitope_2697049_nsp3 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__taxon_name_lower ON public.epitope_2697049_nsp3 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__variant_aa_length ON public.epitope_2697049_nsp3 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__variant_aa_type ON public.epitope_2697049_nsp3 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_host_cell_type ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_host_epi_start ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_host_mhc_allele ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_host_product ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_host_resp_freq ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__vir_n_host_tax_id ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_epi_stop ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp3__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp3__virus_host_is_linear ON public.epitope_2697049_nsp3 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__cell_type ON public.epitope_2697049_nsp4 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_an_nstop ON public.epitope_2697049_nsp4 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_an_start ON public.epitope_2697049_nsp4 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_frag_an_start ON public.epitope_2697049_nsp4 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__epi_frag_an_stop ON public.epitope_2697049_nsp4 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__host_tax_id ON public.epitope_2697049_nsp4 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__host_tax_name ON public.epitope_2697049_nsp4 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__iedb_id ON public.epitope_2697049_nsp4 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__is_linear ON public.epitope_2697049_nsp4 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__mhc_allele ON public.epitope_2697049_nsp4 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__mhc_class_lower ON public.epitope_2697049_nsp4 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__product_lower ON public.epitope_2697049_nsp4 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__response_freq_pos ON public.epitope_2697049_nsp4 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__seq_aa_alt ON public.epitope_2697049_nsp4 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__seq_aa_orig ON public.epitope_2697049_nsp4 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__start_aa_orig ON public.epitope_2697049_nsp4 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__taxon_id ON public.epitope_2697049_nsp4 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__taxon_name_lower ON public.epitope_2697049_nsp4 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__variant_aa_length ON public.epitope_2697049_nsp4 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__variant_aa_type ON public.epitope_2697049_nsp4 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_host_cell_type ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_host_epi_start ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_host_mhc_allele ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_host_product ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_host_resp_freq ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__vir_n_host_tax_id ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_epi_stop ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp4__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp4__virus_host_is_linear ON public.epitope_2697049_nsp4 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__cell_type ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_an_nstop ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_an_start ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_start ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__epi_frag_an_stop ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__host_tax_id ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__host_tax_name ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__iedb_id ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__is_linear ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__mhc_allele ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__mhc_class_lower ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__product_lower ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__response_freq_pos ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__seq_aa_alt ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__seq_aa_orig ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__start_aa_orig ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__taxon_id ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__taxon_name_lower ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__variant_aa_length ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__variant_aa_type ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_cell_type ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_epi_start ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_mhc_allele ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_product ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_host_resp_freq ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__vir_n_host_tax_id ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__virus_host_epi_stop ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp5_3c_like_proteinase__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp5_3c_like_proteinase__virus_host_is_linear ON public.epitope_2697049_nsp5_3c_like_proteinase USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__cell_type ON public.epitope_2697049_nsp6 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_an_nstop ON public.epitope_2697049_nsp6 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_an_start ON public.epitope_2697049_nsp6 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_frag_an_start ON public.epitope_2697049_nsp6 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__epi_frag_an_stop ON public.epitope_2697049_nsp6 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__host_tax_id ON public.epitope_2697049_nsp6 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__host_tax_name ON public.epitope_2697049_nsp6 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__iedb_id ON public.epitope_2697049_nsp6 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__is_linear ON public.epitope_2697049_nsp6 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__mhc_allele ON public.epitope_2697049_nsp6 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__mhc_class_lower ON public.epitope_2697049_nsp6 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__product_lower ON public.epitope_2697049_nsp6 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__response_freq_pos ON public.epitope_2697049_nsp6 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__seq_aa_alt ON public.epitope_2697049_nsp6 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__seq_aa_orig ON public.epitope_2697049_nsp6 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__start_aa_orig ON public.epitope_2697049_nsp6 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__taxon_id ON public.epitope_2697049_nsp6 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__taxon_name_lower ON public.epitope_2697049_nsp6 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__variant_aa_length ON public.epitope_2697049_nsp6 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__variant_aa_type ON public.epitope_2697049_nsp6 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_host_cell_type ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_host_epi_start ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_host_mhc_allele ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_host_product ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_host_resp_freq ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__vir_n_host_tax_id ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_epi_stop ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp6__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp6__virus_host_is_linear ON public.epitope_2697049_nsp6 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__cell_type ON public.epitope_2697049_nsp7 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_an_nstop ON public.epitope_2697049_nsp7 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_an_start ON public.epitope_2697049_nsp7 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_frag_an_start ON public.epitope_2697049_nsp7 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__epi_frag_an_stop ON public.epitope_2697049_nsp7 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__host_tax_id ON public.epitope_2697049_nsp7 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__host_tax_name ON public.epitope_2697049_nsp7 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__iedb_id ON public.epitope_2697049_nsp7 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__is_linear ON public.epitope_2697049_nsp7 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__mhc_allele ON public.epitope_2697049_nsp7 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__mhc_class_lower ON public.epitope_2697049_nsp7 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__product_lower ON public.epitope_2697049_nsp7 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__response_freq_pos ON public.epitope_2697049_nsp7 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__seq_aa_alt ON public.epitope_2697049_nsp7 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__seq_aa_orig ON public.epitope_2697049_nsp7 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__start_aa_orig ON public.epitope_2697049_nsp7 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__taxon_id ON public.epitope_2697049_nsp7 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__taxon_name_lower ON public.epitope_2697049_nsp7 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__variant_aa_length ON public.epitope_2697049_nsp7 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__variant_aa_type ON public.epitope_2697049_nsp7 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_host_cell_type ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_host_epi_start ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_host_mhc_allele ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_host_product ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_host_resp_freq ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__vir_n_host_tax_id ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_epi_stop ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp7__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp7__virus_host_is_linear ON public.epitope_2697049_nsp7 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__cell_type ON public.epitope_2697049_nsp8 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_an_nstop ON public.epitope_2697049_nsp8 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_an_start ON public.epitope_2697049_nsp8 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_frag_an_start ON public.epitope_2697049_nsp8 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__epi_frag_an_stop ON public.epitope_2697049_nsp8 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__host_tax_id ON public.epitope_2697049_nsp8 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__host_tax_name ON public.epitope_2697049_nsp8 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__iedb_id ON public.epitope_2697049_nsp8 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__is_linear ON public.epitope_2697049_nsp8 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__mhc_allele ON public.epitope_2697049_nsp8 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__mhc_class_lower ON public.epitope_2697049_nsp8 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__product_lower ON public.epitope_2697049_nsp8 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__response_freq_pos ON public.epitope_2697049_nsp8 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__seq_aa_alt ON public.epitope_2697049_nsp8 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__seq_aa_orig ON public.epitope_2697049_nsp8 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__start_aa_orig ON public.epitope_2697049_nsp8 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__taxon_id ON public.epitope_2697049_nsp8 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__taxon_name_lower ON public.epitope_2697049_nsp8 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__variant_aa_length ON public.epitope_2697049_nsp8 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__variant_aa_type ON public.epitope_2697049_nsp8 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_host_cell_type ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_host_epi_start ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_host_mhc_allele ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_host_product ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_host_resp_freq ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__vir_n_host_tax_id ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_epi_stop ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp8__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp8__virus_host_is_linear ON public.epitope_2697049_nsp8 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__cell_type ON public.epitope_2697049_nsp9 USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_an_nstop ON public.epitope_2697049_nsp9 USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_an_start ON public.epitope_2697049_nsp9 USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_frag_an_start ON public.epitope_2697049_nsp9 USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__epi_frag_an_stop ON public.epitope_2697049_nsp9 USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__host_tax_id ON public.epitope_2697049_nsp9 USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__host_tax_name ON public.epitope_2697049_nsp9 USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__iedb_id ON public.epitope_2697049_nsp9 USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__is_linear ON public.epitope_2697049_nsp9 USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__mhc_allele ON public.epitope_2697049_nsp9 USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__mhc_class_lower ON public.epitope_2697049_nsp9 USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__product_lower ON public.epitope_2697049_nsp9 USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__response_freq_pos ON public.epitope_2697049_nsp9 USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__seq_aa_alt ON public.epitope_2697049_nsp9 USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__seq_aa_orig ON public.epitope_2697049_nsp9 USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__start_aa_orig ON public.epitope_2697049_nsp9 USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__taxon_id ON public.epitope_2697049_nsp9 USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__taxon_name_lower ON public.epitope_2697049_nsp9 USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__variant_aa_length ON public.epitope_2697049_nsp9 USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__variant_aa_type ON public.epitope_2697049_nsp9 USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_host_cell_type ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_host_epi_start ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_host_mhc_allele ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_host_product ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_host_resp_freq ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__vir_n_host_tax_id ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_epi_stop ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_nsp9__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_nsp9__virus_host_is_linear ON public.epitope_2697049_nsp9 USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__cell_type ON public.epitope_2697049_orf10_protein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__epi_an_nstop ON public.epitope_2697049_orf10_protein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__epi_an_start ON public.epitope_2697049_orf10_protein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__epi_frag_an_start ON public.epitope_2697049_orf10_protein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__epi_frag_an_stop ON public.epitope_2697049_orf10_protein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__host_tax_id ON public.epitope_2697049_orf10_protein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__host_tax_name ON public.epitope_2697049_orf10_protein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__iedb_id ON public.epitope_2697049_orf10_protein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__is_linear ON public.epitope_2697049_orf10_protein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__mhc_allele ON public.epitope_2697049_orf10_protein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__mhc_class_lower ON public.epitope_2697049_orf10_protein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__product_lower ON public.epitope_2697049_orf10_protein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__response_freq_pos ON public.epitope_2697049_orf10_protein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__seq_aa_alt ON public.epitope_2697049_orf10_protein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__seq_aa_orig ON public.epitope_2697049_orf10_protein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__start_aa_orig ON public.epitope_2697049_orf10_protein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__taxon_id ON public.epitope_2697049_orf10_protein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__taxon_name_lower ON public.epitope_2697049_orf10_protein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__variant_aa_length ON public.epitope_2697049_orf10_protein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__variant_aa_type ON public.epitope_2697049_orf10_protein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_host_cell_type ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_host_epi_start ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_host_mhc_allele ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_host_product ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_host_resp_freq ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__vir_n_host_tax_id ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__virus_host_epi_stop ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf10_protein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf10_protein__virus_host_is_linear ON public.epitope_2697049_orf10_protein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__cell_type ON public.epitope_2697049_orf1a_polyprotein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__epi_an_nstop ON public.epitope_2697049_orf1a_polyprotein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__epi_an_start ON public.epitope_2697049_orf1a_polyprotein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__epi_frag_an_start ON public.epitope_2697049_orf1a_polyprotein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__epi_frag_an_stop ON public.epitope_2697049_orf1a_polyprotein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__host_tax_id ON public.epitope_2697049_orf1a_polyprotein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__host_tax_name ON public.epitope_2697049_orf1a_polyprotein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__iedb_id ON public.epitope_2697049_orf1a_polyprotein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__is_linear ON public.epitope_2697049_orf1a_polyprotein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__mhc_allele ON public.epitope_2697049_orf1a_polyprotein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__mhc_class_lower ON public.epitope_2697049_orf1a_polyprotein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__product_lower ON public.epitope_2697049_orf1a_polyprotein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__response_freq_pos ON public.epitope_2697049_orf1a_polyprotein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__seq_aa_alt ON public.epitope_2697049_orf1a_polyprotein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__seq_aa_orig ON public.epitope_2697049_orf1a_polyprotein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__start_aa_orig ON public.epitope_2697049_orf1a_polyprotein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__taxon_id ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__taxon_name_lower ON public.epitope_2697049_orf1a_polyprotein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__variant_aa_length ON public.epitope_2697049_orf1a_polyprotein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__variant_aa_type ON public.epitope_2697049_orf1a_polyprotein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_cell_type ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_epi_start ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_mhc_allele ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_product ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_host_resp_freq ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__vir_n_host_tax_id ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__virus_host_epi_stop ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1a_polyprotein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1a_polyprotein__virus_host_is_linear ON public.epitope_2697049_orf1a_polyprotein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__cell_type ON public.epitope_2697049_orf1ab_polyprotein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_an_nstop ON public.epitope_2697049_orf1ab_polyprotein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_an_start ON public.epitope_2697049_orf1ab_polyprotein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_frag_an_start ON public.epitope_2697049_orf1ab_polyprotein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__epi_frag_an_stop ON public.epitope_2697049_orf1ab_polyprotein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__host_tax_id ON public.epitope_2697049_orf1ab_polyprotein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__host_tax_name ON public.epitope_2697049_orf1ab_polyprotein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__iedb_id ON public.epitope_2697049_orf1ab_polyprotein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__is_linear ON public.epitope_2697049_orf1ab_polyprotein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__mhc_allele ON public.epitope_2697049_orf1ab_polyprotein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__mhc_class_lower ON public.epitope_2697049_orf1ab_polyprotein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__product_lower ON public.epitope_2697049_orf1ab_polyprotein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__response_freq_pos ON public.epitope_2697049_orf1ab_polyprotein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__seq_aa_alt ON public.epitope_2697049_orf1ab_polyprotein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__seq_aa_orig ON public.epitope_2697049_orf1ab_polyprotein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__start_aa_orig ON public.epitope_2697049_orf1ab_polyprotein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__taxon_id ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__taxon_name_lower ON public.epitope_2697049_orf1ab_polyprotein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__variant_aa_length ON public.epitope_2697049_orf1ab_polyprotein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__variant_aa_type ON public.epitope_2697049_orf1ab_polyprotein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_cell_type ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_epi_start ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_mhc_allele ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_product ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_host_resp_freq ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__vir_n_host_tax_id ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__virus_host_epi_stop ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_orf1ab_polyprotein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_orf1ab_polyprotein__virus_host_is_linear ON public.epitope_2697049_orf1ab_polyprotein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__cell_type ON public.epitope_2697049_spike_surface_glycoprotein USING btree (((cell_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__epi_an_nstop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_an_nstop ON public.epitope_2697049_spike_surface_glycoprotein USING btree (epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__epi_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_an_start ON public.epitope_2697049_spike_surface_glycoprotein USING btree (epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__epi_frag_an_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_frag_an_start ON public.epitope_2697049_spike_surface_glycoprotein USING btree (epi_frag_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__epi_frag_an_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__epi_frag_an_stop ON public.epitope_2697049_spike_surface_glycoprotein USING btree (epi_frag_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__host_tax_id ON public.epitope_2697049_spike_surface_glycoprotein USING btree (host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__host_tax_name; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__host_tax_name ON public.epitope_2697049_spike_surface_glycoprotein USING btree (lower((host_taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__iedb_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__iedb_id ON public.epitope_2697049_spike_surface_glycoprotein USING btree (iedb_epitope_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__is_linear ON public.epitope_2697049_spike_surface_glycoprotein USING btree (is_linear) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__mhc_allele ON public.epitope_2697049_spike_surface_glycoprotein USING btree (((mhc_allele)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__mhc_class_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__mhc_class_lower ON public.epitope_2697049_spike_surface_glycoprotein USING btree (lower((mhc_class)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__product_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__product_lower ON public.epitope_2697049_spike_surface_glycoprotein USING btree (lower((product)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__response_freq_pos; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__response_freq_pos ON public.epitope_2697049_spike_surface_glycoprotein USING btree (response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__seq_aa_alt; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__seq_aa_alt ON public.epitope_2697049_spike_surface_glycoprotein USING btree (((sequence_aa_alternative)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__seq_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__seq_aa_orig ON public.epitope_2697049_spike_surface_glycoprotein USING btree (((sequence_aa_original)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__start_aa_orig; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__start_aa_orig ON public.epitope_2697049_spike_surface_glycoprotein USING btree (start_aa_original) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__taxon_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__taxon_id ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__taxon_name_lower; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__taxon_name_lower ON public.epitope_2697049_spike_surface_glycoprotein USING btree (lower((taxon_name)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__variant_aa_length; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__variant_aa_length ON public.epitope_2697049_spike_surface_glycoprotein USING btree (variant_aa_length) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__variant_aa_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__variant_aa_type ON public.epitope_2697049_spike_surface_glycoprotein USING btree (((variant_aa_type)::text)) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_host_cell_type; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_cell_type ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, cell_type) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_host_epi_start; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_epi_start ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, epi_annotation_start) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_host_mhc_allele; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_mhc_allele ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, mhc_allele) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_host_product; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_product ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, product) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_host_resp_freq; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_host_resp_freq ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, response_frequency_pos) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__vir_n_host_tax_id; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__vir_n_host_tax_id ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__virus_host_epi_stop; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__virus_host_epi_stop ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, epi_annotation_stop) WITH (fillfactor='100');


--
-- Name: epi_2697049_spike_surface_glycoprotein__virus_host_is_linear; Type: INDEX; Schema: public; Owner: geco
--

CREATE INDEX epi_2697049_spike_surface_glycoprotein__virus_host_is_linear ON public.epitope_2697049_spike_surface_glycoprotein USING btree (taxon_id, host_taxon_id, is_linear) WITH (fillfactor='100');


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

