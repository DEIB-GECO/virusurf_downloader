CREATE TABLE IF NOT EXISTS lineages (
    accession_id varchar primary key,
    lineage varchar,
    conflict float4,
    pangoLEARN_version varchar,
    pango_version varchar,
    status varchar,
    note varchar
);

\copy lineages FROM 'pathtocsvfile' WITH DELIMITER ',' CSV HEADER;

update sequence s set lineage = l.lineage from lineages l where s.accession_id = l.accession_id and l.lineage != 'None';

DROP TABLE lineages;

