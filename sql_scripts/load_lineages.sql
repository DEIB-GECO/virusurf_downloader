CREATE TABLE IF NOT EXISTS lineages (
    accession_id varchar primary key,
    lineage varchar
);
TRUNCATE lineages RESTART IDENTITY;

-- PROGRAM allows to specify a shell command in place of a standard input.
-- CUT can cut the file to keep some columns of a CSV file.
-- CUT options:
--  -d "x" : split columns separated by character x
-- -f 1,2-4 : selects columns 1 and 2-to-4
-- -- : inverts the selection. It does the complement of the columns specified through -f
\copy lineages FROM PROGRAM 'cut -d "," -f 1,2 pathtocsvfile' WITH DELIMITER ',' CSV HEADER;

update sequence s set lineage = l.lineage from lineages l where s.accession_id = l.accession_id and l.lineage != 'None';

DROP TABLE lineages;

