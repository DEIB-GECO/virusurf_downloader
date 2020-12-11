#!/bin/bash

START=$(date +%s.%N)
echo "writing CSV"
python main.py vcm_dev_bulk_inserts True geco geco78 5432 import coguk
echo "copying CSV into DB"

psql -h localhost -U geco -d vcm_dev_bulk_inserts -c "\copy annotation FROM '/home/alfonsi/virusurf_downloader/generated/COG-UK_sars_cov_2/other/ann.csv' WITH DELIMITER ',' CSV HEADER;"

psql -h localhost -U geco -d vcm_dev_bulk_inserts -c "\copy aminoacid_variant FROM '/home/alfonsi/virusurf_downloader/generated/COG-UK_sars_cov_2/other/aa.csv' WITH DELIMITER ',' CSV HEADER;"

psql -h localhost -U geco -d vcm_dev_bulk_inserts -c "\copy nucleotide_variant FROM '/home/alfonsi/virusurf_downloader/generated/COG-UK_sars_cov_2/other/nuc_var.csv' WITH DELIMITER ',' CSV HEADER;"

psql -h localhost -U geco -d vcm_dev_bulk_inserts -c "\copy variant_impact FROM '/home/alfonsi/virusurf_downloader/generated/COG-UK_sars_cov_2/other/var_imp.csv' WITH DELIMITER ',' CSV HEADER;"

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF

