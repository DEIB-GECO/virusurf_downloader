update sequence
set
    accession_id = CONCAT(accession_id, '_', (select right(virus.taxon_name,1) from virus where virus.virus_id = sequence.virus_id)),
    alternative_accession_id = CONCAT(alternative_accession_id, '_', (select right(virus.taxon_name,1) from virus where virus.virus_id = sequence.virus_id))
from virus
where accession_id in(
    select distinct accession_id
    from sequence
    group by accession_id
    having count(accession_id)>1
);