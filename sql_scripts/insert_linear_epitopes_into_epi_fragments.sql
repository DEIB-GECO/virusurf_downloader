-- This script is called automatically after each epitope import process.
insert into epitope_fragment(epitope_id, epi_fragment_sequence, epi_frag_annotation_start, epi_frag_annotation_stop) 
select e.epitope_id, e.epitope_sequence, e.epi_annotation_start, e.epi_annotation_stop 
from epitope e 
where e.epitope_id not in (select distinct ef.epitope_id from epitope_fragment ef);
