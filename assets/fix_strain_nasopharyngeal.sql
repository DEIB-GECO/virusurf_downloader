-- Corrects strain name == nasopharyngeal
begin;

select distinct s.host_sample_id from "sequence" s natural join host_sample hs where s.strain_name = 'nasopharyngeal' and hs.isolation_source is null;
-- it was only one id: 14753

select distinct s.strain_name, hs.isolation_source from sequence s natural join host_sample hs where s.host_sample_id = 14753;

-- check if host_sample_id = 14753 selects the same sequences as s.strain_name = 'nasopharyngeal' and hs.isolation_source is null
select count(a.sequence_id) from (
select s.sequence_id from sequence s where s.host_sample_id = 14753
except
select s.sequence_id from "sequence" s natural join host_sample hs where s.strain_name = 'nasopharyngeal' and hs.isolation_source is null
) as a; -- if so it should be 0. Then we can use host_id to identify those sequences during the two updates.


update host_sample set isolation_source = 'nasopharyngeal' where host_sample_id = 14753;
update sequence set strain_name = null where host_sample_id = 14753;

-- check result
select distinct s.strain_name, hs.isolation_source from sequence s natural join host_sample hs where s.host_sample_id = 14753;

--rollback;

commit;


