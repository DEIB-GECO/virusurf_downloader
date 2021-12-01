-- CREATE OR REPLACE VIEW public.annotation_cds
-- AS SELECT annotation.annotation_id,
--     annotation.sequence_id,
--     annotation.start,
--     annotation.stop,
--     annotation.gene_name,
--     annotation.product,
--     annotation.external_reference,
--     annotation.aminoacid_sequence
--    FROM annotation
--   WHERE annotation.feature_type::text = 'CDS'::text;
--
--
-- CREATE OR REPLACE VIEW public.annotation_view
-- AS SELECT annotation.sequence_id,
--     annotation.product AS annotation_view_product,
--     annotation.aminoacid_sequence AS annotation_view_aminoacid_sequence,
--     annotation.annotation_nucleotide_sequence AS annotation_view_nucleotide_sequence
--    FROM annotation
--   WHERE annotation.product IS NOT NULL AND
--         (annotation.aminoacid_sequence IS NOT NULL OR annotation.annotation_nucleotide_sequence IS NOT NULL);


CREATE OR REPLACE VIEW public.host_sample_view
AS SELECT host_sample.*,
    host_specie.host_taxon_id, host_specie.host_taxon_name
   FROM host_sample
     JOIN host_specie USING (host_id);


CREATE OR REPLACE VIEW public.nucleotide_variant_limited
AS SELECT nucleotide_variant.*
   FROM nucleotide_variant
  WHERE nucleotide_variant.variant_length <= 20;


DROP MATERIALIZED VIEW IF EXISTS public.nucleotide_variant_annotated;
CREATE MATERIALIZED VIEW public.nucleotide_variant_annotated
WITH (
    FILLFACTOR = 100
)
TABLESPACE default_ts
AS
SELECT nc.nucleotide_variant_id,
    nc.sequence_id,
    nc.variant_type,
    nc.start_original,
    CASE nc.variant_type  when 'SUB' then  nc.sequence_original    ELSE null END as sequence_original,
    CASE nc.variant_type  when 'SUB' then  nc.sequence_alternative ELSE null END as sequence_alternative,
    nc.variant_length,
    ann.gene_name AS n_gene_name
   FROM nucleotide_variant nc
     LEFT JOIN annotation ann ON nc.start_original >= ann.start AND nc.start_original <= ann.stop AND nc.sequence_id = ann.sequence_id AND ann.feature_type = 'gene'
WITH DATA;



-- LOCATION
DROP MATERIALIZED VIEW IF EXISTS public.location;
-- call refresh public.variant_observation before first-time usage!
CREATE MATERIALIZED VIEW IF NOT EXISTS public.location
    TABLESPACE default_ts AS
select distinct province as loc from host_sample where province is not null and province <> '' UNION
select distinct region as loc from host_sample where region is not null and region <> '' UNION
select distinct country as loc from host_sample where country is not null and country <> '' UNION
select distinct geo_group as loc from host_sample where geo_group is not null and geo_group <> '' WITH NO DATA;
-- call refresh public.location before first-time usage!

CREATE EXTENSION IF NOT EXISTS btree_gin;

CREATE INDEX IF NOT EXISTS idx_loc on public.location USING GIN (loc);


-- VARIANT OBSERVATION VIEW
DROP MATERIALIZED VIEW IF EXISTS public.variant_observation;
-- call refresh public.variant_observation before first-time usage!
CREATE MATERIALIZED VIEW IF NOT EXISTS public.variant_observation
TABLESPACE default_ts AS
    select * from (
        WITH
        host_sample_periods as (
                select  host_sample.*,
                left(collection_date,7) as month
                from host_sample
                where coll_date_precision >= 1
        ),
        change_location as (
                select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
                loc as obs_location
                from sequence S
                natural join host_sample H
                join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
                natural join annotation A
                natural join aminoacid_variant V
                --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc
    ),
    time_period as (
            select distinct month as collect_period from host_sample_periods
    ),
    mutated as (
            select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
            loc as obs_location,
            month as collect_period,
            count(*) as mutated_sequences
            from sequence S
            natural join host_sample_periods  H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            natural join annotation A
            natural join aminoacid_variant V
            --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc, month
    ),
    totals as (
            select
            loc as obs_location,
            month as collect_period,
            count(*) as tot_sequences
            from sequence S
            natural join host_sample_periods H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            group by loc, month
    )
    SELECT change_location.change, change_location.obs_location,
    time_period.collect_period, 'month' as period_type,
    coalesce(mutated.mutated_sequences,0) as mutated_sequences,
    coalesce(totals.tot_sequences,0) as tot_sequences
    from change_location
    cross join time_period
    left join mutated on (change_location.change = mutated.change and change_location.obs_location = mutated.obs_location and time_period.collect_period = mutated.collect_period)
    left join totals on (change_location.obs_location = totals.obs_location and time_period.collect_period = totals.collect_period)
    ) A
    UNION
    select * from
    (
            WITH
            host_sample_periods as (
                    select  host_sample.*,
                    case
                    when right(collection_date,2) >= '01' and right(collection_date,2) <= '15'
                            then concat(left(collection_date,8),'01|15')
                    when right(collection_date,2) >= '16' and right(collection_date,2) <= '31'
                            then concat(left(collection_date,8),'16|31')
                    end biweek
                    from host_sample
                    where coll_date_precision >= 2
            ),
            change_location as (
                    select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
                    loc as obs_location
                    from sequence S
                    natural join host_sample H
                    join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
                    natural join annotation A
                    natural join aminoacid_variant V
                    --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc
    ),
    time_period as (
            select distinct biweek as collect_period from host_sample_periods
    ),
    mutated as (
            select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
            loc as obs_location, biweek as collect_period,
            count(*) as mutated_sequences
            from sequence S
            natural join host_sample_periods  H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            natural join annotation A
            natural join aminoacid_variant V
            --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc, biweek
    ),
    totals as (
            select
            loc as obs_location,
            biweek as collect_period,
            count(*) as tot_sequences
            from sequence S
            natural join host_sample_periods H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            group by loc, biweek
    )
    SELECT change_location.change, change_location.obs_location,
    time_period.collect_period, 'biweek' as period_type,
    coalesce(mutated.mutated_sequences,0) as mutated_sequences,
    coalesce(totals.tot_sequences,0) as tot_sequences
    from change_location
    cross join time_period
    left join mutated on (change_location.change = mutated.change and change_location.obs_location = mutated.obs_location and time_period.collect_period = mutated.collect_period)
    left join totals on (change_location.obs_location = totals.obs_location and time_period.collect_period = totals.collect_period)
    ) B
    UNION
    select * from (
            WITH
            host_sample_periods as (
                    select  host_sample.*,
                    case
                    when right(collection_date,2) >= '01' and right(collection_date,2) <= '07'
                            then concat(left(collection_date,8),'01|07')
                    when right(collection_date,2) >= '08' and right(collection_date,2) <= '15'
                            then concat(left(collection_date,8),'08|15')
                    when right(collection_date,2) >= '16' and right(collection_date,2) <= '23'
                            then concat(left(collection_date,8),'16|23')
                    when right(collection_date,2) >= '24' and right(collection_date,2) <= '31'
                            then concat(left(collection_date,8),'24|31')
                    end week
                    from host_sample
                    where coll_date_precision >= 2
            ),
            change_location as (
                    select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
                    loc as obs_location
                    from sequence S
                    natural join host_sample H
                    join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
                    natural join annotation A
                    natural join aminoacid_variant V
                    --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc
    ),
    time_period as (
            select distinct week as collect_period from host_sample_periods
    ),
    mutated as (
            select concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative) as change,
            loc as obs_location, week as collect_period,
            count(*) as mutated_sequences
            from sequence S
            natural join host_sample_periods  H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            natural join annotation A
            natural join aminoacid_variant V
            --where loc in ('United Kingdom','Switzerland')
    --and (concat(split_part(product, ' ', 1), '_', sequence_aa_original, start_aa_original, sequence_aa_alternative) in ('N_R203K'))
    group by concat(split_part(product, ' ', 1),'_', sequence_aa_original, start_aa_original, sequence_aa_alternative), loc, week
    ),
    totals as (
            select
            loc as obs_location,
            week as collect_period,
            count(*) as tot_sequences
            from sequence S
            natural join host_sample_periods H
            join location L on (H.region = L.loc or H.country = L.loc or H.geo_group = L.loc)
            group by loc, week
    )
    SELECT change_location.change, change_location.obs_location,
    time_period.collect_period, 'week' as period_type,
    coalesce(mutated.mutated_sequences,0) as mutated_sequences,
    coalesce(totals.tot_sequences,0) as tot_sequences
    from change_location
    cross join time_period
    left join mutated on (change_location.change = mutated.change and change_location.obs_location = mutated.obs_location and time_period.collect_period = mutated.collect_period)
    left join totals on (change_location.obs_location = totals.obs_location and time_period.collect_period = totals.collect_period)
    ) C
WITH NO DATA;
-- call refresh public.variant_observation before first-time usage!