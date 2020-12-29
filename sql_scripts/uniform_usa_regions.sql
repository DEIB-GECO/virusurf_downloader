-- select country, region,
update host_sample
set region =
case when region ilike 'SLIDELL LA' then 'Louisiana'
	when left(region, strpos(region, ',') - 1) in ('Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming', 'American Samoa', 'District of Columbia', 'Federated States of Micronesia', 'Guam', 'Marshall Islands', 'Northern Mariana Islands', 'Palau', 'Puerto Rico', 'Virgin Islands')
	then left(region, strpos(region, ',') - 1)
	when right(region, length(region)-strpos(region, ',') - 1) in ('Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming', 'American Samoa', 'District of Columbia', 'Federated States of Micronesia', 'Guam', 'Marshall Islands', 'Northern Mariana Islands', 'Palau', 'Puerto Rico', 'Virgin Islands')
	then right(region, length(region)-strpos(region, ',') - 1)
	when length(region) = 2
		then case
			when region ilike 'AL' then 'Alabama' when region ilike 'AK' then 'Alaska' when region ilike 'AZ' then 'Arizona' when region ilike 'AR' then 'Arkansas' when region ilike 'CA' then 'California' when region ilike 'CO' then 'Colorado' when region ilike 'CT' then 'Connecticut' when region ilike 'DE' then 'Delaware' when region ilike 'FL' then 'Florida' when region ilike 'GA' then 'Georgia' when region ilike 'HI' then 'Hawaii' when region ilike 'ID' then 'Idaho' when region ilike 'IL' then 'Illinois' when region ilike 'IN' then 'Indiana' when region ilike 'IA' then 'Iowa' when region ilike 'KS' then 'Kansas' when region ilike 'KY' then 'Kentucky' when region ilike 'LA' then 'Louisiana' when region ilike 'ME' then 'Maine' when region ilike 'MD' then 'Maryland' when region ilike 'MA' then 'Massachusetts' when region ilike 'MI' then 'Michigan' when region ilike 'MN' then 'Minnesota' when region ilike 'MS' then 'Mississippi' when region ilike 'MO' then 'Missouri' when region ilike 'MT' then 'Montana' when region ilike 'NE' then 'Nebraska' when region ilike 'NV' then 'Nevada' when region ilike 'NH' then 'New Hampshire' when region ilike 'NJ' then 'New Jersey' when region ilike 'NM' then 'New Mexico' when region ilike 'NY' then 'New York' when region ilike 'NC' then 'North Carolina' when region ilike 'ND' then 'North Dakota' when region ilike 'OH' then 'Ohio' when region ilike 'OK' then 'Oklahoma' when region ilike 'OR' then 'Oregon' when region ilike 'PA' then 'Pennsylvania' when region ilike 'RI' then 'Rhode Island' when region ilike 'SC' then 'South Carolina' when region ilike 'SD' then 'South Dakota' when region ilike 'TN' then 'Tennessee' when region ilike 'TX' then 'Texas' when region ilike 'UT' then 'Utah' when region ilike 'VT' then 'Vermont' when region ilike 'VA' then 'Virginia' when region ilike 'WA' then 'Washington' when region ilike 'WV' then 'West Virginia' when region ilike 'WI' then 'Wisconsin' when region ilike 'WY' then 'Wyoming' when region ilike 'AS' then 'American Samoa' when region ilike 'DC' then 'District of Columbia' when region ilike 'FM' then 'Federated States of Micronesia' when region ilike 'GU' then 'Guam' when region ilike 'MH' then 'Marshall Islands' when region ilike 'MP' then 'Northern Mariana Islands' when region ilike 'PW' then 'Palau' when region ilike 'PR' then 'Puerto Rico' when region ilike 'VI' then 'Virgin Islands'
		end
	when region ilike '%,%' and length(left(region, strpos(region, ',') - 1)) = 2
		then case
			when left(region, strpos(region, ',') - 1) ilike 'AL' then 'Alabama' when left(region, strpos(region, ',') - 1) ilike 'AK' then 'Alaska' when left(region, strpos(region, ',') - 1) ilike 'AZ' then 'Arizona' when left(region, strpos(region, ',') - 1) ilike 'AR' then 'Arkansas' when left(region, strpos(region, ',') - 1) ilike 'CA' then 'California' when left(region, strpos(region, ',') - 1) ilike 'CO' then 'Colorado' when left(region, strpos(region, ',') - 1) ilike 'CT' then 'Connecticut' when left(region, strpos(region, ',') - 1) ilike 'DE' then 'Delaware' when left(region, strpos(region, ',') - 1) ilike 'FL' then 'Florida' when left(region, strpos(region, ',') - 1) ilike 'GA' then 'Georgia' when left(region, strpos(region, ',') - 1) ilike 'HI' then 'Hawaii' when left(region, strpos(region, ',') - 1) ilike 'ID' then 'Idaho' when left(region, strpos(region, ',') - 1) ilike 'IL' then 'Illinois' when left(region, strpos(region, ',') - 1) ilike 'IN' then 'Indiana' when left(region, strpos(region, ',') - 1) ilike 'IA' then 'Iowa' when left(region, strpos(region, ',') - 1) ilike 'KS' then 'Kansas' when left(region, strpos(region, ',') - 1) ilike 'KY' then 'Kentucky' when left(region, strpos(region, ',') - 1) ilike 'LA' then 'Louisiana' when left(region, strpos(region, ',') - 1) ilike 'ME' then 'Maine' when left(region, strpos(region, ',') - 1) ilike 'MD' then 'Maryland' when left(region, strpos(region, ',') - 1) ilike 'MA' then 'Massachusetts' when left(region, strpos(region, ',') - 1) ilike 'MI' then 'Michigan' when left(region, strpos(region, ',') - 1) ilike 'MN' then 'Minnesota' when left(region, strpos(region, ',') - 1) ilike 'MS' then 'Mississippi' when left(region, strpos(region, ',') - 1) ilike 'MO' then 'Missouri' when left(region, strpos(region, ',') - 1) ilike 'MT' then 'Montana' when left(region, strpos(region, ',') - 1) ilike 'NE' then 'Nebraska' when left(region, strpos(region, ',') - 1) ilike 'NV' then 'Nevada' when left(region, strpos(region, ',') - 1) ilike 'NH' then 'New Hampshire' when left(region, strpos(region, ',') - 1) ilike 'NJ' then 'New Jersey' when left(region, strpos(region, ',') - 1) ilike 'NM' then 'New Mexico' when left(region, strpos(region, ',') - 1) ilike 'NY' then 'New York' when left(region, strpos(region, ',') - 1) ilike 'NC' then 'North Carolina' when left(region, strpos(region, ',') - 1) ilike 'ND' then 'North Dakota' when left(region, strpos(region, ',') - 1) ilike 'OH' then 'Ohio' when left(region, strpos(region, ',') - 1) ilike 'OK' then 'Oklahoma' when left(region, strpos(region, ',') - 1) ilike 'OR' then 'Oregon' when left(region, strpos(region, ',') - 1) ilike 'PA' then 'Pennsylvania' when left(region, strpos(region, ',') - 1) ilike 'RI' then 'Rhode Island' when left(region, strpos(region, ',') - 1) ilike 'SC' then 'South Carolina' when left(region, strpos(region, ',') - 1) ilike 'SD' then 'South Dakota' when left(region, strpos(region, ',') - 1) ilike 'TN' then 'Tennessee' when left(region, strpos(region, ',') - 1) ilike 'TX' then 'Texas' when left(region, strpos(region, ',') - 1) ilike 'UT' then 'Utah' when left(region, strpos(region, ',') - 1) ilike 'VT' then 'Vermont' when left(region, strpos(region, ',') - 1) ilike 'VA' then 'Virginia' when left(region, strpos(region, ',') - 1) ilike 'WA' then 'Washington' when left(region, strpos(region, ',') - 1) ilike 'WV' then 'West Virginia' when left(region, strpos(region, ',') - 1) ilike 'WI' then 'Wisconsin' when left(region, strpos(region, ',') - 1) ilike 'WY' then 'Wyoming' when left(region, strpos(region, ',') - 1) ilike 'AS' then 'American Samoa' when left(region, strpos(region, ',') - 1) ilike 'DC' then 'District of Columbia' when left(region, strpos(region, ',') - 1) ilike 'FM' then 'Federated States of Micronesia' when left(region, strpos(region, ',') - 1) ilike 'GU' then 'Guam' when left(region, strpos(region, ',') - 1) ilike 'MH' then 'Marshall Islands' when left(region, strpos(region, ',') - 1) ilike 'MP' then 'Northern Mariana Islands' when left(region, strpos(region, ',') - 1) ilike 'PW' then 'Palau' when left(region, strpos(region, ',') - 1) ilike 'PR' then 'Puerto Rico' when left(region, strpos(region, ',') - 1) ilike 'VI' then 'Virgin Islands'
		end
	when region ilike '%,%' and (length(right(region, length(region)-strpos(region, ',') - 1)) = 2 or substring(region,length(region)-2,length(region)+1) ilike ',__')
		then case
			when substring(region,length(region)-1,length(region)+1) ilike 'AL' then 'Alabama' when substring(region,length(region)-1,length(region)+1) ilike 'AK' then 'Alaska' when substring(region,length(region)-1,length(region)+1) ilike 'AZ' then 'Arizona' when substring(region,length(region)-1,length(region)+1) ilike 'AR' then 'Arkansas' when substring(region,length(region)-1,length(region)+1) ilike 'CA' then 'California' when substring(region,length(region)-1,length(region)+1) ilike 'CO' then 'Colorado' when substring(region,length(region)-1,length(region)+1) ilike 'CT' then 'Connecticut' when substring(region,length(region)-1,length(region)+1) ilike 'DE' then 'Delaware' when substring(region,length(region)-1,length(region)+1) ilike 'FL' then 'Florida' when substring(region,length(region)-1,length(region)+1) ilike 'GA' then 'Georgia' when substring(region,length(region)-1,length(region)+1) ilike 'HI' then 'Hawaii' when substring(region,length(region)-1,length(region)+1) ilike 'ID' then 'Idaho' when substring(region,length(region)-1,length(region)+1) ilike 'IL' then 'Illinois' when substring(region,length(region)-1,length(region)+1) ilike 'IN' then 'Indiana' when substring(region,length(region)-1,length(region)+1) ilike 'IA' then 'Iowa' when substring(region,length(region)-1,length(region)+1) ilike 'KS' then 'Kansas' when substring(region,length(region)-1,length(region)+1) ilike 'KY' then 'Kentucky' when substring(region,length(region)-1,length(region)+1) ilike 'LA' then 'Louisiana' when substring(region,length(region)-1,length(region)+1) ilike 'ME' then 'Maine' when substring(region,length(region)-1,length(region)+1) ilike 'MD' then 'Maryland' when substring(region,length(region)-1,length(region)+1) ilike 'MA' then 'Massachusetts' when substring(region,length(region)-1,length(region)+1) ilike 'MI' then 'Michigan' when substring(region,length(region)-1,length(region)+1) ilike 'MN' then 'Minnesota' when substring(region,length(region)-1,length(region)+1) ilike 'MS' then 'Mississippi' when substring(region,length(region)-1,length(region)+1) ilike 'MO' then 'Missouri' when substring(region,length(region)-1,length(region)+1) ilike 'MT' then 'Montana' when substring(region,length(region)-1,length(region)+1) ilike 'NE' then 'Nebraska' when substring(region,length(region)-1,length(region)+1) ilike 'NV' then 'Nevada' when substring(region,length(region)-1,length(region)+1) ilike 'NH' then 'New Hampshire' when substring(region,length(region)-1,length(region)+1) ilike 'NJ' then 'New Jersey' when substring(region,length(region)-1,length(region)+1) ilike 'NM' then 'New Mexico' when substring(region,length(region)-1,length(region)+1) ilike 'NY' then 'New York' when substring(region,length(region)-1,length(region)+1) ilike 'NC' then 'North Carolina' when substring(region,length(region)-1,length(region)+1) ilike 'ND' then 'North Dakota' when substring(region,length(region)-1,length(region)+1) ilike 'OH' then 'Ohio' when substring(region,length(region)-1,length(region)+1) ilike 'OK' then 'Oklahoma' when substring(region,length(region)-1,length(region)+1) ilike 'OR' then 'Oregon' when substring(region,length(region)-1,length(region)+1) ilike 'PA' then 'Pennsylvania' when substring(region,length(region)-1,length(region)+1) ilike 'RI' then 'Rhode Island' when substring(region,length(region)-1,length(region)+1) ilike 'SC' then 'South Carolina' when substring(region,length(region)-1,length(region)+1) ilike 'SD' then 'South Dakota' when substring(region,length(region)-1,length(region)+1) ilike 'TN' then 'Tennessee' when substring(region,length(region)-1,length(region)+1) ilike 'TX' then 'Texas' when substring(region,length(region)-1,length(region)+1) ilike 'UT' then 'Utah' when substring(region,length(region)-1,length(region)+1) ilike 'VT' then 'Vermont' when substring(region,length(region)-1,length(region)+1) ilike 'VA' then 'Virginia' when substring(region,length(region)-1,length(region)+1) ilike 'WA' then 'Washington' when substring(region,length(region)-1,length(region)+1) ilike 'WV' then 'West Virginia' when substring(region,length(region)-1,length(region)+1) ilike 'WI' then 'Wisconsin' when substring(region,length(region)-1,length(region)+1) ilike 'WY' then 'Wyoming' when substring(region,length(region)-1,length(region)+1) ilike 'AS' then 'American Samoa' when substring(region,length(region)-1,length(region)+1) ilike 'DC' then 'District of Columbia' when substring(region,length(region)-1,length(region)+1) ilike 'FM' then 'Federated States of Micronesia' when substring(region,length(region)-1,length(region)+1) ilike 'GU' then 'Guam' when substring(region,length(region)-1,length(region)+1) ilike 'MH' then 'Marshall Islands' when substring(region,length(region)-1,length(region)+1) ilike 'MP' then 'Northern Mariana Islands' when substring(region,length(region)-1,length(region)+1) ilike 'PW' then 'Palau' when substring(region,length(region)-1,length(region)+1) ilike 'PR' then 'Puerto Rico' when substring(region,length(region)-1,length(region)+1) ilike 'VI' then 'Virgin Islands'
		end
	when strpos(region, ',') > 0 then left(region, strpos(region, ',') - 1)
	else region
end
where region is not null
and country ilike 'usa'
;
--group by country,region
--order by region;