from lxml import etree
from lxml.etree import ElementTree

from xml_helper import text_at_node

from loguru import logger
from locations import get_local_folder_for, FileType
from os.path import sep
import os

arguments = {
    'sars_cov_1': (f'.{sep}generated{sep}NCBI_sars_cov_1{sep}samples{sep}30271926.xml', f'.{sep}annotations{sep}new_sars_cov_1.tsv'),
    'sars_cov_2': (f'.{sep}generated{sep}New NCBI SARS-Cov-2{sep}samples{sep}sc2-refseq.xml', f'.{sep}annotations{sep}new_ncbi_sars_cov_2.tsv'),
    'dengue_virus_1': (f'.{sep}generated{sep}Dengue Virus 1{sep}samples{sep}NC_001477.1.xml', f'.{sep}annotations{sep}dengue_virus_1.tsv'),
    'dengue_virus_2': (f'.{sep}generated{sep}Dengue Virus 2{sep}samples{sep}NC_001474.2.xml', f'.{sep}annotations{sep}dengue_virus_2.tsv'),
    'dengue_virus_3': (f'.{sep}generated{sep}Dengue Virus 3{sep}samples{sep}NC_001475.2.xml', f'.{sep}annotations{sep}dengue_virus_3.tsv'),
    'dengue_virus_4': (f'.{sep}generated{sep}Dengue Virus 4{sep}samples{sep}NC_002640.1.xml', f'.{sep}annotations{sep}dengue_virus_4.tsv'),
    'mers': (f'.{sep}generated{sep}MERS-CoV{sep}samples{sep}NC_019843.3.xml', f'.{sep}annotations{sep}mers.tsv'),
    'betacoronavirus_england_1': (f'.{sep}generated{sep}Betacoronavirus England 1{sep}samples{sep}NC_038294.1.xml', f'.{sep}annotations{sep}betacoronavirus_england_1.tsv'),
    'zaire_ebolavirus': (f'.{sep}generated{sep}Zaire ebolavirus{sep}samples{sep}NC_002549.1.xml', f'.{sep}annotations{sep}zaire_ebolavirus.tsv'),
    'sudan_ebolavirus': (f'.{sep}generated{sep}Sudan ebolavirus{sep}samples{sep}NC_006432.1.xml', f'.{sep}annotations{sep}sudan_ebolavirus.tsv'),
    'reston_ebolavirus': (f'.{sep}generated{sep}Reston ebolavirus{sep}samples{sep}NC_004161.1.xml', f'.{sep}annotations{sep}reston_ebolavirus.tsv'),
    'bundibugyo_ebolavirus': (f'.{sep}generated{sep}Bundibugyo ebolavirus{sep}samples{sep}NC_014373.1.xml', f'.{sep}annotations{sep}bundibugyo_ebolavirus.tsv'),
    'bombali_ebolavirus': (f'.{sep}generated{sep}Bombali ebolavirus{sep}samples{sep}NC_039345.1.xml', f'.{sep}annotations{sep}bombali_ebolavirus.tsv'),
    'tai_forest_ebolavirus': (f'.{sep}generated{sep}Tai Forest ebolavirus{sep}samples{sep}NC_014372.1.xml', f'.{sep}annotations{sep}tai_forest_ebolavirus.tsv'),
}


def generate_annotation_file(from_reference_sammple_file_path: str, destination_file_path: str):
    def concat_intervals(e):
        intervals = e.xpath(".//INSDInterval")
        intervals2 = []
        for i in intervals:
            start = text_at_node(i, './/INSDInterval_from')
            stop = text_at_node(i, './/INSDInterval_to')
            intervals2.append((start, stop))

        if intervals:
            return ';'.join([','.join(pair) for pair in intervals2])
        else:
            return None, None

    sample_xml: ElementTree = etree.parse(from_reference_sammple_file_path, parser=etree.XMLParser(remove_blank_text=True))
    features_nodes = sample_xml.xpath('/INSDSet/INSDSeq/INSDSeq_feature-table/INSDFeature')
    annotations = []
    for a_feature in features_nodes:
        try:
            chromosmes = a_feature.xpath('INSDFeature_intervals/INSDInterval/INSDInterval_accession')
            chromosome_name = text_at_node(
                chromosmes[0],
                '.',
                mandatory=True)
            for c in chromosmes:
                if text_at_node(c, '.', mandatory=True) != chromosome_name:
                    logger.warning(f'different chromosome names found while generating {destination_file_path}')
            start_stop_string = concat_intervals(a_feature)
            feature_type = text_at_node(
                a_feature,
                './/INSDFeature_key') or '.'
            if feature_type == 'source':
                continue
            feature_type = feature_type.replace('mat_peptide', 'mature_protein_region')
            gene_name = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "gene"]/INSDQualifier_value',
                False) or '.'
            product = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "product"]/INSDQualifier_value',
                False) or '.'
            translation = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "translation"]/INSDQualifier_value',
                False)
            peptide = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "peptide"]/INSDQualifier_value',
                False)
            amino_acid_sequence = translation or peptide or '.'
            protein_id = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "protein_id"]/INSDQualifier_value',
                False) or '.'

            # for a in annotations:
            #     if a[3] == start_stop_string and a[4] == gene_name and a[2] == feature_type:
            #         raise AssertionError()
            annotations.append(
                (chromosome_name, 'RefSeq', feature_type, start_stop_string, gene_name, product, protein_id, amino_acid_sequence)
            )
        except AssertionError as e:
            pass



    annotations_copy = []
    removed = []


    try:
        for i in range(len(annotations)):
            do_not_add =  False
            a = annotations[i]
            a_start = a[3][:a[3].index(',')]
            a_stop = a[3][a[3].rindex(',')+1:]
            for j in range(i+1,len(annotations)):
                a2 = annotations[j]
                a2_start = a2[3][:a2[3].index(',')]
                a2_stop = a2[3][a2[3].rindex(',')+1:]
                # print(f"a: {a[3]} ->  {a_start} - {a_stop} vs a2: {a2[3]} -> {a2_start} - {a2_stop}")
                if a_start == a2_start and a_stop == a2_stop and a[4] == a2[4]:
                    if a[2] == 'gene':
                        do_not_add = True
                        removed.append(a)
                    elif a[5] == a2[5] and a[7] == a2[7]:
                        do_not_add = True
                        removed.append(a)
            if not do_not_add:
                annotations_copy.append(a)
    except ValueError:
        print('ANNOTATIONS')
        for a in annotations:
            print(*a, sep='\t', end='\n')
        print('\n\n')

        print('ANNOTATIONS COPY')
        for a in annotations_copy:
            print(*a, sep='\t', end='\n')
        print('\n\n')

        print('TO REMOVE')
        print(*removed)
    except IndexError:
        logger.exception(f"len annotations: {len(annotations)}, i: {i}, j: {j}")

    sorted(annotations_copy, key=lambda tup: tup[3])

    # for a in annotations_copy:
    #     print(*a, sep='\t', end='\n')
    # print('\n\n')
    # for a in removed:
    #     print(*a, sep='\t', end='\n')

    with open(destination_file_path, mode='w') as ann_file:
        for a in annotations_copy:
            line = '\t'.join(a)
            ann_file.write(line+'\n')


def create_snpeff_folders_for_viruses():
    folders = [
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}betacoronavirus_england_1',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}bombali_ebolavirus',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}bundibugyo_ebolavirus',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}dengue_virus_1',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}dengue_virus_2',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}dengue_virus_3',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}dengue_virus_4',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}mers',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}reston_ebolavirus',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}sudan_ebolavirus',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}tai_forest_ebolavirus',
        f'.{sep}tmp_snpeff{sep}snpEff{sep}data{sep}zaire_ebolavirus'
    ]
    for path in folders:
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)


create_snpeff_folders_for_viruses()
get_local_folder_for('Dengue Virus 1', FileType.SequenceOrSampleData)
get_local_folder_for('Dengue Virus 2', FileType.SequenceOrSampleData)
get_local_folder_for('Dengue Virus 3', FileType.SequenceOrSampleData)
get_local_folder_for('Dengue Virus 4', FileType.SequenceOrSampleData)
get_local_folder_for('MERS-CoV', FileType.SequenceOrSampleData)
get_local_folder_for('Betacoronavirus England 1', FileType.SequenceOrSampleData)
get_local_folder_for('Zaire ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('Sudan ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('Reston ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('Bundibugyo ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('Bombali ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('Tai Forest ebolavirus', FileType.SequenceOrSampleData)
get_local_folder_for('New NCBI SARS-Cov-2', FileType.SequenceOrSampleData)
get_local_folder_for('NCBI_sars_cov_1', FileType.SequenceOrSampleData)
# TODO you have to download the xml files listed in arguments and put them in the correct path manually in order for this to work
for x in arguments.keys():
    generate_annotation_file(*arguments[x])

