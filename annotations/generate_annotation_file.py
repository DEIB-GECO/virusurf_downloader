from lxml import etree
from lxml.etree import ElementTree

from xml_helper import text_at_node
from data_sources.ncbi_any_virus.settings import known_settings as ncbi_knonw_settings
from data_sources.ncbi_services import get_samples_accession_ids, download_or_get_ncbi_sample_as_xml
from loguru import logger
from locations import get_local_folder_for, FileType
from os.path import sep
from os import makedirs, chdir
import os
import sys


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
            # get chromosome
            chromosmes = a_feature.xpath('INSDFeature_intervals/INSDInterval/INSDInterval_accession')
            chromosome_name = text_at_node(
                chromosmes[0],
                '.',
                mandatory=True)
            # warn if more than one chromosome
            for c in chromosmes:
                if text_at_node(c, '.', mandatory=True) != chromosome_name:
                    logger.warning(f'different chromosome names found while generating {destination_file_path}')
            # interval position
            start_stop_string = concat_intervals(a_feature)
            # feature type (CDS/ UTR / etc.)
            feature_type = text_at_node(
                a_feature,
                './/INSDFeature_key') or '.'
            if feature_type == 'source':
                continue
            feature_type = feature_type.replace('mat_peptide', 'mature_protein_region')
            # gene
            gene_name = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "gene"]/INSDQualifier_value',
                False) or '.'
            gene_name = gene_name.replace('orf', 'ORF')
            # protein
            product = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "product"]/INSDQualifier_value',
                False) or '.'
            product = product.replace('orf', 'ORF')
            # AA sequence (one of translation or peptide)
            translation = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "translation"]/INSDQualifier_value',
                False)
            peptide = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "peptide"]/INSDQualifier_value',
                False)
            amino_acid_sequence = translation or peptide or '.'
            # protein ID
            protein_id = text_at_node(
                a_feature,
                './/INSDQualifier[./INSDQualifier_name/text() = "protein_id"]/INSDQualifier_value',
                False) or '.'

            annotations.append(
                (chromosome_name, 'RefSeq', feature_type, start_stop_string, gene_name, product, protein_id, amino_acid_sequence)
            )
        except AssertionError as e:
            pass

    # filter annotations (remove duplicates)
    annotations_copy = []
    removed = []

    try:
        for i in range(len(annotations)):
            # decide which annotations to consider
            do_not_add =  False
            a = annotations[i]  # pick one annotation
            # separate start_stop_string
            a_start = a[3][:a[3].index(',')]
            a_stop = a[3][a[3].rindex(',')+1:]
            # check if in the following annotations, there is one having the same start and stop coordinates
            for j in range(i+1,len(annotations)):
                a2 = annotations[j]
                a2_start = a2[3][:a2[3].index(',')]
                a2_stop = a2[3][a2[3].rindex(',')+1:]
                # print(f"a: {a[3]} ->  {a_start} - {a_stop} vs a2: {a2[3]} -> {a2_start} - {a2_stop}")

                # if same coordinates and same gene:
                #   ignore this one if the other one has same protein name and same AA sequence
                #   (this is necessary because there are identical annotations (e.g. of mature protein region) except for
                #   the protein_id.)
                if a_start == a2_start and a_stop == a2_stop and a[4] == a2[4]:
                    if a[5] == a2[5] and a[7] == a2[7]:
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


uniformed_protein_names = {
    'envelope protein': 'E (envelope protein)',
    'membrane glycoprotein': 'M (membrane glycoprotein)',
    'nucleocapsid phosphoprotein': 'N (nucleocapsid phosphoprotein)',
    "2'-O-ribose methyltransferase": "NSP16 (2'-O-ribose methyltransferase)",
    '3C-like proteinase': "NSP5 (3C-like proteinase)",
    "3'-to-5' exonuclease": "NSP14 (3'-to-5' exonuclease)",
    'endoRNAse': 'NSP15 (endoRNAse)',
    'helicase': 'NSP13 (helicase)',
    'leader protein': 'NSP1 (leader protein)',
    'nsp10': 'NSP10',
    'nsp11': 'NSP11',
    'nsp2': 'NSP2',
    'nsp3': 'NSP3',
    'nsp4': 'NSP4',
    'nsp6': 'NSP6',
    'nsp7': 'NSP7',
    'nsp8': 'NSP8',
    'nsp9': 'NSP9',
    'RNA-dependent RNA polymerase': 'NSP12 (RNA-dependent RNA polymerase)',
    'ORF3a protein': 'NS3 (ORF3a protein)',
    'ORF6 protein': 'NS6 (ORF6 protein)',
    'ORF7a protein': 'NS7a (ORF7a protein)',
    'ORF7b': 'NS7b (ORF7b)',
    'ORF8 protein': 'NS8 (ORF8 protein)',
    'surface glycoprotein': 'Spike (surface glycoprotein)'
}


def uniform_protein_names_sc2(ann_file_path: str):
    new_file_path = f"{ann_file_path}_2"
    with open(ann_file_path, mode='r') as original, open(new_file_path, mode="w+") as new:
        for original_line in original:

            chromosome_name, source, feature_type, start_stop_string, gene_name, product, protein_id, amino_acid_sequence = original_line.split('\t')
            # replace product ( == protein)
            product = uniformed_protein_names.get(product, product)
            # write
            new_line = '\t'.join([chromosome_name, source, feature_type, start_stop_string, gene_name, product, protein_id, amino_acid_sequence])
            new.write(new_line)     # trailing \n is already included in amino_acid_sequence
    os.remove(ann_file_path)
    os.rename(new_file_path, ann_file_path)


def download_refseq_of_viruses():
    for virus_key_name, import_parameter in ncbi_knonw_settings.items():
        virus_dir_name = import_parameter["generated_dir_name"]

        virus_dir_path = get_local_folder_for(virus_dir_name, FileType.SequenceOrSampleData)
        refseq_query = import_parameter["reference_sample_query"]
        refseq_sample_acc_id = get_samples_accession_ids(refseq_query)
        assert len(refseq_sample_acc_id) == 1, f'Invalid reference sequence for virus {import_parameter["log_with_name"]}'
        path_of_refseq = download_or_get_ncbi_sample_as_xml(virus_dir_path, refseq_sample_acc_id[0])
        yield path_of_refseq


def get_paths_for_annotation_files():
    for import_parameter in ncbi_knonw_settings.values():
        path_of_annotation_file = import_parameter["annotation_file_path"]
        yield path_of_annotation_file


if __name__ == '__main__':
    print("This module expects the path to the project directory as argument")
    try:
        proj_root_dir = sys.argv[1]
    except IndexError:
        print("Missing mandatory argument: path of project root dir")
        sys.exit(-1)
    chdir(proj_root_dir)
    refseq_paths = list(download_refseq_of_viruses())
    annotations_paths = list(get_paths_for_annotation_files())
    for i in range(len(refseq_paths)):
        refseq_path = refseq_paths[i]
        target_annotation_file = annotations_paths[i]

        generate_annotation_file(refseq_path, target_annotation_file)
        if ('sars' in target_annotation_file and '2' in target_annotation_file) or 'covid' in target_annotation_file:
            uniform_protein_names_sc2(target_annotation_file)