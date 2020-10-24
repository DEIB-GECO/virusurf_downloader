from Bio import Align, Seq, pairwise2
from Bio.Data import CodonTable
import numpy as np
import os

from loguru import logger

def parse_annotated_variants(annotated_variants):
    result = []
    for variant in annotated_variants:
        try:
            _, start_original, _, _, _, _, others, snpeff_ann = variant.split("\t")
        except ValueError:
            # happened with variant = protein_coding||c.*2479_*2480insCGCC|||||2480|,ANNNN|downstream_gene_variant|MODIFIER|ORF6|GU280_gp06|transcript|GU280_gp06|...
            logger.warning(f'nuc_variant skipped because variant was: {variant}')
            continue

        annotations = []
        for ann in snpeff_ann.split(","):
            try:
                s = ann.split("|")
                annotations.append([s[1], s[2], s[3]])
            except:
                pass

        variant_type, start_alternative, variant_length, sequence_original, sequence_alternative = others.split(',')

        result.append({'sequence_original': sequence_original,
                       'sequence_alternative': sequence_alternative,
                       'start_original': start_original,
                       'start_alternative': start_alternative,
                       'variant_length': variant_length,
                       'variant_type': variant_type,
                       'annotations': annotations
                       })
    return result


def call_annotation_variant(annotation_file, ref_aligned, seq_aligned, ref_positions, seq_positions, sequence_id):
    table = CodonTable.ambiguous_dna_by_id[1]

    list_annotations = []

    class Ann:
        def __init__(self, ann_type, ann_pos, gene, protein, protein_id, aa_seq):
            self.ann_type = ann_type
            self.ann_pos = ann_pos
            self.gene = gene
            self.protein = protein
            self.protein_id = protein_id
            self.aa_seq = aa_seq

    def parse_pos(l):
        return [(int(pos.strip().split(",")[0]), int(pos.strip().split(",")[1])) for pos in l.strip().split(";")]

    ref_annotations = []
    with open(annotation_file) as f:
        for line in f:
            s = line.strip().split("\t")
            ann_type = s[2]
            ann_pos = parse_pos(s[3])
            gene = None if s[4] == "." else s[4]
            protein = None if s[5] == "." else s[5]
            protein_id = None if s[6] == "." else s[6]
            aa_seq = None if s[7] == "." else s[7]
            ref_annotations.append(Ann(ann_type, ann_pos, gene, protein, protein_id, aa_seq))

    for annotation in ref_annotations:

        gene = annotation.gene
        protein = annotation.protein
        protein_id = annotation.protein_id
        atype = annotation.ann_type
        nuc_start = annotation.ann_pos[0][0]
        nuc_stop = annotation.ann_pos[-1][1]

        nuc_seq = "".join(
            [x[1] for x in zip(ref_positions, seq_aligned) if nuc_start <= x[0] and nuc_stop >= x[0]]).replace("-", "")

        #if nuc_seq is not None and aa_seq is None:
        #    logger.warning('nuc_seq is not None and aa_seq is None (' + gene + ',' + protein + ')')


        list_mutations = []
        if annotation.ann_type == 'mature_protein_region' or annotation.ann_type == 'CDS':
            dna_ref = ''
            for (start, stop) in annotation.ann_pos:
                dna_ref += "".join([x[1] for x in zip(ref_positions, seq_aligned) if start <= x[0] and stop >= x[0]]).replace("-", "")
                            
            if len(dna_ref)%3 == 0 and len(dna_ref) > 0:
                aa_seq = Seq._translate_str(dna_ref, table, cds=False).replace("*", "")

                alignment_aa = pairwise2.align.globalms(annotation.aa_seq, aa_seq, 2, -1, -1, -.5)

                try:
                    ref_aligned_aa = alignment_aa[0][0]
                    seq_aligned_aa = alignment_aa[0][1]
                except IndexError:
                    logger.error(f'solito Index out of range in function "call_annotation_variant". This annotation will be skipped. '
                                 f'Sequence id: {sequence_id}'
                                 f'Function args were\n.'
                                 f'ref aligned: {ref_aligned}\n'
                                 f'seq aligned: {seq_aligned}\n'
                                 f'ref positions: {ref_positions}\n'
                                 f'seq positions: {seq_positions}'
                                 )
                    continue

                ref_positions_aa = np.zeros(len(seq_aligned_aa), dtype=int)
                pos = 0
                for i in range(len(ref_aligned_aa)):
                    if ref_aligned_aa[i] != '-':
                        pos += 1
                    ref_positions_aa[i] = pos

                seq_positions_aa = np.zeros(len(seq_aligned_aa), dtype=int)
                pos = 0
                for i in range(len(ref_aligned_aa)):
                    if seq_aligned_aa[i] != '-':
                        pos += 1
                    seq_positions_aa[i] = pos



                list_mutations = []

                ins_open = False
                ins_len = 0
                ins_pos = None
                ins_seq = ""
                for i in range(len(ref_aligned_aa)):
                    if ref_aligned_aa[i] == '-':
                        ins_open = True
                        ins_len += 1
                        ins_pos = ref_positions_aa[i]
                        ins_seq += seq_aligned_aa[i]
                    else:
                        if ins_open:
                            v = (gene, protein, protein_id, ins_pos, "-" * ins_len, ins_seq, "INS")
                            list_mutations.append(v)
            
                            ins_open = False
                            ins_len = 0
                            ins_pos = None
                            ins_seq = ""
                if ins_open:
                    v = (gene, protein, protein_id, ins_pos, "-" * ins_len, ins_seq, "INS")
                    list_mutations.append(v)

                del_open = False
                del_len = 0
                del_pos = None
                del_seq = ""
                for i in range(len(ref_aligned_aa)):
                    if seq_aligned_aa[i] == '-':
                        if not del_open:
                            del_pos = ref_positions_aa[i]
                        del_pos_seq = seq_positions_aa[i]
                        del_open = True
                        del_len += 1
                        del_seq += ref_aligned_aa[i]
                    else:
                        if del_open:
                            v = (gene, protein, protein_id, del_pos, del_seq, "-" * del_len, "DEL")
                            list_mutations.append(v)

                            del_open = False
                            del_len = 0
                            del_pos = None
                            del_pos_seq = None
                            del_seq = ""

                if del_open:
                    v = (gene, protein, protein_id, del_pos, del_seq, "-" * del_len, "DEL")
                    list_mutations.append(v)

                mut_open = False
                mut_len = 0
                mut_pos = None
                mut_pos_seq = None
                mut_seq_original = ""
                mut_seq_mutated = ""
                for i in range(len(ref_aligned_aa)):
                    if ref_aligned_aa[i] != '-' and seq_aligned_aa[i] != '-' and ref_aligned_aa[i] != seq_aligned_aa[i]:
                        if not mut_open:
                            mut_pos = ref_positions_aa[i]
                            mut_pos_seq = seq_positions_aa[i]
                        mut_open = True
                        mut_len += 1
                        mut_seq_original += ref_aligned_aa[i]
                        mut_seq_mutated += seq_aligned_aa[i]
                    else:
                        if mut_open:
                            v = (gene, protein, protein_id, mut_pos, mut_seq_original, mut_seq_mutated, "SUB")
                            list_mutations.append(v)

                            mut_open = False
                            mut_len = 0
                            mut_pos = None
                            mut_pos_seq = None
                            mut_seq_original = ""
                            mut_seq_mutated = ""

                if mut_open:
                    v = (gene, protein, protein_id, mut_pos, mut_seq_original, mut_seq_mutated, "SUB")
                    list_mutations.append(v)

                list_annotations.append(
                    (gene, protein, protein_id, atype, nuc_start, nuc_stop, nuc_seq, aa_seq, list_mutations))

            elif len(dna_ref) == 0:
                list_annotations.append(
                    (gene, protein, protein_id, atype, nuc_start, nuc_stop, None, None, []))

            else:
                list_annotations.append(
                    (gene, protein, protein_id, atype, nuc_start, nuc_stop, nuc_seq, None, []))

    return list_annotations


def add_variant(sequence_id, pos_ref, pos_seq, length, original, mutated, variant_type, reference, sequence, chr_name):
    # return [sequence_id, pos_ref, pos_seq, length, original,  mutated, variant_type]
    if variant_type == "INS":
        return "\t".join(
            map(str, [chr_name,
                      max(1, pos_ref),
                      ".",
                      reference[max(1, pos_ref) - 1],
                      sequence[0:length + 1] if (pos_ref == 0) else sequence[pos_seq - 2:pos_seq - 1 + length],
                      ".",
                      ",".join(map(str, [variant_type, pos_seq, length, original, mutated]))]))
    elif variant_type == "DEL":
        return "\t".join(
            map(str, [chr_name,
                      pos_ref,
                      ".",
                      reference[0:length + 1] if (pos_seq == 0) else reference[pos_ref - 2:pos_ref - 1 + length],
                      sequence[max(1, pos_seq) - 1],
                      ".",
                      ",".join(map(str, [variant_type, pos_seq, length, original, mutated]))]))

    else:
        return "\t".join(
            map(str, [chr_name,
                      pos_ref,
                      ".",
                      original,
                      mutated,
                      ".",
                      ",".join(map(str, [variant_type, pos_seq, length, original, mutated]))]))


def call_nucleotide_variants(sequence_id, reference, sequence, ref_aligned, seq_aligned, ref_positions, seq_positions,
                             chr_name, snpeff_database_name: str):

    variants = []

    ins_open = False
    ins_len = 0
    ins_pos = None
    ins_pos_seq = None
    ins_seq = ""
    for i in range(len(ref_aligned)):
        if ref_aligned[i] == '-':
            if not ins_open:
                ins_pos_seq = seq_positions[i]
            ins_open = True
            ins_len += 1
            ins_pos = ref_positions[i]
            ins_seq += seq_aligned[i]
        else:
            if ins_open:
                v = add_variant(sequence_id, ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference,
                                sequence, chr_name)
                variants.append(v)

                ins_open = False
                ins_len = 0
                ins_pos = None
                ins_pos_seq = None
                ins_seq = ""
    if ins_open:
        v = add_variant(sequence_id, ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference, sequence,
                        chr_name)
        variants.append(v)

    del_open = False
    del_len = 0
    del_pos = None
    del_pos_seq = None
    del_seq = ""
    for i in range(len(ref_aligned)):
        if seq_aligned[i] == '-':
            if not del_open:
                del_pos = ref_positions[i]
            del_pos_seq = seq_positions[i]
            del_open = True
            del_len += 1
            del_seq += ref_aligned[i]
        else:
            if del_open:
                v = add_variant(sequence_id, del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference,
                                sequence, chr_name)
                variants.append(v)

                del_open = False
                del_len = 0
                del_pos = None
                del_pos_seq = None
                del_seq = ""

    if del_open:
        v = add_variant(sequence_id, del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference, sequence,
                        chr_name)
        variants.append(v)

    mut_open = False
    mut_len = 0
    mut_pos = None
    mut_pos_seq = None
    mut_seq_original = ""
    mut_seq_mutated = ""
    for i in range(len(ref_aligned)):
        if ref_aligned[i] != '-' and seq_aligned[i] != '-' and ref_aligned[i] != seq_aligned[i]:
            if not mut_open:
                mut_pos = ref_positions[i]
                mut_pos_seq = seq_positions[i]
            mut_open = True
            mut_len += 1
            mut_seq_original += ref_aligned[i]
            mut_seq_mutated += seq_aligned[i]
        else:
            if mut_open:
                v = add_variant(sequence_id, mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB",
                                reference, sequence, chr_name)
                variants.append(v)

                mut_open = False
                mut_len = 0
                mut_pos = None
                mut_pos_seq = None
                mut_seq_original = ""
                mut_seq_mutated = ""

    if mut_open:
        v = add_variant(sequence_id, mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB", reference,
                        sequence, chr_name)
        variants.append(v)

    variant_file = "./tmp_snpeff/{}.vcf".format(sequence_id)
    with open(variant_file, "w") as f:
        for m in variants:
            f.write(m + '\n')

    if variants:
        os.system("java -jar ./tmp_snpeff/snpEff/snpEff.jar  {}  {} > ./tmp_snpeff/output_{}.vcf".format(snpeff_database_name,
                                                                                                         variant_file,
                                                                                                         sequence_id))

        os.remove("./tmp_snpeff/{}.vcf".format(sequence_id))
        try:
            with open("./tmp_snpeff/output_{}.vcf".format(sequence_id)) as f:
                annotated_variants = [line for line in f if not line.startswith("#")]
            os.remove("./tmp_snpeff/output_{}.vcf".format(sequence_id))
        except FileNotFoundError:
            annotated_variants = list()
            pass
    else:
        annotated_variants = list()

    return parse_annotated_variants(annotated_variants)


def sequence_aligner(sequence_id, reference, sequence, chr_name, annotation_file, snpeff_database_name):
    aligner = Align.PairwiseAligner()
    aligner.match_score = 3.0  # the documentation states we can pass the scores in the constructor of PairwiseAligner but it doesn't work
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    alignments = str(aligner.align(reference, sequence)[0]).strip().split('\n')

    ref_aligned = alignments[0]
    seq_aligned = alignments[2]

    ref_positions = np.zeros(len(seq_aligned), dtype=int)
    pos = 0
    for i in range(len(ref_aligned)):
        if ref_aligned[i] != '-':
            pos += 1
        ref_positions[i] = pos

    seq_positions = np.zeros(len(seq_aligned), dtype=int)
    pos = 0
    for i in range(len(seq_aligned)):
        if seq_aligned[i] != '-':
            pos += 1
        seq_positions[i] = pos

    annotated_variants = call_nucleotide_variants(sequence_id, reference, sequence, ref_aligned, seq_aligned,
                                                  ref_positions, seq_positions, chr_name, snpeff_database_name)

    annotations = call_annotation_variant(annotation_file, ref_aligned, seq_aligned, ref_positions, seq_positions, sequence_id)

    return annotations, annotated_variants




