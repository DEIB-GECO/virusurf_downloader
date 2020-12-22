from Bio import SeqIO
from Bio import Align, Seq, pairwise2
from Bio.Data import CodonTable
from io import StringIO
import numpy as np
import os
# from .prepared_parameters import parameters
import json
import datetime
from loguru import logger


class BlastResult:
    def __init__(self, matched_id, lenght, pident):
        self.matched_id = matched_id
        self.lenght = lenght
        self.pident = pident

    def items(self):
        return [("lenght", self.lenght), ("Percentage of identical matches", self.pident)]


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)


def is_int(n):
    try:
        int(n)
        return True
    except ValueError:
        return False


def to_int(n):
    try:
        return int(n)
    except Exception:
        return None


def is_float(n):
    try:
        int(n)
        return True
    except ValueError:
        return False


def to_float(n):
    try:
        return float(n)
    except Exception:
        return None


def is_date(v):
    try:
        datetime.datetime.strptime(v, '%Y-%m-%d')
        return True
    except ValueError:
        return False


def to_date(n):
    try:
        return str(datetime.datetime.strptime(n, '%Y-%m-%d').date())
    except Exception:
        return None


def get_metadata_schema(metadata):
    schema = []
    for name in metadata[list(metadata.keys())[0]].keys():
        metadatum_dict = {
            "name" : name,
            "forPopulationDescription": True,
            "forFiltering": True,
            "type": name if name == "lineage" else "categorical"
        }
        schema.append(metadatum_dict)

    for meta in schema:
        name = meta['name']
        if name != "lineage":
            values = [v[name] for a,v in metadata.items()]
            if all([is_int(x) for x in values if x != ""]):
                meta['type'] = 'numerical'
                for sid in metadata:
                    metadata[sid][name] = to_int(metadata[sid][name])
            elif all([is_float(x) for x in values if x != ""]):
                meta['type'] = 'numerical'
                for sid in metadata:
                    metadata[sid][name] = to_float(metadata[sid][name])
            elif all([is_date(x) for x in values if x != ""]):
                meta['type'] = 'date'
                for sid in metadata:
                    metadata[sid][name] = to_date(metadata[sid][name])
            else:
                meta['type'] = 'categorical'
                for sid in metadata:
                    if metadata[sid][name] == '':
                        metadata[sid][name] = None

    return schema


def extract_nuc_mut_for_json(mut):
    return [mut["start_original"],
            mut["sequence_original"],
            mut['sequence_alternative'],
            mut["variant_type"],
            mut["annotations"]]


def filter_ann_and_variants(annotations_w_aa_variants):
    """
    Transforms SUBs and DELs so that they're all of length 1
    Removes
    - substitutions whose alternative sequence is X (aligner error)
    """
    new_annotations_w_aa_variants = []
    for gene_name, product, protein_id, feature_type, start, stop, nuc_seq, amino_acid_seq, aa_variants in annotations_w_aa_variants:
        # filter variants
        new_aa_variants = []
        for gene, protein_name, protein_code, mutpos, ref, alt, mut_type in aa_variants:
            # transform variants
            if mut_type == 'DEL':
                for i in range(len(ref)):
                    new_mutpos = mutpos + i if mutpos is not None else None
                    new_aa_variants.append(
                        (gene, protein_name, protein_code, new_mutpos, ref[i], '-', mut_type))
            elif mut_type == 'SUB':
                for i in range(len(ref)):
                    if alt[i] == 'X':
                        continue
                    else:
                        new_mutpos = mutpos + i if mutpos is not None else None
                        new_aa_variants.append(
                            (gene, protein_name, protein_code, new_mutpos, ref[i], alt[i], mut_type))
            else:
                new_aa_variants.append((gene, protein_name, protein_code, mutpos, ref, alt, mut_type))
        # include annotations of type gene and only the other ones coding for amino acids
        if feature_type == 'gene' or amino_acid_seq is not None:
            new_annotations_w_aa_variants.append((
                gene_name, product, protein_id, feature_type, start, stop, nuc_seq, amino_acid_seq,
                new_aa_variants
            ))
    return new_annotations_w_aa_variants


def call_annotation_variant(annotation_file, ref_aligned, seq_aligned, ref_positions, seq_positions, sequence_id = 666):
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

                alignment_aa = pairwise2.align.globalms(annotation.aa_seq, aa_seq, 3, -1, -3, -1)

                try:
                    ref_aligned_aa = alignment_aa[0][0]
                    seq_aligned_aa = alignment_aa[0][1]
                except IndexError:
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
        elif atype == 'gene':
            list_annotations.append(
                (gene, protein, protein_id, atype, nuc_start, nuc_stop, nuc_seq, None, []))

    return list_annotations


def filter_nuc_variants(nuc_variants):
    """
    Transforms nucleotide variants of type SUB longer than 1 in multiple variants of length 1.
    Ignore insertions or substitutions of alternative sequence 'n'.
    Removes duplicated variant impacts.
    """
    new_nuc_variants = []
    for n in nuc_variants:
        seq_original = n['sequence_original']
        seq_alternative = n['sequence_alternative']
        start_original = int(n['start_original'])
        start_alternative = int(n['start_alternative'])
        variant_length = int(n['variant_length'])
        variant_type = n['variant_type']
        impacts = n['annotations']
        impacts_set = set(tuple(values) for values in impacts if not values[0].startswith('GU280'))
        # split substituions into single point mutations + ignore SUBs or INS to 'n'
        if (variant_type == 'SUB' and variant_length > 1) or seq_alternative.lower()=='n':
            for i in range(variant_length):
                if seq_alternative[i].lower()!= 'n':
                    new_nuc_variants.append({
                        'sequence_original': seq_original[i],
                        'sequence_alternative': seq_alternative[i],
                        'start_original': start_original + i,
                        'start_alternative': start_alternative + i,
                        'variant_length': 1,
                        'variant_type': variant_type,
                        'annotations': impacts_set
                })
        else:
            n['annotations'] = impacts_set
            new_nuc_variants.append(n)
    return new_nuc_variants


def parse_annotated_variants(annotated_variants):
    result = []
    for variant in annotated_variants:
        try:
            _, start_original, _, _, _, _, others, snpeff_ann = variant.split("\t")
        except ValueError:
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


class InputException(Exception):
    def __init__(self, message):
        self.msg = message


def add_variant_factory(chr_name):
    def add_variant(pos_ref, pos_seq, length, original, mutated, variant_type, reference, sequence):
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

    return add_variant


def call_nucleotide_variants(sequence_id, reference, sequence, ref_aligned, seq_aligned, ref_positions, seq_positions,
                             chr_name, snpeff_database_name):

    add_variant = add_variant_factory(chr_name)
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
                v = add_variant(ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference,
                                sequence)
                variants.append(v)

                ins_open = False
                ins_len = 0
                ins_pos = None
                ins_pos_seq = None
                ins_seq = ""
    if ins_open:
        v = add_variant(ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference, sequence)
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
                if del_pos != 1:
                    v = add_variant(del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference,
                                    sequence)
                    variants.append(v)

                del_open = False
                del_len = 0
                del_pos = None
                del_pos_seq = None
                del_seq = ""

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
                v = add_variant(mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB", reference, sequence)
                variants.append(v)

                mut_open = False
                mut_len = 0
                mut_pos = None
                mut_pos_seq = None
                mut_seq_original = ""
                mut_seq_mutated = ""


    if mut_open:
        v = add_variant(mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB", reference, sequence)
        variants.append(v)

    variant_file = "./tmp_snpeff/{}.vcf".format(sequence_id)
    with open(variant_file, "w") as f:
        for m in variants:
            f.write(m + '\n')

    if variants:
        os.system("java -jar ./tmp_snpeff/snpEff/snpEff.jar  {}  {} > ./tmp_snpeff/output_{}.vcf".format(snpeff_database_name,
                                                                                                         variant_file,
                                                                                                         sequence_id))

        try:
            with open("./tmp_snpeff/output_{}.vcf".format(sequence_id)) as f:
                annotated_variants = [line for line in f if not line.startswith("#")]
            os.remove("./tmp_snpeff/output_{}.vcf".format(sequence_id))
        except FileNotFoundError:
            annotated_variants = list()
            pass
    else:
        annotated_variants = list()

    try:
        os.remove("./tmp_snpeff/{}.vcf".format(sequence_id))
    except:
       pass
    return filter_nuc_variants(parse_annotated_variants(annotated_variants))


def choose_alignment(alignments: Align.PairwiseAlignments):
    try:
        first_alignment = str(next(alignments))[:-1]    # remove trailing "\n"
    except StopIteration as e:
        logger.error('No alignments available for this sequence')
        raise e
    ref_aligned, _, seq_aligned = first_alignment.split('\n')
    min_length_without_gaps = len(seq_aligned) - seq_aligned.count('-')
    # compare length of this alignment with next 10K alignments
    num_alignements = 1
    for i in range(10000):
        try:
            next_alignment = str(next(alignments))[:-1]
            next_ref_aligned, _, next_seq_aligned = next_alignment.split('\n')
            next_length_without_gaps = len(next_seq_aligned) - next_seq_aligned.count('-')
            if next_length_without_gaps < min_length_without_gaps:
                min_length_without_gaps = next_length_without_gaps
                ref_aligned = next_ref_aligned
                seq_aligned = next_seq_aligned
            num_alignements += 1
        except:
            break
    if num_alignements > 1000:
        logger.trace(f"More than 1000 alignments have been compared ({num_alignements})")
    return ref_aligned, seq_aligned


def sequence_aligner(sequence_id, reference, sequence, chr_name, annotation_file, snpeff_database_name):
    aligner = Align.PairwiseAligner()
    aligner.match_score = 3.0  # the documentation states we can pass the scores in the constructor of PairwiseAligner but it doesn't work
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -2.5
    aligner.extend_gap_score = -1

    ref_aligned, seq_aligned = choose_alignment(aligner.align(reference, sequence))
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

    annotated_variants = call_nucleotide_variants(sequence_id,
                                                  reference,
                                                  sequence,
                                                  ref_aligned,
                                                  seq_aligned,
                                                  ref_positions,
                                                  seq_positions,
                                                  chr_name,
                                                  snpeff_database_name
                                                  )

    annotations = filter_ann_and_variants(
        call_annotation_variant(annotation_file,
                                ref_aligned,
                                seq_aligned,
                                ref_positions,
                                seq_positions
                                )
    )

    return annotations, annotated_variants

#
# def parse_inputs(input_fasta, input_metadata):
#     fasta_sequences = SeqIO.parse(StringIO(input_fasta), 'fasta')
#     sequences = {x.id: x.seq.lower() for x in fasta_sequences}
#
#     metadata = {}
#     meta_rows = input_metadata.strip().split("\n")
#
#     header = meta_rows[0].strip().split(",")
#     for line in meta_rows[1:]:
#         s = line.strip().split(",")
#         sid = s[0]
#         seq_metadata = {a: v.strip() for a, v in list(zip(header, s))[1:]}
#         metadata[sid] = seq_metadata
#
#     if len(set(metadata.keys()).union(set(sequences.keys()))) > len(set(metadata.keys())):
#         raise InputException("Some sequences in the FASTA file do not have a corresponding entry in the metadata.")
#
#     if len(set(metadata.keys()).union(set(sequences.keys()))) > len(set(sequences.keys())):
#         raise InputException("Some metadata rows in do not have a corresponding sequence in the FASTA.")
#
#     return sequences, metadata
#
#
# def pipeline(sequences, metadata, pid, species = 'sars_cov_2'):
#     ref_fasta_file_name,\
#     annotation_file_name,\
#     chr_name,\
#     snpeff_db_name,\
#     blast_meta_file,\
#     product_json_file, \
#     blast_db_name = parameters[species]
#
#     print(f'#\n#\n#Pipeline: {"load parameters"}\n#\n#')
#
#     #read reference FASTA of the species
#     reference_sequence = SeqIO.parse(open(ref_fasta_file_name),
#                                      'fasta').__next__().seq
#     reference_sequence = reference_sequence.lower()
#     print(f'#\n#\n#Pipeline: {"loaded reference"}\n#\n#')
#
#     ## load blast metadata
#     blast_meta_dict = {}
#     with open(blast_meta_file) as f:
#         header = f.readline().strip().split("\t")
#         for line in f:
#             s = line.strip().split("\t")
#             blast_meta_dict[s[0]] = {a: v for a, v in zip(header[1:], s[1:])}
#
#     ## Call Pangolin for lineage assignement
#     pangolin_fasta = f"pangolin_tmp/pango_{pid}.fast"
#     pangolin_output = f"pango_{pid}.pan"
#     with open(pangolin_fasta, "w") as f:
#         for sid, seq in sequences.items():
#             f.write(f">{sid}\n")
#             f.write(f'{str(seq)}\n')
#     os.system(f"bash pangolin_script.sh {pangolin_fasta} pangolin_tmp {pangolin_output}")
#     with open("pangolin_tmp/"+ pangolin_output) as f:
#         f.readline()
#         for line in f:
#             sid, lineage, _, _, status, _ = tuple(line.strip().split(","))
#             if status != "passed_qc":
#                 lineage = "unknown"
#             metadata[sid]['lineage'] = lineage
#     os.remove("pangolin_tmp/"+ pangolin_output)
#
#     print(f'#\n#\n#Pipeline: {"Pangolin executed"}\n#\n#')
#
#     #call for blast
#     blast_out_file = f'{pid}.blast'
#     os.system(f'blastn -query {pangolin_fasta}  \
#     -db blast_db/{blast_db_name} \
#     -num_alignments 20 \
#     -num_threads 5 \
#     -outfmt "7" \
#     -out {blast_out_file}')
#
#     blast_matching_sids = {}
#     with open(blast_out_file) as f:
#         for line in f:
#             if not line.startswith("#"):
#                 s = line.strip().split("\t")
#                 query_sid = s[0]
#                 matching_sid = s[1]
#                 pident = float(s[2])
#                 length = int(s[3])
#                 bres = BlastResult(matching_sid, length, pident)
#                 blast_matching_sids[query_sid] = blast_matching_sids.get(query_sid, list())
#                 blast_matching_sids[query_sid].append(bres)
#     os.remove(blast_out_file)
#     os.remove(pangolin_fasta)
#
#     print(f'#\n#\n#Pipeline: {"Blast executed"}\n#\n#')
#
#     #read product json
#     with open(product_json_file) as json_file:
#         product_json = json.load(json_file)
#
#     #initialize json.results
#     result_json = {
#         "sequencesCount": len(sequences.keys()),
#         "chrom": chr_name,
#         "referenceSequence": str(reference_sequence),
#         "schema": get_metadata_schema(metadata),
#         "products": product_json['products']
#     }
#
#     print(f'#\n#\n#Pipeline: {"Initialize json"}\n#\n#')
#
#     annotated_variants = {}
#     annotations = {}
#
#     for sid, sequence in sequences.items():
#         print(f'#\n#\n#Pipeline: {"Analizing sequence "} {sid}\n#\n#')
#         annotated_variants[sid], annotations[sid] = sequence_aligner(sid,
#                                                                      reference_sequence,
#                                                                      sequence,
#                                                                      chr_name,
#                                                                      snpeff_db_name,
#                                                                      annotation_file_name)
#     sequences_json = {}
#     for sid in sequences.keys():
#         json_muts_nc = []
#         for mut in annotated_variants[sid]:
#             json_mut_nc = extract_nuc_mut_for_json(mut)
#             json_muts_nc.append(json_mut_nc)
#
#         json_anns = {}
#         for ann in annotations[sid]:
#             prot = ann[1]
#             aamut = [[x[3], x[4], x[5], x[6]] for x in ann[-1]]
#             json_anns[prot] = aamut
#         sequence_json = {"id": sid,
#                          "meta": metadata[sid],
#                          "closestSequences": [[mid.matched_id,
#                                                {x[0]:x[1] for x in mid.items()+list(blast_meta_dict[mid.matched_id].items())}
#                                                ] for mid in list(blast_matching_sids[sid])],
#                          "variants": {"N": {"schema": ["position",
#                                                       "from",
#                                                       "to",
#                                                       "type",
#                                                       ["effect", "putative_impact", "gene"]],
#                                             "variants": json_muts_nc},
#                                       "A": {"schema": ["position", "from", "to", "type"],
#                                             "variants": json_anns}
#                                       },
#                          "sequence": str(sequences[sid])}
#
#         sequences_json[sid] = sequence_json
#
#     result_json["sequences"] = sequences_json
#
#     output_json = {"ready": True,
#                    "result": result_json}
#
#     return output_json