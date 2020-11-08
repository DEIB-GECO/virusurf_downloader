from Bio import SeqIO
from Bio import Align
import sys
import numpy as np
import os
from pipeline_nuc_variants__annotations__aa import \
    parse_annotated_variants, \
    filter_nuc_variants, \
    call_annotation_variant, \
    filter_ann_and_variants
from data_sources.ncbi_any_virus.ncbi_importer import prepared_parameters

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
                v = add_variant(del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference,
                                sequence)
                variants.append(v)

                del_open = False
                del_len = 0
                del_pos = None
                del_pos_seq = None
                del_seq = ""

    if del_open:
        v = add_variant(del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference, sequence)
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

    return parse_annotated_variants(annotated_variants)



#def sequence_aligner(sequence_id, reference, sequence, chr_name, annotation_file, snpeff_database_name):
def sequence_aligner(sequence_id, reference, sequence, chr_name, snpeff_database_name, annotation_file):
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

    annotated_variants = filter_nuc_variants(
        call_nucleotide_variants(sequence_id,
                                 reference,
                                 sequence,
                                 ref_aligned,
                                 seq_aligned,
                                 ref_positions,
                                 seq_positions,
                                 chr_name,
                                 snpeff_database_name
                                 )
    )

    annotations = filter_ann_and_variants(
        call_annotation_variant(annotation_file,
                                ref_aligned,
                                seq_aligned,
                                ref_positions,
                                seq_positions
                                )
    )

    return annotations
    #return annotated_variants


def main():
    print("Usage: python virusviz_pipeline.py input.fasta input.csv [species (Default = new_sars_cov_2)]")

    ##parameters
    fasta_file_name = sys.argv[1]
    meta_file_name = sys.argv[2]


    try:
        species = sys.argv[3]
    except:
        species = "new_ncbi_sars_cov_2"
    reference_fasta_file_name = "./annotations/{}.fa".format(species)
    reference_sequence = SeqIO.parse(open(reference_fasta_file_name),
                                         'fasta').__next__().seq

    _,_,_,annotation_file_name,chr_name,snpeff_db_name = prepared_parameters[species]

    ##read metadata
    metadata = {}
    with open(meta_file_name) as f:
        header = f.readline().strip().split(",")

        for line in f.readlines():
            s = line.strip().split(",")
            sid = s[0]
            seq_metadata = {a: v for a, v in list(zip(header, s))[1:]}
            metadata[sid] = seq_metadata

    # read sequences
    fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
    sequences = {x.id: x.seq for x in fasta_sequences}

    pangolin_output_file = fasta_file_name + ".pan"
    os.system("bash pangolin_script.sh {} pangolin_tmp {}".format(fasta_file_name, pangolin_output_file))
    with open("pangolin_tmp/"+pangolin_output_file) as f:
        f.readline()
        for line in f:
            sid, lineage, _, _, status, _ = tuple(line.strip().split(","))
            if status != "passed_qc":
                lineage = "unknown"
            metadata[sid]['lineage'] = lineage
    os.remove("pangolin_tmp/"+pangolin_output_file)

    print("Done with lineage assignment")

    for sid, sequence in sequences.items():
        annotated_variants = sequence_aligner(sid, reference_sequence, sequence, chr_name, snpeff_db_name, annotation_file_name)

        print("Analizing sequence: {}".format(sid))




    ##todo: check that ids are same in both metadata and sequences

if __name__ == "__main__":
    main()
