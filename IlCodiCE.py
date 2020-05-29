import numpy as np
from Bio import pairwise2
import os

def add_variant(sequence_id, pos_ref, pos_seq,  length, original,  mutated, variant_type, reference, sequence):
    #return [sequence_id, pos_ref, pos_seq, length, original,  mutated, variant_type]
    if variant_type == "INS":
        return "\t".join(
            map(str, ["NC_045512", 
                      max(1,pos_ref), 
                      ".", 
                      reference[max(1,pos_ref)-1], 
                      sequence[0:length+1] if (pos_ref == 0) else sequence[pos_seq-2:pos_seq-1+length], 
                      ".", 
                      ",".join(map(str, [variant_type,pos_seq,length, original, mutated])) ]))
    elif variant_type == "DEL":
        return "\t".join(
            map(str, ["NC_045512", 
                      pos_ref, 
                      ".", 
                      reference[0:length+1] if (pos_seq == 0) else reference[pos_ref-2:pos_ref-1+length], 
                      sequence[max(1,pos_seq)-1], 
                      ".", 
                      ",".join(map(str, [variant_type,pos_seq,length, original, mutated])) ]))
        
    else:
        return "\t".join(
            map(str, ["NC_045512", 
                      pos_ref, 
                      ".", 
                      original, 
                      mutated, 
                      ".", 
                      ",".join(map(str, [variant_type,pos_seq,length, original, mutated])) ]))


def create_aligner_to_reference(reference):
    
    def inner_fun(sequence, sequence_id = 1):
        variants = []
        
        alignments = pairwise2.align.globalms(reference, sequence, 2, -1, -1, -.5)
        ref_aligned = alignments[0][0]
        seq_aligned = alignments[0][1]
        
        ref_positions = np.zeros(len(seq_aligned), dtype=int)
        pos = 0
        for i in range(len(ref_aligned)):
            if ref_aligned[i] != '-':
                pos +=1
            ref_positions[i] = pos
            
        seq_positions = np.zeros(len(seq_aligned), dtype=int)
        pos = 0
        for i in range(len(seq_aligned)):
            if seq_aligned[i] != '-':
                pos += 1
            seq_positions[i] = pos
           
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
                    v = add_variant(sequence_id, ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference, sequence)
                    variants.append(v)

                    ins_open = False
                    ins_len = 0
                    ins_pos = None
                    ins_pos_seq = None
                    ins_seq = ""
        if ins_open:
            v = add_variant(sequence_id, ins_pos, ins_pos_seq, ins_len, "-" * ins_len, ins_seq, "INS", reference, sequence)
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
                    v = add_variant(sequence_id, del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference, sequence)
                    variants.append(v)

                    del_open = False
                    del_len = 0
                    del_pos = None
                    del_pos_seq = None
                    del_seq = ""

        if del_open:
            v = add_variant(sequence_id, del_pos, del_pos_seq, del_len, del_seq, "-" * del_len, "DEL", reference, sequence)
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
                    v = add_variant(sequence_id, mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB", reference, sequence)
                    variants.append(v)

                    mut_open = False
                    mut_len = 0
                    mut_pos = None
                    mut_pos_seq = None
                    mut_seq_original = ""
                    mut_seq_mutated = ""

        if mut_open:
            v = add_variant(sequence_id, mut_pos, mut_pos_seq, mut_len, mut_seq_original, mut_seq_mutated, "SUB", reference, sequence)
            variants.append(v)
            
        
        with open("./tmp_snpeff/{}.vcf".format(sequence_id), "w") as f:
            for m in variants:
                f.write(m+'\n')
            
        os.system("java -jar ./tmp_snpeff/snpEff/snpEff.jar  covid  ./tmp_snpeff/MN000998.vcf > ./tmp_snpeff/output_{}.vcf".format(sequence_id, sequence_id))
        
        
        with open("./tmp_snpeff/output_{}.vcf".format(sequence_id)) as f:
            annotated_variants = [line for line in f if not line.startswith("#")]
        
        return annotated_variants
        
    
    return inner_fun
