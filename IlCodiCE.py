import numpy as np
from Bio import pairwise2
from Bio import Seq
from Bio.Data import CodonTable
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


def create_aligner_to_reference(reference, annotation_file):
    
    table = CodonTable.ambiguous_dna_by_id[1]
    reference_annotations = []
    
    with open(annotation_file) as f:
        for line in f:
            s = line.strip().split("\t")
            name = s[0]
            start = int(s[1])
            stop = int(s[2])
            reference_annotations.append((name, start, stop))
    
    #for a in reference_annotations:
        #print(a[0], reference[a[1]])
    
    def inner_fun(sequence, sequence_id = 1):
        variants = []
        
        #name, type, start, stop, Aminoacidsequence
        annotations = []
        
        #name, pos, original, alt, type
        aa_variants = []
        
        alignments = pairwise2.align.globalms(reference, sequence, 2, -1, -1, -.5)
        ref_aligned = alignments[0][0]
        seq_aligned = alignments[0][1]
        
        #print(ref_aligned)
        #print(seq_aligned)
        
        
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
        
        paired_positions = list(zip(ref_positions, seq_positions))
        for a in reference_annotations:
            name = a[0]
            ann_type = "CDS"
            start_ref = a[1]
            stop_ref = a[2]
            ann_pos = [s for (r,s) in paired_positions if r >= start_ref and r <= stop_ref]
            seq_start = ann_pos[0]
            seq_stop = ann_pos[-1]
            dna_seq = sequence[seq_start-1:seq_stop].replace("-", "")
            dna_ref = reference[start_ref-1:stop_ref]
            aa_variants_gene = []
            
            try:
                aa_ref = Seq._translate_str(dna_ref, table, cds=True)
                aa_seq = Seq._translate_str(dna_seq, table, cds=True)
                alignment_aa = pairwise2.align.globalms(aa_ref, aa_seq, 2, -1, -1, -.5)
                
                ref_aligned_aa = alignment_aa[0][0]
                seq_aligned_aa = alignment_aa[0][1]
                                
                ref_positions_aa = np.zeros(len(seq_aligned_aa), dtype=int)
                pos = 0
                for i in range(len(ref_aligned_aa)):
                    if ref_aligned_aa[i] != '-':
                        pos +=1
                    ref_positions_aa[i] = pos
                    
                aligned_aa = list(zip(ref_positions_aa, ref_aligned_aa, seq_aligned_aa))
                mut_set = set()
                for t in aligned_aa:
                    if t[2] == "-":
                        mutpos = t[0]
                        original = t[1]
                        alternative = t[2]
                        mut_type = "DEL"
                        mut_set.add((name, mutpos, original, alternative, mut_type))
                    elif t[1] != t[2] :
                        mutpos = t[0]
                        original = [r for (p,r,s) in aligned_aa if p == mutpos if r!="-"][0]
                        alternative = "".join([s for (p,r,s) in aligned_aa if p == mutpos])
                        mut_type = "SUB"
                        mut_set.add((name, mutpos, original, alternative, mut_type))
                    elif t[1] == "-":
                        mutpos = t[0]
                        original = t[1]
                        alternative = t[2]
                        mut_type = "SUB"
                        mut_set.add((name, mutpos, original, alternative, mut_type))
                for mut in mut_set:
                    aa_variants.append(mut)
                    aa_variants_gene.append(mut)
                        
                    
                
                
            except:
                seq_start = None
                seq_stop = None
                aa_seq = None
                
            annotations.append((name, ann_type, seq_start, seq_stop, aa_seq, aa_variants_gene))
           
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
            
        
        variant_file = "./tmp_snpeff/{}.vcf".format(sequence_id)
        with open(variant_file, "w") as f:
            for m in variants:
                f.write(m+'\n')
            
        os.system("java -jar ./tmp_snpeff/snpEff/snpEff.jar  covid  {} > ./tmp_snpeff/output_{}.vcf".format(variant_file, sequence_id))
        
        
        with open("./tmp_snpeff/output_{}.vcf".format(sequence_id)) as f:
            annotated_variants = [line for line in f if not line.startswith("#")]
        os.remove("./tmp_snpeff/{}.vcf".format(sequence_id))
        

        return (annotated_variants, annotations, aa_variants)
        
    
    return inner_fun


def parse_annotated_variants(annotated_variants):
    result = []
    for variant in annotated_variants:
        _, start_original, _, _, _, _, others, _ = variant.split("\t")

        variant_type, start_alternative, variant_length, sequence_original, sequence_alternative = others.split(',')

        result.append({'sequence_original': sequence_original,
                       'sequence_alternative': sequence_alternative,
                       'start_original': start_original,
                       'start_alternative': start_alternative,
                       'variant_length': variant_length,
                       'variant_type': variant_type,
                       })
    return result
