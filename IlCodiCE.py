import numpy as np

def print_matrix(m):
    righe = m.shape[0]
    colonne = m.shape[1]
    for r in range(righe):
        print("\t".join([str(m[r,c]) for c in range(colonne)]))

def add_variant(sequence_id, pos, length, original,  mutated, variant_type):
#     print(sequence_id, pos, length, original,  mutated, variant_type)
    return [sequence_id, pos, length, original,  mutated, variant_type]

def alignment(sequence_id = 1, seq1 = "attaaaggtttataccttcccaggtaaca", seq2 = "acccdccdcdccaacaca"):
    variants=[]
    indel_score = -1
    match_score = 2
    mismatch_score = -1

    limit = int(1.5*(abs(len(seq1) - len(seq2)) + 1))
    limit = max(limit, 100)
    limit = min(limit, len(seq1), len(seq2))

    use_both = 0
    use_seq1 = 1
    use_seq2 = 2
    empty = 3

    seq1 = seq1
    seq2 = seq2



    scores = np.zeros((len(seq2)+1,len(seq1)+1))
    scores.fill(-1000000)
    traces=np.zeros((len(seq2)+1,len(seq1)+1),dtype=object)

    #print("Matrices created")

    #for i in range(len(seq2)+1):
    #    for j in range(len(seq1)+1):
    #        traces[i,j] = []
    traces.fill(empty)

    #print("Trace initialized")

    for i in range(len(seq1)+1):
        scores[0,i] = i*indel_score
        if (i > 0):
            traces[0,i] = use_seq1

    for i in range(len(seq2)+1):
        scores[i,0] = i*indel_score
        if (i > 0):
            traces[i,0] = use_seq2

    #print("Scores initialized")

    #print_matrix(scores)

    for rs  in range(len(seq2)):
    #    if rs % 100 == 0:
    #        print(rs)
        for cs in range(max(0, rs-limit), min(len(seq1), rs+limit)):


            rm = rs+1
            cm = cs+1

            score_use_seq1 = scores[rm, cm - 1] + indel_score
            score_use_seq2 = scores[rm - 1, cm] + indel_score
            score_use_both = scores[rm - 1, cm - 1] + (match_score if seq1[cs]==seq2[rs] else mismatch_score)

            if score_use_both >= max(score_use_seq1, score_use_seq2):
                scores[rm, cm] = score_use_both
                traces[rm,cm] = use_both
            if score_use_seq1 >= max(score_use_seq2, score_use_both):
                scores[rm, cm] = score_use_seq1
                traces[rm,cm] = use_seq1
            if score_use_seq2 >= max(score_use_seq1, score_use_both):
                scores[rm, cm] = score_use_seq2
                traces[rm,cm] = use_seq2


    #print_matrix(scores)

    seq1_aligned = ""
    seq2_aligned = ""

    cm = len(seq1)
    rm = len(seq2)

    while not(cm == 0 and rm == 0):
        action = traces[rm, cm]
        if action == use_both:
            if seq1[cm-1] == seq2[rm-1]:
                seq1_aligned += seq1[cm-1]
                seq2_aligned += seq2[rm-1]
            else:
                seq1_aligned += seq1[cm - 1].upper()
                seq2_aligned += seq2[rm - 1].upper()
            cm -= 1
            rm -= 1
        elif action == use_seq1:
            seq1_aligned += seq1[cm - 1]
            seq2_aligned += "-"
            cm -= 1
        else:
            seq1_aligned += "-"
            seq2_aligned += seq2[rm - 1]
            rm -=1

    seq1_aligned = seq1_aligned[::-1]
    seq2_aligned = seq2_aligned[::-1]

    #print(len(seq1_aligned))
    #print(len(seq2_aligned))

    #for i in range(0, (len(seq1_aligned)//50)+1):
    #    print(seq1_aligned[i*50:(i+1)*50])
    #    print(seq2_aligned[i * 50:(i + 1) * 50])
    #    print()

    ref_positions = np.zeros(len(seq1_aligned), dtype=int)
    pos = 0
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] != '-':
            pos +=1
        ref_positions[i] = pos

    #print(ref_positions)
    ins_open = False
    ins_len = 0
    ins_pos = None
    ins_seq = ""
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] == '-':
            ins_open = True
            ins_len += 1
            ins_pos = ref_positions[i]
            ins_seq += seq2_aligned[i]
        else:
            if ins_open:
                v = add_variant(sequence_id, ins_pos, ins_len, "-"*ins_len, ins_seq, "INS")
                variants.append(v)
                
                ins_open = False
                ins_len = 0
                ins_pos = None
                ins_seq = ""
    if ins_open:
        v = add_variant(sequence_id, ins_pos, ins_len, "-"*ins_len, ins_seq, "INS")
        variants.append(v)

    del_open = False
    del_len = 0
    del_pos = None
    del_seq = ""
    for i in range(len(seq1_aligned)):
        if seq2_aligned[i] == '-':
            if not del_open:
                del_pos = ref_positions[i]
            del_open = True
            del_len += 1
            del_seq += seq1_aligned[i]
        else:
            if del_open:
                v = add_variant(sequence_id, del_pos, del_len, del_seq, "-"*del_len, "DEL")
                variants.append(v)
                
                del_open = False
                del_len = 0
                del_pos = None
                del_seq = ""

    if del_open:
        v = add_variant(sequence_id, del_pos, del_len, del_seq, "-"*del_len, "DEL")
        variants.append(v)


    mut_open = False
    mut_len = 0
    mut_pos = None
    mut_seq_original = ""
    mut_seq_mutated = ""
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] != '-' and seq2_aligned[i] != '-' and seq1_aligned[i] != seq2_aligned[i]:
            if not mut_open:
                mut_pos = ref_positions[i]
            mut_open = True
            mut_len += 1
            mut_seq_original += seq1_aligned[i]
            mut_seq_mutated += seq2_aligned[i]
        else:
            if mut_open:
                v = add_variant(sequence_id, mut_pos, mut_len, mut_seq_original.lower(), mut_seq_mutated.lower(), "SUB")
                variants.append(v)
                
                mut_open = False
                mut_len = 0
                mut_pos = None
                mut_seq_original = ""
                mut_seq_mutated = ""

    if mut_open:
        v = add_variant(sequence_id, mut_pos, mut_len, mut_seq_original.lower(), mut_seq_mutated.lower(), "SUB")
        variants.append(v)
    return variants
alignment()