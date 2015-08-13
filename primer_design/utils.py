import sys



DNA_NN3 = { 'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8), 
            'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 
            'sym': (0, -1.4),
            'AA/TT': (-7.9, -22.2), 'AT/TA': (-7.2, -20.4), 
            'TA/AT': (-7.2, -21.3), 'CA/GT': (-8.5, -22.7), 
            'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0),
            'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 
            'GC/CG': (-9.8, -24.4), 'GG/CC': (-8.0, -19.9) }


DNA_IMM1 = { 'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 
             'CG/GT': (-4.1, -11.7), 'CT/GG': (-2.8, -8.0), 
             'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3), 
             'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 
             'TG/AT': (-0.1, -1.7), 'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3), 
             'AA/TG': (-0.6, -2.3), 'AG/TA': (-0.7, -2.3), 
             'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2), 
             'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7), 
             'TG/AA': (3.0, 7.4), 'AC/TT': (0.7, 0.2), 
             'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5), 
             'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5), 
             'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7), 
             'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7), 
             'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8), 
             'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2), 
             'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8), 
             'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2), 
             'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5), 
             'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6), 
             'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4), 
             'TT/AT': (0.2, -1.5), 'AI/TC': (-8.9, -25.5), 
             'TI/AC': (-5.9, -17.4), 'AC/TI': (-8.8, -25.4), 
             'TC/AI': (-4.9, -13.9), 'CI/GC': (-5.4, -13.7), 'GI/CC': (-6.8, -19.1),
             'CC/GI': (-8.3, -23.8), 'GC/CI': (-5.0, -12.6), 
             'AI/TA': (-8.3, -25.0), 'TI/AA': (-3.4, -11.2), 'AA/TI': (-0.7, -2.6), 
             'TA/AI': (-1.3, -4.6), 'CI/GA': (2.6, 8.9), 'GI/CA': (-7.8, -21.1), 
             'CA/GI': (-7.0, -20.0), 'GA/CI': (-7.6, -20.2), 
             'AI/TT': (0.49, -0.7), 'TI/AT': (-6.5, -22.0), 'AT/TI': (-5.6, -18.7), 
             'TT/AI': (-0.8, -4.3), 'CI/GT': (-1.0, -2.4), 'GI/CT': (-3.5, -10.6), 
             'CT/GI': (0.1, -1.0), 'GT/CI': (-4.3, -12.1), 
             'AI/TG': (-4.9, -15.8), 'TI/AG': (-1.9, -8.5), 'AG/TI': (0.1, -1.8), 
             'TG/AI': (1.0, 1.0), 'CI/GG': (7.1, 21.3), 'GI/CG': (-1.1, -3.2), 
             'CG/GI': (5.8, 16.9), 'GG/CI': (-7.6, -22.0), 
             'AI/TI': (-3.3, -11.9), 'TI/AI': (0.1, -2.3), 'CI/GI': (1.3, 3.0), 
             'GI/CI': (-0.5, -1.3)}
DNA_TMM1 = { 'AA/TA': (-3.1, -7.8), 'TA/AA': (-2.5, -6.3), 
             'CA/GA': (-4.3, -10.7), 'GA/CA': (-8.0, -22.5), 
             'AC/TC': (-0.1, 0.5), 'TC/AC': (-0.7, -1.3), 
             'CC/GC': (-2.1, -5.1), 'GC/CC': (-3.9, -10.6), 
             'AG/TG': (-1.1, -2.1), 'TG/AG': (-1.1, -2.7), 
             'CG/GG': (-3.8, -9.5), 'GG/CG': (-0.7, -19.2), 
             'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9), 
             'CT/GT': (-6.1, -16.9), 'GT/CT': (-7.4, -21.2), 
             'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 
             'CA/GC': (-2.6, -5.9), 'CC/GA': (-2.7, -6.0), 
             'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1), 
             'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0), 
             'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 
             'CC/GT': (-3.2, -8.0), 'CT/GC': (-3.9, -10.6), 
             'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8), 
             'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2), 
             'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 
             'CA/GG': (-3.9, -9.6), 'CG/GA': (-6.0, -15.5), 
             'GA/CG': (-4.3, -11.1), ' GG/CA': (-4.6, -11.4),
             'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8), 
             'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 
             'CG/GT': (-3.8, -9.0), 'CT/GG': (-6.6, -18.7), 
             'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1), 
             'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)} 

DNA_DE1 = { 'AA/.T': (0.2, 2.3), 'AC/.G': (-6.3, -17.1),
            'AG/.C': (-3.7, -10.0), 'AT/.A': (-2.9, -7.6), 
            'CA/.T': (0.6, 3.3), 'CC/.G': (-4.4, -12.6), 
            'CG/.C': (-4.0, -11.9), 'CT/.A': (-4.1, -13.0), 
            'GA/.T': (-1.1, -1.6), 'GC/.G': (-5.1, -14.0), 
            'GG/.C': (-3.9, -10.9), 'GT/.A': (-4.2, -15.0), 
            'TA/.T': (-6.9, -20.0), 'TC/.G': (-4.0, -10.9), 
            'TG/.C': (-4.9, -13.8), 
            'TT/.A': (-0.2, -0.5), '.A/AT': (-0.7, -0.8), 
            '.C/AG': (-2.1, -3.9), '.G/AC': (-5.9, -16.5), 
            '.T/AA': (-0.5, -1.1), '.A/CT': (4.4, 14.9), 
            '.C/CG': (-0.2, -0.1), '.G/CC': (-2.6, -7.4), 
            '.T/CA': (4.7, 14.2), '.A/GT': (-1.6, -3.6), 
            '.C/GG': (-3.9, -11.2), '.G/GC': (-3.2, -10.4), 
            '.T/GA': (-4.1, -13.1), '.A/TT': (2.9, 10.4), 
            '.C/TG': (-4.4, -13.1), '.G/TC': (-5.2, -15.0),
            '.T/TA': (-3.8, -12.6)}


def stability(seq1,seq2, temp = 60.0, shift = 0,
             nn_table = DNA_NN3,
             imm_table = DNA_IMM1,
             tmm_table = DNA_TMM1,
             de_table = DNA_DE1):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    backup_seq1 = seq1
    backup_seq2 = seq2
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    dH = 0
    dS = 1
    deltaH = 0
    deltaS = 0

    if shift or (len(seq1) != len(seq2)):
        if shift > 0:
            seq1 = '.' * shift + seq1
        if shift < 0:
            seq2 = '.' * abs(shift) + seq2
        if len(seq2) > len(seq1):
            seq1 += (len(seq2) - len(seq1)) * '.'
        elif len(seq2) < len(seq1):
            seq2 += (len(seq1) - len(seq2)) * '.'

        while seq1[:2] == '..' or seq2[:2] == '..':
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        while seq1[-2:] == '..' or seq2[-2:] == '..':
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]


        # dangling ends
        # Note: Differ from biopython in that, if there is 
        # dangling ends AND mismatch, it is scored as if the mismatch isn't ther
        if seq1[0] == '.' or seq2[0] == '.':
            left_de = seq1[:2] + '/' + seq2[:2]
            try:
                deltaH += de_table[left_de][dH]
                deltaS += de_table[left_de][dS]
            except KeyError:
                left_de = left_de[0] + complement_dic[left_de[-1]] + left_de[2:]
                deltaH += de_table[left_de][dH]
                deltaS += de_table[left_de][dS]
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        if seq1[-1] == '.' or seq2[-1] == '.':
            right_de = seq1[-2:][::-1] + '/' + seq2[-2:][::-1]
            try:
                deltaH += de_table[right_de][dH]
                deltaS += de_table[right_de][dS]
            except KeyError:
                right_de = right_de[0] + complement_dic[right_de[-1]] + right_de[2:]
                deltaH += de_table[right_de][dH]
                deltaS += de_table[right_de][dS]
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]


    # terminal mismatch
    left_tmm = seq2[:2][::-1] + '/' + seq1[:2][::-1]
    if left_tmm in tmm_table:
        deltaH += tmm_table[left_tmm][dH]
        deltaS += tmm_table[left_tmm][dS]
        seq1 = seq1[1:]
        seq2 = seq2[1:]

    right_tmm = seq1[-2:] + '/' + seq2[-2:]
    if right_tmm in tmm_table:
        deltaH += tmm_table[right_tmm][dH]
        deltaS += tmm_table[right_tmm][dS]
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]


    ends = backup_seq1[0] + backup_seq1[-1]
    AT = ends.count('A') + ends.count('T')
    GC = ends.count('G') + ends.count('C')
    deltaH += nn_table['init_A/T'][dH] * AT
    deltaS += nn_table['init_A/T'][dS] * AT
    deltaH += nn_table['init_G/C'][dH] * GC
    deltaS += nn_table['init_G/C'][dS] * GC

    for basenumber in range(len(seq1) -1):
        neighbors = seq1[basenumber:basenumber +2] + '/' + \
                    seq2[basenumber:basenumber +2]
        if neighbors in imm_table:
            deltaH += imm_table[neighbors][dH]
            deltaS += imm_table[neighbors][dS]
        elif neighbors[::-1] in imm_table:
            deltaH += imm_table[neighbors[::-1]][dH]
            deltaS += imm_table[neighbors[::-1]][dS]
        elif neighbors in nn_table:
            deltaH += nn_table[neighbors][dH]
            deltaS += nn_table[neighbors][dS]
        elif neighbors[::-1] in nn_table:
            deltaH += nn_table[neighbors[::-1]][dH]
            deltaS += nn_table[neighbors[::-1]][dS]
        else:
            print("Warning: data for %s not found"%neighbors)

    deltaG = deltaH + temp * deltaS
    return deltaG

def weighted_num_complement(seq1,seq2,gc_weight = 1.5):
    weighted_num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            b1 = b1.upper()
            if b1 == 'G' or b1 == 'C':
                weighted_num_complement += gc_weight
            else:
                weighted_num_complement +=1                                                 
    return weighted_num_complement
    

def is_complement(base1,base2):
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    if base1 in complement_dic:
        return complement_dic[base1] == base2
    else:
        return False



def handle_bedfile(bedfile):
    file = open(bedfile)
    bed_dictionary = {}
    for line in file:
        chromo, start, end = line.strip().split('\t')[:3]
        start, end = map(int, [start, end])
        if chromo in bed_dictionary:
            bed_dictionary[chromo].append(  (start, end)  )
        else:
            bed_dictionary[chromo] = [(start,end)]
    return bed_dictionary



complement = {'A':'T','G':'C','C':'G','T':'A'}
def rev_complement(seq):
    seq = seq.upper()
    return ''.join(complement[base] for base in seq)[::-1]







