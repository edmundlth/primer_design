"""
Utilities functions that will be called by the primer_design.py program
Author : Edmund Lau (elau1@student.unimelb.edu.au)
------------------------------------------------------------------------
Description:
Provide the utilitis which enable
- visualisation
- file handling
------------------------------------------------------------------------


"""
from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_GC, Tm_Wallace
import math
import itertools

def generate_output_pair(template, tile, f_r_lens):
    start = tile[0]
    end = tile[0] + tile[1]
    f_len, r_len = f_r_lens
    relevant_template = str(template[start - f_len: end + r_len])
    f_primer = relevant_template[:f_len]
    r_primer = rev_complement(relevant_template[-r_len:])
    output = []
    for primer in [f_primer, r_primer]:
        info = map(str,(primer, raw_tm(primer), raw_entropy(primer), 
              raw_hairpin(primer), raw_gc(primer),
              raw_gc_clamp(primer)))
        output.append(info)
    return tuple(output)


def visualise_tile(template, tile, f_r_lens):
    template = str(template)
    start = tile[0]
    end = tile[0] + tile[1]
    f_len, r_len = f_r_lens
    relevant_template = template[start - f_len:end + r_len]
    complement_template = complement(relevant_template)
    f_primer ="5'>" +  relevant_template[:f_len] + ">3'"
    f_bond =' '*3 + '|' * f_len
    r_primer = ' ' * (end - start + f_len)+"3'<"\
                 + complement(relevant_template[- r_len:])+"<5'"
    r_bond = ' ' * (end - start + f_len+3) + '|' * r_len
    relevant_template = "5'>" + relevant_template + ">3'"
    complement_template = "3'<" + complement_template + "<5'"
    print(relevant_template)
    print(r_bond)
    print(r_primer)
    print(f_primer)
    print(f_bond)
    print(complement_template)
    return (relevant_template, r_bond, r_primer,
            f_primer, f_bond, complement_template)




def bed_coords(file_name):
    bed_file = open(file_name)
    start_len_coords = []
    for line in bed_file:
        coords = line.split('\t')[:3]
        start = int(coords[1])
        length = int(coords[2]) - start
        start_len_coords.append((coords[0],start,length))
    bed_file.close()
    return start_len_coords



def write_primer(primer_info, file_handle):
    """ Write the information of the chosen primer
    to the output file. 
    """
    file_handle.write('\t'.join(primer_info) + '\n')
    file_handle.flush()

def raw_entropy(primer):
    primer = primer.upper()
    length = float(len(primer))
    prob = [primer.count(base)/length for base in set(primer)]
    entropy = - sum( p * math.log(p) for p in prob)
    return entropy

def raw_tm(primer):
    return Tm_NN(primer)

def raw_hairpin(primer):
    num_match = 0
    length = len(primer)
    for i in range(1,length):
        top = primer[:i][::-1]
        bottom = primer[i:]
        new_num_match = weighted_num_complement(top,bottom)
        if new_num_match > num_match:
            num_match = new_num_match
    return num_match

def raw_self_dimer(primer):
    length1 = len(seq1)
    length2 = len(seq2)
    total_length = len(seq1) + len(seq2)
    score = 0
    for i in range(1,total_length):
        top = seq1[-i:total_length -i]
        bottom = seq2[i- total_length: i]
        new_score = num_complement(top,bottom)
        score = new_score if new_score > score else score
    return score

def raw_gc(primer):
    primer = primer.upper()
    count = float(primer.count('G') + primer.count('C'))
    return count/ len(primer)

def raw_gc_clamp(primer):
    primer = primer.upper()
    if primer[-1] == 'G' or primer[-1] == 'C':
        return 1
    return 0
    

def raw_run(primer):
    run =  max(len(list(group[1]))-1
               for group in itertools.groupby(primer.upper()))
    return run
############## Utils for scoring #######################
def weighted_num_complement(seq1,seq2):
    weighted_num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            b1 = b1.upper()
            if b1 == 'G' or b1 == 'C':
                weighted_num_complement +=1.5
            else:
                weighted_num_complement +=1
    return weighted_num_complement

def num_complement(seq1,seq2):
    num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            num_complement +=1
    return num_complement


            
            

def is_complement(base1,base2):
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    if base1 in complement_dic:
        return complement_dic[base1] == base2
    else:
        return False

def get_primer(template, tile, primer_length, which_primer = 'f'):
    which_primer = which_primer.lower()
    if which_primer == 'f':
        f_coord = tile[0]
        return template[f_coord - primer_length : f_coord]
    elif which_primer == 'r':
        r_coord = tile[0] + tile[1]
        return template[r_coord : r_coord + primer_length]
    else:
        raise ValueError('Input f or r to specify forward or reverse primer')

def complement(seq):
    seq = seq.upper()
    complement = {'A':'T','G':'C','T':'A','C':'G'}
    return ''.join(complement[base] for base in seq)

def rev_complement(seq):
    seq = seq.upper()
    complement = {'A':'T','G':'C','C':'G','T':'A'}
    return ''.join(complement[base] for base in seq)[::-1]
