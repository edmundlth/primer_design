"""

Utilities functions that will be called by the primer_design program

"""

import math
import itertools
import random



def visualise(template,tiling,left,right,primer_length):
    template = template[:left].lower()+\
               template[left:right+1].upper()+\
               template[right+1:].lower()
    upper = ' '* (tiling[0] - primer_length) +\
            template[tiling[0]- primer_length : tiling[0]]
    lower = ' ' * (tiling[0]+ primer_length + 1)
    template = template[:tiling[0]] + '|' + template[tiling[0]:]
    
    last_tile = tiling[0]
    for inc in tiling[1:]:
        last_tile += inc +1
        lower += ' '* (inc - primer_length+1) + \
                 template[last_tile: last_tile+primer_length]
        upper += ' '* (inc -primer_length +1) + \
                 template[last_tile - primer_length: last_tile]
        template = template[:last_tile] + '|'+template[last_tile:]
    upper = upper[:-primer_length]
    print(upper)
    print(template)
    print(lower)


def get_primer_pair(template,tile, primer_length):
    f_coord = tile[0]
    r_coord = tile[1]+1
    f_primer = template[f_coord- primer_length : f_coord]
    r_primer = template[r_coord : r_coord + primer_length + 1]
    return (f_primer, r_primer)





######################## Scoring #######################################
def score(template,tile,primer_length):
    pair = get_primer_pair(template,tile,primer_length)
    f_primer = pair[0]
    r_primer = pair[1]

    repeat = repeat_score(f_primer) + repeat_score(r_primer)
    gc = gc_score(f_primer) + gc_score(r_primer)
    entropy = entropy_score(f_primer) + entropy_score(r_primer)
    hairpin = hairpin_score(f_primer) + hairpin_score(r_primer)
    primer_pair_dimer = dimer_score(f_primer,r_primer)
    self_dimer = self_dimer_score(f_primer) + self_dimer_score(r_primer)
    
    score = entropy + gc - repeat - hairpin - primer_pair_dimer - self_dimer
    return score


def entropy_score(primer):
    return sum(
        -p*math.log(p) for p in
        [primer.count(i)/float(len(primer)) for i in set(primer.upper())])
    
def gc_score(primer):
    primer = primer.upper()
    return primer.count('G') + primer.count('C')


def repeat_score(primer):
    return sum(len(list(group[1]))-1
               for group in itertools.groupby(primer.upper()))


def hairpin_score(primer):
    return 0
    #this function is not fully thought of yet!!
    length = len(primer)
    num_match = 0
    for i in range(1, math.ceil(length/2.0) +1):
        top = primer[:i]
        bottom = primer[i:]
        new_num_match = num_complement(top,bottom)
        if new_num_match > num_match:
            num_match = new_num_match
    return num_match



def dimer_score(seq1,seq2):
    total_length = len(seq1) + len(seq2)
    score = 0
    for i in range(1,total_length):
        top = seq1[-i:total_length -i]
        bottom = seq2[i- total_length: i]
        new_score = num_complement(top,bottom)
        score = new_score if new_score > score else score
    return score

def self_dimer_score(primer):
    return dimer_score(primer,primer)



############## Utils for scoring ########################3
def num_complement(seq1,seq2):
    assert len(seq1) == len(seq2)
    num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            num_complement +=1
    return num_complement
            
            

def is_complement(base1,base2):
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    return complement_dic[base1] == base2
