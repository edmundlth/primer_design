"""
This module contain the functions that give a normalised score
for various criteria.

Previous comments had been deleted since they are incompatible with 
the current version anymore.

Comment will be added when the version is stable again.

"""



from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_Wallace, Tm_GC
from utils import (rev_complement, get_primer,
                   weighted_num_complement, num_complement,
                   is_complement, complement)
import math
import itertools
import random






######################## Scoring utils ##################################

def entropy_score(primer):
    """ H = sum( -pi * log pi) for all i) """
    primer = primer.upper()
    length = float(len(primer))
    prob = [primer.count(base)/length for base in set(primer)]
    entropy = - sum([p*math.log(p) for p in prob])
    # The range of entropy is 0 to log(4) since max 
    # max entropy is when all base are equaly likely to occur
    # ie max_S = -4 * (1/4*log(1/4)) = log(4)
    # we are using natural log at the moment, but normalising
    # cancel out the factors
    normalised_entropy = entropy * 100 / math.log(4.0)
    return normalised_entropy

def tm_score(primer,target_tm):
    """We assume the input is temperature in Celsius
    greater than 0 degC!!!
    From experiment with random sequences, the reverse
    complement of a sequence has the same Tm to 
    at least 5 decimal places."""
    if primer: #if primer sequence is not empty
        del_tm = Tm_NN(primer) - target_tm
    else:
        del_tm = -target_tm
    # normalise the score by scaling with range of 
    # Tm. Which is 20 - 80 degC.
    # Upper bound is given by nearest neighbor prediction
    # of a pure GC sequence with length = 30.
    # Lower bound is given by nearest neighbor prediction
    # of pure AT sequence with length = 15
    lower = 20
    upper = 80
    del_range = (upper - target_tm) - (lower - target_tm)
    normalised = del_tm / del_range
    if normalised < 0:
        # lower than target get a "doubly" lower score,
        # ie faster drop off, steeper -ve gradient
        return 100 + 120 * normalised # normalise is negative
    elif normalised > 0:
        #higher than target is "kindof" ok,
        # ie the score has less steep gradient
        return 100 - 80 * normalised
    else: # 0 is great !!
        return 100

    
def gc_score(primer):
    primer = primer.upper()
    count = float(primer.count('G') + primer.count('C'))
    normalised_gc = count/len(primer) * 100.0
    return normalised_gc


def run_score(primer):
    # run = max continuous run

    run =  max(len(list(group[1]))-1
               for group in itertools.groupby(primer.upper()))
    #normalise number of runs by the length of primer
    normalised_run = 100 * (1 - float(run) / len(primer))
    return normalised_run


def hairpin_score(primer):
    num_match = 0
    length = len(primer)
    for i in range(1,length):
        top = primer[:i][::-1]
        bottom = primer[i:]
        new_num_match = weighted_num_complement(top,bottom)
        if new_num_match > num_match:
            num_match = new_num_match
    # normalise based on the range of hairpin binding sites
    # the maximum number of binding is primer length /2
    normalised_num_match = 100 * (1 - num_match / (length/2.0))
    return normalised_num_match

def gc_clamp_score(primer, which_primer = 'f'):
    if which_primer == 'f':
        end_with = primer[-1].upper()
    elif which_primer == 'r':
        end_with = primer[0].upper()
    
    if end_with == 'G' or end_with == 'C':
        return 100
    else:
        return 0



def dimer_score(seq1,seq2):
    length1 = len(seq1)
    length2 = len(seq2)
    total_length = len(seq1) + len(seq2)
    score = 0
    for i in range(1,total_length):
        top = seq1[-i:total_length -i]
        bottom = seq2[i- total_length: i]
        new_score = num_complement(top,bottom)
        score = new_score if new_score > score else score
    # normalised base on the range of the number of dimer
    # binding sites possible, ie thelength of the shorter
    # of the 2 primer
    normalised_score = 100 * (1 - score / min(length1,length2))
    return score

def self_dimer_score(primer):
    return dimer_score(primer,primer)







############ TEST CASES #################################

if __name__ == '__main__':
    bases = ['A','T','G','C']
    for i in range(5):
        primer = ''.join([bases[random.randint(0,3)] 
                       for i in range(23)])

        target_tm = 64.0
        print('gc_score:%s'%gc_score(primer))
        print('run_score:%s'%run_score(primer))
        print('entropy_score:%s'%entropy_score(primer))
        print('hairpin_score:%s'%hairpin_score(primer))
        print('tm_score:%s'%tm_score(primer,target_tm))

