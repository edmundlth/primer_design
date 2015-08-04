"""
This module contain all the scoring metric considered during the project.


Previous documentation is ommited because it's no longer compatible with 
current state of the program. The comment will be rethink of and added when
the version is stable again.

"""
import random

from score_utils import (tm_score, entropy_score,
                         run_score, gc_score,
                         hairpin_score, dimer_score,
                         self_dimer_score, gc_clamp_score)





######################## Scoring #######################################


def score_primer_Lp(primer,target_tm, which_primer = 'f',p = 2):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    clamp = gc_clamp_score(primer, which_primer)
    score = (Tm ** p + run ** p \
           + gc ** p + hairpin ** p \
           + entropy ** p + clamp ** p) ** (1.0/p)
    return score

def score_primer_Lp_3(primer,target_tm, which_primer = 'f',p = 2):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    if which_primer == 'f':
        primer_3prime = primer[len(primer)//2:]
    elif which_primer == 'r':
        primer_3prime = primer[:len(primer)//2]
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    entropy_3prime = entropy_score(primer_3prime)
    gc_3prime = gc_score(primer_3prime)
    run_3prime = run_score(primer_3prime)
    clamp = gc_clamp_score(primer, which_primer)

    score = (Tm ** p + run ** p \
            + gc ** p + hairpin ** p \
            + entropy ** p \
            + entropy_3prime ** p \
            + gc_3prime ** p \
            + run_3prime ** p \
            + clamp ** p) ** (1.0/p)

    return score

def score_primer_linear(primer,target_tm, which_primer = 'f'):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    clamp = gc_clamp_score(primer, which_primer)
    score =  Tm \
            + entropy\
            + clamp \
            + 0.8 * hairpin\
            + 0.8 * gc
    return score

def score_primer_linear_3(primer,target_tm, which_primer = 'f'):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    if which_primer == 'f':
        primer_3prime = primer[len(primer)//2:]
    elif which_primer == 'r':
        primer_3prime = primer[:len(primer)//2]
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    entropy_3prime = entropy_score(primer_3prime)
    gc_3prime = gc_score(primer_3prime)
    run_3prime = run_score(primer_3prime)
    clamp = gc_clamp_score(primer, which_primer)

    score = 2 * entropy_3prime \
            + 2 * gc_3prime \
            + 2 * run_3prime \
            + Tm \
            + entropy\
            + clamp \
            + 0.8 * hairpin\
            + 0.8 * gc
    return score


def score_primer_sum(primer,target_tm, which_primer = 'f'):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    clamp = gc_clamp_score(primer, which_primer)

    return Tm + run + gc + hairpin + entropy + clamp

def score_primer_sum_3(primer,target_tm, which_primer = 'f'):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    if which_primer == 'f':
        primer_3prime = primer[len(primer)//2:]
    elif which_primer == 'r':
        primer_3prime = primer[:len(primer)//2]
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    entropy_3prime = entropy_score(primer_3prime)
    gc_3prime = gc_score(primer_3prime)
    run_3prime = run_score(primer_3prime)
    clamp = gc_clamp_score(primer, which_primer)
    
    
    score = Tm + gc + entropy\
            + entropy_3prime\
            + gc_3prime + run\
            + run_3prime + hairpin + clamp
    return score
    



##############  Test #####################

if __name__ == '__main__':
    score_primer = score_primer_Lp_3
    bases = ['A','T','G','C']
    for i in range(5):
        primer = ''.join([bases[random.randint(0,3)] 
                       for i in range(23)])

        target_tm = 64.0
        print('total_score: %s\n'%score_primer(primer,target_tm))


