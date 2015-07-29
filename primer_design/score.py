"""
Scoring functions that will be called by the primer_design.py program
Author : Edmund Lau (elau1@student.unimelb.edu.au)
------------------------------------------------------------------------
Description:
Provide the utilies need by the primer design algorithm.
The primary utility is the scoring function score() which itself calls
other functions which score a particular criteria:
Example (not exhausive since the file is being modified):
- entropy score
- repeats
- dimer and self-dimer
- hairpin
- gc content
------------------------------------------------------------------------


Coventions:
1) Using 0-based indexing
2) A region is specified by its STARTING POSITION and its LENGTH
   -- a_region = (start, length)
   -- this convention is applied to all tiles or template, ie.
      any contiguous sequence on the background DNA
   -- similarly, any particular tiling scheme is denoted by
      [ staring position (0-based index) , list of tile sizes ]
      e.g. [10, 100,108,105,101,100] is a tiling that starts at
           position 10 and has 5 tiles.
3) The score is normalised to a range of 0 - 100
    


Definitions: (using Haskell notation for type)

template :: String
-- The string of DNA that the region of interest is embedded in.

region :: String
-- a contigous substring of the background_DNA


coord :: (Int,Int)
-- specification of where the template is in the background DNA
   It obey all convention
   eg. (50,1000) is a template
        starting at the 51st character to the 1050th inclusive.
template_coord :: (Int,Int)
-- see coord
tile_coord :: (Int,Int)
-- see coord
primer_coord :: (Int,Int)
-- see coord

tile_sizes :: [Int]
-- The list of different tile sizes to choose from. The tile size should be
   within the tolerance of HiPlex size selection.


primer_combination :: high level definition only
-- The set of primers (or the set of primer pairs) that will ultimately be
   added to the multiplex PCR reaction vessel.

PRIMER_LENGTH :: Int
-- the length of all primer. currently a fixed constant.

    
primer_pair :: (String,String)
-- The tuple of (forward primer, reverse primer)
   at present stage, each primer is simply a substring of
   the original sequence where forward primer has 3' end right before the
   left cut of the tile and extend towards the 5' end as far as primer_length.

"""
import random

from score_utils import (tm_score, entropy_score,
                         run_score, gc_score,
                         hairpin_score, dimer_score,
                         self_dimer_score, gc_clamp)





######################## Scoring #######################################


def score_primer_Lp(primer,target_tm, which_primer = 'f',p = 2):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    clamp = gc_clamp(primer, which_primer)
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
    clamp = gc_clamp(primer, which_primer)

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
    clamp = gc_clamp(primer, which_primer)
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
    clamp = gc_clamp(primer, which_primer)

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
    clamp = gc_clamp(primer, which_primer)

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
    
    
    score = Tm + gc + entropy\
            + entropy_3prime\
            + gc_3prime + run\
            + run_3prime + hairpin 
    return score
    


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
        print('total_score: %s\n'%score_primer(primer,target_tm))


