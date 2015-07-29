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
from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_Wallace, Tm_GC
import math
import itertools
import random






######################## Scoring utils ##################################



def score_primer_linear(primer,target_tm, which_primer = 'f'):
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




def score_primer_sum(primer,target_tm, which_primer = 'f'):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    
    Tm = tm_score(primer,target_tm)
    run = run_score(primer)
    gc = gc_score(primer)
    hairpin =  hairpin_score(primer)
    entropy = entropy_score(primer)
    return Tm + run + gc + hairpin + entropy

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

def gc_clamp(primer, which_primer = 'f'):
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






############## Utils for scoring #######################
def weighted_num_complement(seq1,seq2):
    weighted_num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            b1 = b1.upper()
            if b1 == 'G' or b1 == 'C':
                weighted_num_complement +=2
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
        print('total_score: %s\n'%score_primer(primer,target_tm))


######## Depreciated functions #######################
def score1(template,tile,primer_length):
    """
    DEPRECIATED !!!!
    The scoring function which return the primer score of the tile.
    
    NOTE that it is on a function of template, tile and primer_length only.
    Those are "local" variable which is independent of the choice of other
    primer pairs in the ultimate combination. (a primer pair amplifying
    a tile is taken into consideration however).
    """
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




def score_tile(template, tile, f_len, r_len,target_tm):
    """ Return the score of the tile with the specified f and r
primer length. The score of the tile is simply the sum of the score
of the forward and reverse primer. Consider primer pair dimer only """
    f_primer = get_primer(template, tile, f_len, which_primer = 'f')
    r_primer = get_primer(template, tile, r_len, which_primer = 'r')
    score =   score_primer(f_primer,target_tm)\
            + score_primer(r_primer,target_tm)\
            - dimer_score(f_primer, r_primer)
    return score

def get_primer_pair(template,tile, primer_length):
    """
    template : a string of DNA
    tile : a 2-tuple (start_pos, length) specifying the tile of interest
           on the template.
    primer_length : integer
    
    Return a 2-tuple of string (f_primer, r_primer)
    This function simple return the sequences flanking the tile with
    length = primer_length
    """
    f_coord = tile[0]
    r_coord = tile[1]+1
    f_primer = template[f_coord- primer_length : f_coord]
    r_primer = template[r_coord : r_coord + primer_length + 1]
    return (f_primer, r_primer)


