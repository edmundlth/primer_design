"""
Utilities functions that will be called by the primer_design3.py program
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







def visualise(template,tiling,left,right,primer_lengths):
    """
    This function generate a visual of the position of the primers base
    on the tiling scheme inputed.
    Primer pair of a tile is simply the sequences of length = primer_length
    which flank the tile.
    Forward primers are shown above the sequence at the right position
    Reverse primers are shown below the sequence at the right position
    The region of interest is in upper case, lower case otherwise
    """
    
    template = template[:left].lower()+\
               template[left:right+1].upper()+\
               template[right+1:].lower()
    upper = ' '* (tiling[0] - primer_lengths[0][0]) +\
            template[tiling[0]- primer_lengths[0][0] : tiling[0]]
    lower = ' ' * (tiling[0]+ primer_lengths[0][1] + 1)
    template = template[:tiling[0]] + '|' + template[tiling[0]:]
    
    last_tile = tiling[0]
    for i in range(1,len(primer_lengths)):
        inc = tiling[i]
        f_primer_length = primer_lengths[i][0]
        r_primer_length = primer_lengths[i][1]
        last_tile += inc +1
        lower += ' '* (inc - r_primer_length+1) + \
                 template[last_tile: last_tile+r_primer_length]
        upper += ' '* (inc -f_primer_length +1) + \
                 template[last_tile - f_primer_length: last_tile]
        template = template[:last_tile] + '|'+template[last_tile:]
    upper = upper[:-f_primer_length]

    chunk_size = 60
    for i in range(0,len(template),chunk_size):
        print('   ',upper[i: i+chunk_size],'   ')
        print("5'>",template[i: i+chunk_size],">3'")
        print('   ',lower[i: i+chunk_size],'   ')

        


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

def rev_complement(seq):
    seq = seq.upper()
    complement = {'A':'T','G':'C','C':'G','T':'A'}
    return ''.join(complement[base] for base in seq)[::-1]


######################## Scoring #######################################

def score_tile(template, tile, f_len, r_len):
    """ Return the score of the tile with the specified f and r
primer length. The score of the tile is simply the sum of the score
of the forward and reverse primer. Consider primer pair dimer only """
    f_primer = get_primer(template, tile, f_len, which_primer = 'f')
    r_primer = get_primer(template, tile, r_len, which_primer = 'r')
    score =   score_primer(f_primer)\
            + score_primer(r_primer)\
            - dimer_score(f_primer, r_primer)
    return score

def score_primer(primer):
    """Scoring primer's local property only. Cross dimer not considered
This is a function primer sequence(string) -> Real number """
    Tm = Tm_score(primer)
    repeat = repeat_score(primer)
    gc = gc_score(primer)
    entropy = entropy_score(primer)
    hairpin =  hairpin_score(primer)
    self_dimer = 0#self_dimer_score(primer)
    
    score = Tm + gc + entropy - repeat - hairpin - entropy
    return score
    

def score1(template,tile,primer_length):
    """
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


def entropy_score(primer):
    """ H = sum( -pi * log pi) for all i) """
    primer = primer.upper()
    return sum(
        -p*math.log(p) for p in
        [primer.count(i)/float(len(primer)) for i in set(primer)])
    
def gc_score(primer):
    primer = primer.upper()
    return primer.count('G') + primer.count('C')


def repeat_score(primer):
    return sum(len(list(group[1]))-1
               for group in itertools.groupby(primer.upper())) 


def hairpin_score(primer):
    num_match = 0
    for i in range(1,len(primer)):
        top = primer[:i][::-1]
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

def Tm_score(primer):
    if primer:
        return Tm_NN(primer)
    else:
        return 0



############## Utils for scoring #######################
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



if __name__ == "__main__":
   print(tile_score("ATGCTGAAGGCTAGAGGCTAGTAGCTAGATCAGGCTGGAGACTAGCTGCGGGGGCTAGAGGACTAGC", (23,10), 15,18))



