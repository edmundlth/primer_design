"""
Primer design tool for Hiplex MPS-PCR
Author: Edmund Lau (elau1@student.unimelb.edu.au)
----------------------------------------------------------------
Description:
At the moment, the program takes in a template DNA, and the coordinates
specifiying the region of interest in the template which is to be amplified,
and try to find the best tiling which produce the best scored primer
combination to perform multiplex PCR.
**All terms and conventions are defined below.
----------------------------------------------------------------


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

"""

from utils import score, visualise
from time import time
from random import randint

import sys
sys.setrecursionlimit(10000)

# Currently primer length is fixed. 
PRIMER_LENGTH = 20

# GLOBAL VARIABLES
# MEMO is used in dynamic programming's memoization
# MEMO_CALLS and CALLS keep track of the times memo is looked up and
# the number of recursive calls.
# This should implies : len(MEMO) + MEMO_CALLS = CALLS
MEMO = {}
MEMO_CALLS = 0
CALLS = 0

def design(template, tile_sizes, pos, length):
    """
    This is the primary function of this program. It takes in a
    DNA template a specifed region of interest and tile it so that:
    1) Each tiles have sizes as specified in tile_sizes
    2) A tile has to overlap the region of interest by at least one base.
       **implying left and right boundary of the leftmost and rightmost
         tiles may exceed the region of interest.
    3) The resultant primer combination that will amplify the whole region
       has maximal score.

    template : the string of DNA where the region of interest is embedded
    tile_sizes : a list of integer specify legal tile size
    pos : a integer specify the 0-based index of the starting position of
          the region of interest
    length : an integer specify the length of the region of interest.
    """
    # Initialisation
    # best_score is initialise to -inf since current scoring allows -ve score
    best_combination = []
    best_score = -10000000
    # loop through each possible tile sizes and
    # all its possible overlap configuration
    # then choose the best corresponding "suffix" - the rest of the the seq
    # return the combination where
    # score = first_tile_score + best_suffix_score      is a maximum.
    for size in tile_sizes:
        for overlap in range(1, size +1):
            prefix_pos = pos - (size - overlap)
            suffix_pos = pos + overlap
            suffix_length = length - overlap
            
            first_tile = (prefix_pos, size)
            first_tile_score = score(template, first_tile, PRIMER_LENGTH)

            best_suffix_design = suffix_design(template,
                                               tile_sizes,
                                               suffix_pos,
                                               suffix_length)
            
            new_score = first_tile_score + best_suffix_design[1]
            if new_score > best_score:
                best_combination = [prefix_pos , size] +\
                                   best_suffix_design[0][1:]
                best_score = new_score
    print("Best score = %f" %best_score)
    return best_combination
                                               

def suffix_design(template, tile_sizes, suffix_pos, length):
    """
    Recursively tile the template starting from suffix_pos,
    with appropriate tile sizes so that the resultant primer combination
    has maximal score.
    Uses memoization (Dynamic programming) to make this algorithm linear
    in the variable length, len(tile_sizes)
    """
    # keeping track of how many calls are made to this function
    global CALLS, MEMO_CALLS
    CALLS +=1
    
    if length <= 0: #non-positive length not tilable.
        return ([],0)
    
    # Memoization
    if (suffix_pos, length) in MEMO:
        MEMO_CALLS +=1
        return MEMO[(suffix_pos,length)]

    # initialisation : best_score intialise to -inf since allow -ve score
    best_tiling = [suffix_pos]
    best_score = -10000
    for size in tile_sizes:
        prefix_score = score(template, (suffix_pos, size), PRIMER_LENGTH)
        if size >= length: # base case
            new_tiling = [suffix_pos,size]
            new_score = prefix_score
        else:
            # recursion: solve same problem for the suffix of suffix
            rest = suffix_design(template,
                                 tile_sizes,
                                 suffix_pos + size,
                                 length - size,)
            new_tiling = [suffix_pos,size] + rest[0][1:]
            new_score = prefix_score + rest[1]
        if new_score > best_score:
            best_tiling = new_tiling
            best_score = new_score
    result = (best_tiling, best_score)
    
    # memoize the result
    MEMO[ (suffix_pos, length) ] = result
    return result


def main():
    print("With PRIMER_LENGTH = %i and region size ~ 1000\n"
          "This will take ~15 seconds\n\n"%PRIMER_LENGTH)
    bases = ['A','T','G','C']
    
##    template1 = 'AAATGCACGAAAAATCGTGGCGTTTTATTTGTGCTAGTCGTGCGTGAAAATTCGTCCC'
##    template_length1 = len(template1) - 11
##    tile_sizes1 = [8,9,10]
##    
##    template2 = "AATGCGTCAGTTGAC"
##    template_length2 = len(template2) -6
##    tile_sizes2 = [3,4,5]
    
    tile_sizes3 = [100,101,102,103,104,105,106,107,108,109,110]
    template_length3 = 3000
    template3 = ''
    for b in [bases[randint(0,3)] for i in range(template_length3)]:
        template3 +=b
    template_length3 -= max(tile_sizes3)

    template = template3
    tile_sizes = tile_sizes3
    template_length = template_length3
    pos = max(tile_sizes) + PRIMER_LENGTH
    
    before = time()
    tiling = design(template, tile_sizes, pos, template_length)
    after = time()

    visualise(template,tiling, pos, pos + template_length, PRIMER_LENGTH )
    print("The corresponding tiling in the\n"
          "format [start_pos, tile sizes] is: \n\n",
          tiling, '\n')
    print("CALLS = %i "%CALLS)
    print("MEMO call = %i"%MEMO_CALLS)
    print("MEMO length = number of new CALLS = %i" %len(MEMO))
    print("Time_taken = %f" %(after - before))
    input("\n\n\nPress any key to quit")


if __name__ == "__main__":
    main()

