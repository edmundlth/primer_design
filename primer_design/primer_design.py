"""
Primer design tool for Hiplex MPS-PCR
Author: Edmund Lau (elau1@student.unimelb.edu.au)
----------------------------------------------------------------
Description:
At the moment, the program takes in a template DNA, and the coordinates
specifiying the region of interest in the template which is to be
amplified, and try to find the best tiling which produce the
best scored primer combination to perform multiplex PCR.
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
-- The list of different tile sizes to choose from.
   The tile size should be  within the tolerance of
   HiPlex size selection.


primer_combination :: high level definition only
-- The set of primers (or the set of primer pairs) that will
   ultimately be added to the multiplex PCR reaction vessel.

PRIMER_LENGTH :: Int
-- the length of all primer. currently a fixed constant.

"""

from score import (score_tile, score_primer,
                   get_primer,
                   dimer_score)
from utils import visualise
from time import time
from random import randint
from argparse import ArgumentParser
import logging

import sys
sys.setrecursionlimit(10000)

#PRIMER_LENGTH = 20
#LENGTH_VAR = 2

# GLOBAL VARIABLES
# MEMO is used in dynamic programming's memoization
# MEMO_CALLS and CALLS keep track of the times memo is looked up and
# the number of recursive calls.
# This should implies : len(MEMO) + MEMO_CALLS = CALLS
MEMO = {}
MEMO_CALLS = 0
CALLS = 0

def design(template, tile_sizes, pos, length, primer_length, length_var):
    """
    This is the primary function of this program. It takes in a
    DNA template a specifed region of interest and tile it so that:
    1) Each tiles have sizes as specified in tile_sizes
    2) A tile has to overlap the region of interest by at
       least one base.
       **implying left and right boundary of the leftmost and rightmost
         tiles may exceed the region of interest.
    3) The resultant primer combination that will amplify the whole
       region has maximal score.

    template : the string of DNA where the region of interest
               is embedded
    tile_sizes : a list of integer specify legal tile size
    pos : a integer specify the 0-based index of the
          starting position of the region of interest
    length : an integer specify the length of the region of interest.
    """
    # Initialisation
    # best_score is initialise to -inf since
    # current scoring allows -ve score
    best_combination = []
    best_score = -10000000
    best_lengths = []
    # loop through each possible tile sizes and
    # all its possible overlap configuration
    # then choose the best corresponding "suffix"
    # - the rest of the the seq
    # return the combination where
    # score = first_tile_score + best_suffix_score      is a maximum.
    for size in tile_sizes:
        for overlap in range(1, size +1):
            prefix_pos = pos - (size - overlap)
            suffix_pos = pos + overlap
            suffix_length = length - overlap
            
            first_tile = (prefix_pos, size)
            first_tile_primers = tile_primer_lens(template,
                                                  first_tile,
                                                  primer_length,
                                                  length_var) 
            first_tile_score = first_tile_primers[0]
            first_primer_lengths = [first_tile_primers[1]] 

            best_suffix_design = suffix_design(template,
                                               tile_sizes,
                                               suffix_pos,
                                               suffix_length)
            
            new_score = first_tile_score + best_suffix_design[1]
            if new_score > best_score:
                best_combination = [prefix_pos , size] +\
                                   best_suffix_design[0][1:]
                best_score = new_score
                best_lengths = first_primer_lengths\
                        + best_suffix_design[2] 
    print("Best score = %f" %best_score)
    return best_combination,best_lengths
                                               

def suffix_design(template, tile_sizes, suffix_pos, length,primer_length,length_var):
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
    
    #if length <= 0: #non-positive length not tilable.
    #    return ([],0)
    
    # Memoization
    if suffix_pos in MEMO:
        MEMO_CALLS +=1
        return MEMO[suffix_pos]

    # initialisation : best_score intialise to -inf
    # since allow -ve score
    best_tiling = [suffix_pos]
    best_score = -10000
    best_primer_lens = [(PRIMER_LENGTH,PRIMER_LENGTH)]
    for size in tile_sizes:
        prefix_primers = tile_primer_lens(template,
                                          (suffix_pos, size),
                                           primer_length,
                                           length_var)
        prefix_score = prefix_primers[0]
        prefix_primer_lengths = [prefix_primers[1]] 
        if size >= length: # base case
            new_tiling = [suffix_pos,size]
            new_score = prefix_score
            new_primer_lengths = prefix_primer_lengths
        else:
            # recursion: solve same problem for the suffix of suffix
            rest = suffix_design(template,
                                 tile_sizes,
                                 suffix_pos + size,
                                 length - size,)
            new_tiling = [suffix_pos,size] + rest[0][1:]
            new_score = prefix_score + rest[1]
            new_primer_lengths = prefix_primer_lengths + rest[2]

        if new_score > best_score:
            best_tiling = new_tiling
            best_score = new_score
            best_primer_lengths = new_primer_lengths
    result = (best_tiling, best_score, best_primer_lengths)
    
    # memoize the result
    MEMO[suffix_pos] = result
    return result

def tile_primer_lens(template, tile, primer_length, var):
    '''This function returns the best scored 
    (score, (forward primer length, reverse primer length))
    of the specified tile (start pos, tile length), in the
    template (string). The function loop through all
    possible forward and reverse primer length 
    (independently, hence not considering dimer) within 
    the primer length variance specified in var to find
    the best scored combination
    The score is simply:
        score_primer(f_primer) + score_primer(r_primer)'''
    best_f_score = -100000
    best_f_length = primer_length
    best_r_score = -100000
    best_r_length = primer_length
    vary_range = range(-var, var + 1)
    for vary in vary_range:
        length = primer_length + vary
        f_primer = get_primer(template, tile,
                              length, which_primer = 'f')
        r_primer = get_primer(template, tile,
                              length, which_primer = 'r')
        new_f_score = score_primer(f_primer)
        new_r_score = score_primer(r_primer)

        if new_f_score > best_f_score:
            best_f_length = length
            best_f_score = new_f_score
        if new_r_score > best_r_score:
            best_r_length = length
            best_r_score = new_r_score
    best_score = best_f_score + best_r_score
    best_lengths = (best_f_length, best_r_length)
    return (best_score, best_lengths)




def tile_primer_lens1(template, tile, primer_length, var):
    '''Depreciated! 
    This the orginal function that returns the best scored
    primer lengths choices of a particular tile considering
    tile primer dimerisation. Because of that consideration, 
    it has be check all v^2 number of primers, where
    v = var = primer_length_variation.'''
    best_score = -100000 
    best_lengths = (primer_length, primer_length)
    vary_range = range(-var, var +1)
    for f_vary in vary_range:
        for r_vary in vary_range:
            f_len = primer_length + f_vary
            r_len = primer_length + r_vary
            new_score = score_tile(template, tile, f_len, r_len)  
            if new_score > best_score:
                best_score = new_score
                best_lengths = (f_len, r_len) 
    return (best_score, best_lengths) 

DEFAULT_LOG_FILE = "primer_design.log"

def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Design primers for Hi-Plex")
    #parser.add_argument(
    #    '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE, 
        help='Log progress in FILENAME. Defaults to {}'.format(DEFAULT_LOG_FILE))
    parser.add_argument(
        '--primer-len', metavar = 'primer_len', type = int, default = 20,
        help = 'an integer specifying optimal primer length')
    parser.add_argument(
        '--length_var', metavar = 'len_var', type = int, default = 5,
        help = 'an integer specifying primer length variance')
    parser.add_argument(
        "--tile_sizes", metavar = "tile_sizes", type = int,
        help = 'a list of integers specifying the list of tile sizes desired')
    parser.add_argument(
        '--tiles', metavar = 'tiles', type = int, nargs=2, 
    return parser.parse_args()


def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOGFILE.'''
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    # Log the command line that was used to run the program
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

def main():
    args = parse_args()
    start_log(args.log)
    bed_file = open(args.bed)
    fasta_file = open(args.fa)


    bases = ['A','T','G','C']
    
    
    tile_sizes = [100,101,102,103,104,105,106,107,108,109,110]
    template_length = 1000
    template = ''
    for b in [bases[randint(0,3)] for i in range(template_length3)]:
        template3 +=b
    template_length -= max(tile_sizes)

    pos = max(tile_sizes) + PRIMER_LENGTH
    print("region size ~ %i\n"
          %(template_length))

    before = time()
    tiling = design(template, tile_sizes, pos, template_length)
    after = time()

    visualise(template,tiling[0], pos, pos + template_length, tiling[1] )
    print("The corresponding tiling in the\n"
          "format [start_pos, tile sizes] is: \n\n",
          tiling, '\n')
    print("CALLS = %i "%CALLS)
    print("MEMO call = %i"%MEMO_CALLS)
    print("MEMO length = number of new CALLS = %i" %len(MEMO))
    print("Time_taken = %f" %(after - before))
    quiting = raw_input("\n\n\nPress return/enter to quit")
    print("Thanks :) ")


if __name__ == "__main__":
    main()

