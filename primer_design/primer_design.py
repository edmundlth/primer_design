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

from Bio import SeqIO
from score import (score_tile, score_primer,
                   get_primer, dimer_score,
                   rev_complement)
from utils import visualise, visualise_tile, bed_coords
from time import time
from random import randint
from argparse import ArgumentParser
import logging

import sys
sys.setrecursionlimit(10000)


# GLOBAL VARIABLES
# MEMO is used in dynamic programming's memoization
MEMO = {}

def design(template, tile_sizes, pos, length, primer_length, length_var,target_tm):
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
    best_combination = []
    best_score = 0
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
                                                  length_var,
                                                  target_tm) 
            first_tile_score = first_tile_primers[0]
            first_primer_lengths = [first_tile_primers[1]] 

            best_suffix_design = suffix_design(template,
                                               tile_sizes,
                                               suffix_pos,
                                               suffix_length,
                                               primer_length,
                                               length_var,
                                               target_tm)
            
            new_score = first_tile_score + best_suffix_design[1]
            if new_score > best_score:
                best_combination = [prefix_pos , size] +\
                                   best_suffix_design[0][1:]
                best_score = new_score
                best_lengths = first_primer_lengths\
                        + best_suffix_design[2] 
    # print("Best score = %f" %best_score)
    return best_combination,best_lengths, best_score
                                               

def suffix_design(template, tile_sizes, suffix_pos, length,primer_length,length_var,target_tm):
    """
    Recursively tile the template starting from suffix_pos,
    with appropriate tile sizes so that the resultant primer combination
    has maximal score.
    Uses memoization (Dynamic programming) to make this algorithm linear
    in the variable length, len(tile_sizes)
    """
    
    # Memoization
    if suffix_pos in MEMO:
        return MEMO[suffix_pos]

    best_tiling = [suffix_pos]
    best_score = 0
    best_primer_lens = [(primer_length,primer_length)]
    for size in tile_sizes:
        prefix_primers = tile_primer_lens(template,
                                          (suffix_pos, size),
                                           primer_length,
                                           length_var,
                                           target_tm)
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
                                 length - size,
                                 primer_length,
                                 length_var,
                                 target_tm)
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

def tile_primer_lens(template, tile, primer_length, var,target_tm):
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
    best_f_score, best_r_score = -100000, -100000
    best_f_length, best_r_length = primer_length, primer_length
    vary_range = range(-var, var + 1)
    for vary in vary_range:
        length = primer_length + vary
        f_primer = get_primer(template, tile,
                              length, which_primer = 'f')
        r_primer = get_primer(template, tile,
                              length, which_primer = 'r')
        new_f_score = score_primer(f_primer,target_tm, which_primer = 'f')
        new_r_score = score_primer(r_primer,target_tm, which_primer = 'r')

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
    parser = ArgumentParser(description = 'Hiplex primer design tools')
    parser.add_argument('--log', metavar = "LOG_FILE", type = str,
                    default = DEFAULT_LOG_FILE, help = 'Log file. Default to %s'%DEFAULT_LOG_FILE)
    parser.add_argument('--primer_len', metavar = 'p', type = int,
                    default = 20, help = 'An integer specifying the optimal primer length, default to 20')
    parser.add_argument('--outfile', type = str, default = 'primer_out.txt', 
                    help = 'Specify the file name where the records of the program output will be written')
    parser.add_argument('--len_var', type = int, default = 5,
                    help = 'An integer specifying the optimal primer length, default to 5')
    parser.add_argument('--tiles', type = int, required = True, nargs = 2,
                    help = '''A pair of integers that specify the minimum tile size and 
                    the range of tile sizes. eg --tiles 100 5 means the tile sizes are
                    100, 101, 102, 103, 104, 105''')
    parser.add_argument('--tm', metavar = 'tm', type = float, default = 60.0, 
                    help = '''A floating point number specifying the optimal melting temperature,default to 60.0 degC''')
    parser.add_argument('--fa', type = str, nargs='*', 
                    help = '''The fasta file containing the DNA sequence where the regions of interest are embeded in''')
    parser.add_argument('--bed', type = str,
                    help = '''The BED file where the coordinates of the regions of interest are written in''')
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
    before_time = time()
    args = parse_args()
    start_log(args.log)
    coords = bed_coords(args.bed)
    target_tm = args.tm
    outfile = open(args.outfile,'w')
    primer_len = args.primer_len
    len_var = args.len_var
    tile_sizes = [args.tiles[0] + i for i in range(args.tiles[1] +1) ]
    sequences = {}
    for chromo, start_pos, length in coords:
        MEMO = {} # memo for dynamic programing. Refresh so that it's empty for every exon
        if chromo not in sequences:
            #currently we only have one sequence per fasta file
            sequences[chromo] = next(SeqIO.parse('fasta/%s.fa'%chromo, 'fasta')).seq
        template = sequences[chromo]
        tile_lengths, primer_lengths,score = design(template,
                                                    tile_sizes,
                                                    start_pos,
                                                    length,
                                                    primer_len,
                                                    len_var,
                                                    target_tm)
        pos = tile_lengths[0]
        for i in range(len(primer_lengths)):
            f_r_len = primer_lengths[i]
            length = tile_lengths[i+1]
            tile = (pos, length)
            output = visualise_tile(template, tile, f_r_len)

            outfile.write('Chromosome=%s start=%i end=%i score=%s\n'
                    %(chromo,start_pos,start_pos +length -1,score))
            outfile.write('Forward primer: %s\n'%(output[0]))
            outfile.write('Reverse primer: %s\n'%(rev_complement(output[2][-f_r_len[1]:])))

            for line in output:
                outfile.write(line+'\n')
            outfile.write('\n\n')
        outfile.flush()
    after_time = time()

    print("Time_taken = %f" %(after_time - before_time))




if __name__ == "__main__":
    main()

