"""
        Hi-Plex Primer design Tool
Authors: Bernard J. Pope, 
         Daniel J. Park,
         Tu Nguyen,
         Edmund Lau
Email: elau1@student.unimelb.edu.au
###########################################################################

One line description:
    This program design primers for HiPlex PCR reaction

Description:
    This program takes a set of DNA coordinates specifying the regions in
    the DNA that is to be amplified, a list command line arguments (see below)
    specifying various parameters and output a set of primers that will
    amplify those regions in a multiplex PCR reaction.

    The DNA coordinates are specified in BED files format which should
    be in the working directory. A reference genome contained in a
    directory called fasta containing the sequence of each chromosome
    in separate fasta files has to be provided in the working directory.

    The outputs will be recored into 3 tsv files: (using default names)
    - primer_out.txt - recording primer names, sequence, length and
                       scores with various raw data about the primer
    - primer_design.log - logging information containing the
                          commandline inputs
    - auxiliary_primer_out.txt - contain auxiliary informations
                                 about the search process including
                                 run time, tiles size mean, std_deviation,
                                 overlap mean, std_deviation

    This program uses bottom up dynamic programming (see description below)
    to sift through the search space of all possible primer set
    which will completely tile all regions specified. We allow tiles
    to overlap and primers to vary in length up to the corresponding
    user specified constrains.


Conventions:
    - Indexing are 0-based (similar to BED format and Python string)
    - output files are tab separated (.tsv)
    - Variable names ending with 's' implies a list or tuples of objects


Definitions:
    coords :: (String, Int, Int)
    - (chromosome, start pos, end pos) 0-based index
    - eg. ('chr1', 50, 100) is the downstream sequence at chromosome 1
      with 0 based index 50,51,52,....,99

    reference :: String or Bio.Seq object
    - The reference sequence in which out region of interest or
      primers are embedded in.
    - In this program, it's the whole chromosome.

    region :: (Int, Int)
    - (start pos, end pos) 0-based index, open end points
    - this is part of the reference which is intended to be amplified

    primer :: (Int, Int, String)
    - (3' position, length, direction)
    - direction is taken with respect to the sense DNA strand.
      'f' = forward = 5' -> 3' = left to right
      'r' = reverse = 3' -> 5' = right to left
    - only applicable in the context where the reference is known
    - 0-based indexing.
    - eg. (50,5, 'f') is the forward primer found by taking
      the sequence with 0-based index
                   46,47,48,49,50
      (50, 5, 'r') is the reverse primer found by taking
      the sequence with 0-based index
                   50,51,52,53,54

    primer_set :: [primer]
    - a set of primers which, tile the entire region

    sequence :: String
    - the actual sequence of nucleic acid alphabets

    tile :: (Int, Int)
    - (Start pos, length), 0-based index relative to the reference 
    - eg. (10, 5) is the tile starting at 10 with length 5, ie
      it cover the segment with index 10,11,12,13,14
    - associated with a tile is the forward and revserse primers
      (f_primer and r_primer) which, in a PCR reaction will
      ideally amplify the tile.
    - note: max_tile and min_tile aren't exactly tiles, they
            are integers representing tile size

    overlap :: Int
    - The amount intersection (shared positions) between 2
      adjacent tiles.
    - eg. The tiles (10,5) and (13,4) has overlap=2


    user_inputs :: argparse.Namespace object
    - a record of all user commandline inputs



Brief Description of the Bottom Up Dynamic Programming:
    Pseudocode:

    memo = some recording data structure

    for each possible position:
        for each possible tile choice at this position:
            score of this pos = score of this choice +
                                   score all previous choices
                                   compatible to this choice
                                   (this is the entry of memo
                                   at position - tile size)

        record best score in the memo along with the choices
        made which enable a post-process to back-trace
        the optimum choices

    # not considering primer length variation and overlaps yet
    # see comments in code body for adaptation
    # now that the memo is prepared
    search through memo:
       the best scored final position will be our solution


Details of Dynamic Programming:
    Legal tiling pattern:
    - All tiles are of specified sizes and has at least 1 base
      intersection with the region
    - The union of the tiles must at least cover the region
    - Tiles can overlap but the number of overlapping bases is at most
      'allowed_overlap'.

There are 2 memoization in this Dp, 
    - pos_memo 
         -> memoize, as a list, 
            the best (score, tile size, overlap, f_primer, r_primer)
         -> f_primer, r_primer are unique representation of the 
            forward and reverse primer respectively in the context of a
            reference sequence, (see 'primer' in definitions section).
         -> Number of entry of the memo 
             = number of bases that can possibly tiled by legal tiling
             = region_length + 2 * (max_tile -1)
         -> Each entry correspond to a position and will be updated
            during the Dp so that the entry represent the best
            tile to choose at this position given all the previous
            position's scores while obeying legal tiling rule.

    - primer_memo
          -> memoize the processed primers and their scores as a 
             dictionary with primers as keys and tuple of scores 
             and primers' data as value.

Example:
    Let say we are at position P in the middle of Dp, with the
    following reference with * representing the region and - repesenting
    the flanking slack space.
    The numbers below each base represent the best score as determined
    by previous loops.
    
     - - - - - * * * * * * * * * - - - - -
     1 3 3 6 8 5 7 3 6 8 9 2 1 P

    Now, let tile sizes given be [2,3] and allowed_overlap=1
    (all 'tile_score' entries below are hypothetical)
    
    Check 1: current best_score = 0
         t_size = 2
         overlap = 0
         tile_score = 2
         This would give total_score = 2 + previous_pos_score
                                     = 2 + 2
                                     = 4
                (greater than best_score thus update best_score)
    
     - - - - - * * * * * * * * * - - - - -
     1 3 3 6 8 5 7 3 6 8 9 2 1 P
                             ^=^      # tile of size 2 overlaping p
                      ...==^
    Check 2: current best_score = 4
        t_size = 2
        overlap = 1
        tile_score = 4
        Total_score = 4 + 1  # The immediate previous position since
                    = 5      # previous tile to overlap it by 1 base.
        (update best_score)


     - - - - - * * * * * * * * * - - - - -
     1 3 3 6 8 5 7 3 6 8 9 2 1 P
                             ^=^
                      .....==^
    
    Check 3: current best_score = 5
        t_size = 3
        overlap = 0
        tile_score = 1
        total_score = 1 + 9 = 10 (update best_score)
    
     - - - - - * * * * * * * * * - - - - -
     1 3 3 6 8 5 7 3 6 8 9 2 1 P
                           ^===^
                    ...==^

    Check 4: current best_score = 10
        t_size = 3
        overlap = 1
        tile_score = 5
        total_score = 3 + 2   (no update needed)


     - - - - - * * * * * * * * * - - - - -
     1 3 3 6 8 5 7 3 6 8 9 2 1 P
                           ^===^
                     ....==^


    Finish checking ...
        best_score = 10
        memoize best_score = 10
                best_tile_size = 3
                best_overlap = 0
                best_f and best_f are determined by best_primer_in_tile()

    go to next position...
"""
from argparse import ArgumentParser
import logging
import sys
import os
import time
import csv
from multiprocessing import Pool

from score import Score
from primer_search import Dp_search, regions_ref_seqs_generator, write_primer, write_aux
###############################################################################

DEFAULT_LOG_FILE = 'primer_design.log'
DEFAULT_OUTFILE = 'primer_out.tsv'
DEFAULT_PRIMER_LENGTH = 20
DEFAULT_PRIMER_LENGTH_VAR = 8
DEFAULT_ALLOWED_OVERLAP = 10
DEFAULT_TM_FUNC = 'Tm_NN'
DEFAULT_SCORE_FUNC = 1      # this is the p-value in the Lp-norm
DEFAULT_SCORE_WEIGHTS = [1, 1, 1, 1, 1, 1]    # default to no weighting
DEFAULT_ADJUST_SCORE = 1.0   # this is the exponent in Score/numtile^exponent
DEFAULT_TARGET_TM = 64.0
DEFAULT_TM_UNDERACHIEVE = 1.0
DEFAULT_SALTCORR = 3
DEFAULT_CONC = [50, 0, 0, 0, 0, 25, 25]
DEFAULT_GC_WEIGHT = 1.5
DEFAULT_AUXFILE = 'auxiliary_primer_out.tsv'


def parse_args():
    'Parse command line arguments for the program'
    parser = ArgumentParser(description='Hi-Plex primer design tool')
    parser.add_argument('--log', metavar="LOG_FILE", type=str,
                        default=DEFAULT_LOG_FILE,
                        help='Log file. Default to %s'%DEFAULT_LOG_FILE)
    parser.add_argument('--bed', metavar="BED_FILE", type=str,
                        required=True,
                        help='''Path to the BED file specifying all the 
                        coordinates of the regions of interest ''')
##    parser.add_argument('--merge', metavar='BOOL', type=bool,
#                        default=True,
#                        help='''Boolean value that specify if user wants
#                        to get warning about regions that are too close
#                        to chromosome's boundary and
#                        merged regions in the bedfiles that are too 
#                        close to each other and tile those regios instead.
#                        ''')
    parser.add_argument('--fa', metavar='FASTA_FILE', type=str,
                        required=True,
                        help='''The path to the reference fasta files''')
    parser.add_argument('--outfile', metavar="OUTPUT_FILE", type=str,
                        default=DEFAULT_OUTFILE,
                        help=''' A string specifying the name of the output
                        file''')
    parser.add_argument('--sense_heel', metavar='SENSE_HEEL', type=str,
                        required=True,
                        help='''A string of nucleic acid specifying the
                        sense heel strand''')
    parser.add_argument('--antisense_heel', metavar='ANTISENSE_HEEL', type=str,
                        required=True,
                        help='''A string of nucleic acid specifying the
                        antisense heel strand''')
    parser.add_argument('--tiles', metavar="MAX_MIN_TILE_SIZE", type=int,
                        required=True, nargs=2,
                        help=''' A pair of integers specifying the maximum
                        and the minimum tile sizes inclusive''')
    parser.add_argument('--primer_length', metavar='PRIMER_LENGTH', type=int,
                        default=DEFAULT_PRIMER_LENGTH,
                        help='''An integer specifying the optimal primer length
                        Defaulted to %s'''%DEFAULT_PRIMER_LENGTH)
    parser.add_argument('--primer_length_var', metavar='PRIMER_LENGTH_VARIATION',
                        type=int, default=DEFAULT_PRIMER_LENGTH_VAR,
                        help='''An integer specifying the amount of variation
                        from the optimal primer length.
                        Defaulted to %s
                        Eg: optimal length of 20 and variation of 5 gives
                            primer length ranging from 15 to 25 inclusive'''
                        %DEFAULT_PRIMER_LENGTH_VAR)
    parser.add_argument('--allowed_overlap', metavar="ALLOWED_TILE_OVERLAP",
                        type=int, default=DEFAULT_ALLOWED_OVERLAP,
                        help='''An integer specifying the allowed overlaping
                        between successive tiles.
                        Defaulted to %s'''%DEFAULT_ALLOWED_OVERLAP)
    parser.add_argument('--score_func', metavar='P_NORM_VAL',
                        type=int, default=DEFAULT_SCORE_FUNC,
                        choice=[1,2],
                        help='''A number specifying the p value in 
                        the Lp-norm used. In particular, p=1 reduces to
                        a linear sum. 
                        Default to %s'''%DEFAULT_SCORE_FUNC)
    parser.add_argument('--score_weights', metavar='WEIGHT',
                        type=float, nargs='*',
                        default=DEFAULT_SCORE_WEIGHTS,
                        help='''A list of weights for scoring, 
                        the order of the list is given by,
                        score =   a0 * tm 
                                + a1 * entropy
                                + a2 * hairpin
                                + a3 * gc_content
                                + a4 * gc_clamp
                                + a5 * run
                        Default to:
                        [a0, a1, a2, a3, a4, a5] = %s'''%DEFAULT_SCORE_WEIGHTS)
    parser.add_argument('--adjust_score', metavar='EXPONENT', type=float,
                        default=DEFAULT_ADJUST_SCORE,
                        help='''The exponent, e, in the formula,
                        adjusted_score = score / (tile_count)^e.
                        This is correction is added so that the algorithm
                        will not bias towards high tile count since 
                        'score' is additive.
                        Default to %s'''%DEFAULT_ADJUST_SCORE)
    parser.add_argument('--tm_func', metavar='TM_FUNCTION', type=str,
                        default=DEFAULT_TM_FUNC,
                        choices=['Tm_NN', 'Tm_GC', 'Tm_Wallace'],
                        help='''A string specifying the name of the
                        melting temperature prediction algorithm used.
                        Defaulted to %s'''%DEFAULT_TM_FUNC)
    parser.add_argument('--target_tm', metavar="TARGET_TM", type=float,
                        default=DEFAULT_TARGET_TM,
                        help='''A floating point number specifying the
                        target melting point of the reaction mixture
                        in degree Celsius.
                        Defaulted to %s degC'''%DEFAULT_TARGET_TM)
    parser.add_argument('--tm_underachieve', metavar="TM_UNDERACHIEVE_PENALTY",
                        type=float, default=DEFAULT_TM_UNDERACHIEVE,
                        help='''A floating point number specifying the
                        penalty weight given to an underachieving primer
                        in term of tm, ie has lower melting temperature than
                        target_tm. Lower melting temperature should be score
                        harsher than higher ones, there for the value should
                        be > 1.''')
    parser.add_argument('--NN_table', metavar='NN_TABLE',
                        type=int, default=3,
                        choices=range(1,5),
                        help='''An integer specifying the nearest neighbour
                        thermodynamic table to be used.\n
                        1 -- Breslauer et al. (1986)\n
                        2 -- Sugimoto et al. (1996)\n
                        3 -- Allawi and SantaLucia (1997)\n
                        4 -- SantaLucia & Hicks (2004)\n
                        ''')
    parser.add_argument('--saltcorr', metavar="SALTCORR", type=int,
                        default=DEFAULT_SALTCORR,
                        choices=range(0,8),
                        help=''' An integer from 0 to 7 inclusive indicating
                        the saltcorrection method to be used during
                        melting temperature prediction.
                        0 - no salt corrections.
                        To see individual method, refer to
                        BioPython Bio.SeqUtils.MeltingTemp module
                        Default to %s.'''%DEFAULT_SALTCORR)
    parser.add_argument('--conc', metavar="CONC", type=float,
                        nargs=7, default=DEFAULT_CONC,
                        help=''' 7 floating point numbers which give the
                        the concentration of the following chemical species:
                            - Na (in mM)\n
                            - K  (in mM)\n
                            - Mg (in mM)\n
                            - Tris-HCl (in mM)\n
                            - dNTPs (in mM)\n
                            - dnac1 (in nM) note the nM. This is the
                                            concentration of the more abundant
                                            DNA species (likely to be primer)\n
                            - dnac2 (in nM) Concentration of the less abundant
                                            DNA species.\n

                            Default to %s
                            '''%DEFAULT_CONC)
    parser.add_argument('--gc_weight', metavar="GC_WEIGHT", type=float,
                        default=DEFAULT_GC_WEIGHT,
                        help='''A floating point number >1 specifying the
                        weight given to G-C binding over A-T binding during
                        the scoring of hairpins.
                        Default to %s'''%DEFAULT_GC_WEIGHT)
    parser.add_argument('--auxfile', metavar="AUXILIARY_FILE", type=str,
                        default=DEFAULT_AUXFILE,
                        help='''A string specifying the name of an auxiliary
                        output files that will record some additional data
                        and statistics''')
    return parser.parse_args()

def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOG_FILE'''
    logging.basicConfig(filename=log,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


def main():
    ''' Main function: Entry point of primer_design program'''
    #parse commandline inputs
    user_inputs = parse_args()
    start_log(user_inputs.log)
    primer_design(user_inputs)


def primer_design(user_inputs):
    """
    The function that takes in the all command line parameter 
    and execute the primer design process.
    """
    # handle the user input regarding scoring
    # the score function that will be used
    score_func = Score(user_inputs).score_func
    with open(user_inputs.outfile, 'w') as outfile:
        with open(user_inputs.auxfile, 'w') as auxfile:
            aux_header = header = ['chrom', 'region_start', 'region_end', 
                                   'region_length', 'time_taken', 
                                   'mean_tile', 'tile_count',
                                   'std_deviation_tile', 'mean_overlap',
                                   'num_primers_scored']
            outfile_header = ['name', 'start', 'end', 'sequence', 'length',
                              'score', 'tm', 'entropy', 'hairpin',
                              'gc', 'gc_clamp', 'run']
            out_writer = csv.writer(outfile,delimiter='\t')
            aux_writer = csv.writer(auxfile, delimiter='\t')
            out_writer.writerow(outfile_header)
            aux_writer.writerow(aux_header)


            # Now we loop through all specified regions and their reference sequence
            # search for optimal primer set (prepare the memos),
            # pick the optimal primer set,
            # write primers and their datas to output files
            bedfilename = user_inputs.bed
            correction_exponent = user_inputs.adjust_score
            # use parallel processing to design primers for independent regions
            process = lambda region_ref_seq : process_region(region_ref_seq, 
                                                            user_inputs, 
                                                            score_func)
            #multi = Pool(2)
            # This will become multiprocess.map when the pickling has been handled
            processed = map(process, regions_ref_seqs_generator(bedfilename,
                                                                      user_inputs))
            map(lambda searcher: write_primer(out_writer, searcher, correction_exponent),
                processed)
            map(lambda searcher: write_aux(aux_writer, searcher), processed)
            

def process_region(region_ref_seq, user_inputs, score_func):
    searcher = Dp_search(user_inputs, region_ref_seq, score_func)
    #write_primer(out_writer, searcher)
    
    coords = region_ref_seq[0]
    #write_aux(aux_writer, searcher)
    sys.stderr.write(
                     "Time taken for %s %s %s with length %s : %s\n"
                     %(coords[0],coords[1],coords[2], 
                       searcher.region_length, searcher.time_taken)
                    )
    return searcher


if __name__ == '__main__':
    main()
