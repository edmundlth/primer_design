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

from Bio import SeqIO

from utils import rev_complement, mean, std_deviation
from score import Score
###############################################################################

DEFAULT_LOG_FILE = 'primer_design.log'
DEFAULT_OUTFILE = 'primer_out.tsv'
DEFAULT_PRIMER_LENGTH = 20
DEFAULT_PRIMER_LENGTH_VAR = 8
DEFAULT_ALLOWED_OVERLAP = 5
DEFAULT_TM_FUNC = 'Tm_NN'
DEFAULT_SCORE_FUNC = 'score_Lp'
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
    parser.add_argument('--score_func', metavar="SCORING_FUNCTION",
                        type=str, default=DEFAULT_SCORE_FUNC,
                        choices=['score_Lp', 'score_linear'],
                        help=''' A string specifying the name of the
                        scoring function to be used.''')
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
    scoring = Score(user_inputs) # handle the user input regarding scoring
    score_func = scoring.score_func # the score function that will be used

    # parse bed file and return a a list in the form [(region_name, start, end)]
    # region_name is in the format: <chrom>_<name from 3rd column in bedfile>
    # eg: chr1_name1
    # the name might be an empty string
    bed_coords = handle_bedfile(user_inputs.bed)

    outfile_header = ['name', 'start', 'end', 'sequence', 'length',
                      'score', 'tm', 'entropy', 'hairpin',
                      'gc', 'gc_clamp', 'run']
    with open(user_inputs.outfile, 'w') as outfile:
        out_writer = csv.writer(outfile,delimiter='\t')
        out_writer.writerow(outfile_header)

        # initialise a auxiliary data dictionary that will collect data
        # from each run which will be written into an auxiliary file
        # a post process statistics is then done on it.
        auxiliary_data = {}

        # Now we loop through all specified regions,
        # search for optimal primer set (prepare the memos),
        # pick the optimal primer set,
        # write primers and their datas to output files
        for bed_data in bed_coords:
            before_time = time.time() # time each run

            region_name, start, end = bed_data
            coords = (region_name, start, end)
            #coords = (chromo, region_coord[0], region_coord[1])
            searcher = Dp_search(user_inputs, coords, score_func)
            #primer_set = pick_primer_set(searcher)
            write_primer(out_writer, searcher)
            outfile.flush()

            after_time = time.time()
            time_taken = after_time - before_time
            logging.info("%s %s %s"%coords)
            logging.info("length %i"%searcher.region_length)
            print("Time taken for %s %s %s: %s\n"%(coords[0],coords[1],coords[2],time_taken))

            # assemble aux_data
            auxiliary_data[coords] = [list(searcher.aux_data['tiles']),
                                      list(searcher.aux_data['overlap']),
                                      searcher.region_length, time_taken,
                                      len(searcher.primer_memo)]

        # The primer design process is now completed
        # write auxiliary data into auxfile
        write_aux(auxiliary_data, user_inputs.auxfile)

###################################

class Dp_search(object):
    """
    Dp_search object initialise the 'search for best primer set'
    process. It takes in all command-line arguments
       *may or may not be usefull
    , a scoring function (* see Score() class) and
    region_coords :: (chrom, start, end), eg (chr1, 0, 100).

    It initialise a class level pos_memo and primer_memo for
    Dynamic programming purpose.
    """
    def __init__(self, user_inputs, region_coords, score_func):
        self.score_primer = score_func

        self.min_tile, self.max_tile = user_inputs.tiles
        self.tile_sizes = [self.min_tile + extend
                           for extend in
                           range(self.max_tile - self.min_tile + 1)]
        self.region_name, self.region_start, self.region_end = region_coords
        self.chrom = self.region_name.split('_')[0]
        self.region_length = self.region_end - self.region_start

        self.reference = _get_reference(user_inputs.fa, self.chrom)

        # tiling_range is the maximum number of bases that will be
        # covered by any tiling pattern
        # this will be the number of postion needed to be memoized
        self.tiling_range = self.region_length + 2 * (self.max_tile -1)
        self.allowed_overlap = user_inputs.allowed_overlap
        self.primer_length = user_inputs.primer_length
        self.primer_length_var = user_inputs.primer_length_var

        # position memo is a list of 5-tuple recording 
        # (score, tile size, overlap, f_primer, r_primer)
        self.pos_memo = [(0, 0, 0, None, None)
                         for index in range(self.tiling_range)]
        self.primer_memo = {}
        self.aux_data = {'tiles':[], 'overlap':[]}

        logging.info('Initialised new Dp_search object')
        logging.info('Region name, start, end, length')
        logging.info('%s, %s, %s, %s'
                     %(self.region_name,
                       self.region_start,
                       self.region_end,
                       self.region_length))
        logging.info('Tilling range = %s'%self.tiling_range)
        self.dp_search()
        logging.info('dp_search() ended')
        logging.info('%s'%len(self.primer_memo))




    def dp_search(self):
        """
        This method will try all legal tiling patterns, it satisfy:
            The region is a subset of the union of the tiles
            All tiles are of the sizes specified
            All tiles has at least 1bp intersection with the region 
            Tiles can overlap each other at most allowed_overlap bp.
            (specified by user_inputs, see __init__)

        DP memoization:
        The best scored tiling pattern for each positions are recorded
        in the pos_memo and all scored primers are recored in the
        primer_memo

        """
        logging.info('dp_search() started')

        # the bottom up DP start at the first base of the region
        # since legal tile must intersect at least 1 base
        # but ends PREMATURELY at the end of the region + min_tile
        # so as to ensure that all tiles intersect with the region
        # Though we can allow tile go beyond, that case is
        # handled in a separate for loop
        # since for each position, an extra check has to be in place
        # to ensure the tile still intersect the region
        # This separation should be handled with a separate private method
        # (because i am repeating myself in codes)

        start_search = self.region_start
        end_search = self.region_end-1 + self.max_tile
        # notice 0 based indexing for the end of region is self.region_end -1
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = 0
            for t_size in self.tile_sizes:
                if pos - t_size < self.region_end: # if tile still overlap
                    # a tile is defined as (start pos, tile size)
                    tile = (pos - t_size +1, t_size)
                    # select the best primer with respect to primer length
                    # and return the sum of their scores and primers themselves.
                    tile_score, f_primer, r_primer = self.best_primers_in_tile(tile)
                    for overlap in range(self.allowed_overlap):
                        # Access the score of the previous tiling in pos_memo
                        # if this tile is chosen and take the best one,
                        # allowing for overlap of this and previous tile
                        # To access previous pos:
                        # 1) bring pos back to 0 if pos is at self.region_start
                        # 2) add self.max_tile -1 to go to current pos_memo entry
                        # 3) shift back by t_size to get previous tiles score
                        # 4) but allow for the overlaping.
                        prefix_pos = pos - self.region_start \
                                     + self.max_tile -1 \
                                      - t_size + overlap
                        prefix_score = self.pos_memo[prefix_pos][0]
                        total_score = prefix_score + tile_score

                        if total_score > best_score:
                            best_score = total_score
                            best_f = f_primer
                            best_r = r_primer
                            best_overlap = overlap
                            best_tile_size = t_size
                            logging.info('Best score, f, r, overlap, tile_size')
                            logging.info('%s %s %s %s %s'
                                         %(best_score,
                                           best_f,
                                           best_r,
                                           best_overlap,
                                           best_tile_size))
            position_in_memo = pos - self.region_start + self.max_tile -1
            #logging.info('Position in memo: %s'%position_in_memo)
            self.pos_memo[position_in_memo] = (best_score,
                                               best_tile_size,
                                               best_overlap,
                                               best_f,
                                               best_r)

    def best_primers_in_tile(self, tile):
        '''
        Take a tile in the reference region and look for the best
        forward and reverse primers (best with respect to their length)
        '''
        logging.info('Start choosing best primers in tile')
        best_f_score = 0
        best_r_score = 0
        tile_start, t_size = tile
        tile_end = tile_start + t_size -1
        var = self.primer_length_var
        primer_length = self.primer_length


        for vary in range(-var, var +1):
            this_primer_length = primer_length + vary
            # deal with forward primer
            f_primer = (tile_start -1, this_primer_length, 'f')
            if f_primer in self.primer_memo:
                f_scores = self.primer_memo[f_primer]
            else:
                f_sequence = get_primer_seq(self.reference, f_primer)
                f_scores = self.score_primer(f_sequence, direction='f')
                self.primer_memo[f_primer] = f_scores
                #logging.info('new primer: %s'%str(f_primer))
                #logging.info('Size of primer memo: %s'%len(self.primer_memo))

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length, 'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = get_primer_seq(self.reference, r_primer)
                r_scores = self.score_primer(r_sequence, direction='r')
                self.primer_memo[r_primer] = r_scores
                #logging.info('new primer: %s'%str(r_primer))
                #logging.info('Size of primer memo: %s'%len(self.primer_memo))


            if f_scores[0] > best_f_score:
                best_f_score = f_scores[0]
                best_f = f_primer
            if r_scores[0] > best_r_score:
                best_r_score = r_scores[0]
                best_r = r_primer
        total_score = best_f_score + best_r_score
        best_result = (total_score, best_f, best_r)
        logging.info('Finished choosing best primers in tile')
        logging.info('Chosen primer:\n Score:%s\nForward:%s\nReverse:%s'
                     %best_result)

        return best_result


##################

def get_primer_seq(reference, primer_specification):
    '''
    Return the sequence of primer from the reference based
    on its specification.

    primer_specification : (3' position, length, direction)
    '''
    three_prime_pos, length, direction = primer_specification
    if direction == 'f':
        return reference[three_prime_pos - length+1: three_prime_pos +1]
    elif direction == 'r':
        seq = reference[three_prime_pos: three_prime_pos + length]
        # return the reverse complement if the primer is a reverse primer
        return rev_complement(seq)


def _get_reference(fa_path, chrom):
    """
    Intended as a function used when initialising the Dp_search class.
    (see __init__) It returns the whole sequence of chromosome
    of concern in this Dp_search object
    """
    file_path = os.path.join(fa_path, chrom + '.fa')
    seq_read = SeqIO.read(file_path, 'fasta')

    # This function requires that each fasta file contains
    # a single chromosome only and the file name is consitent
    # with the chromosome name (the first column of BED-file)
    return seq_read.seq

def pick_primer_set(searcher):
    '''
    Look through starting from the end of pos_memo, pick the best scored
    "starting position" and then trace back using the tile_size and overlap
    chosen at each position to recover the optimal primer set.

    '''
    if not filter(lambda x: x != (0,0,0, None, None),searcher.pos_memo):
        logging.info('Warning the search has not been done')
        return
    # look from the end of pos_memo until it hits
    # the end of the region and determine the best starting position
    best_reverse_start = None # this will be a negative integer
    best_start_score = 0
    for pos in range(-1, -searcher.max_tile -1, -1):
        pos_info = searcher.pos_memo[pos]
        pos_score = pos_info[0]
        if pos_score > best_start_score:
            best_start_score = pos_score
            best_reverse_start = pos

    region_size_remaining = searcher.region_length \
                            + searcher.max_tile \
                            + best_reverse_start
                             # reverse_start is negative

    position = best_reverse_start
    # while we are still in the region
    # notice the definition of "region" is loosen to include
    # the overhang of the last tile
    while abs(position) < region_size_remaining:
        (pos_score, tile_size,
         overlap, f_primer, r_primer) = searcher.pos_memo[position]
        yield f_primer, r_primer

        #searcher.primer_set.append(f_primer)
        #searcher.primer_set.append(r_primer)

        # Beware of tile_size, overlap = 0 case, an infinite loop result
        if tile_size == 0 and overlap == 0:
            logging.info("There is error and this will go into an infinite loop")
            return

        position -= (tile_size - overlap)
        searcher.aux_data['tiles'].append(tile_size)
        searcher.aux_data['overlap'].append(overlap)



def write_aux(data_dic, auxfile_name):
    with open(auxfile_name,'w') as auxfile:
        header = ['chrom', 'region_start', 'region_end', 'region_length', 'time_taken',
                  'mean_tile', 'num_tile','std_deviation_tile', 'mean_overlap',
                  'num_primers_scored']
        auxfile.write('\t'.join(header) + '\n')
        for coords, data in data_dic.items():
            chrom, start, end = coords
            tiles, overlap, length, time_taken, num_primers = data
            output = map(str,
                            [chrom, start, end, length, time_taken,
                            mean(tiles), len(tiles), std_deviation(tiles),
                            mean(overlap), num_primers]
                        )
            auxfile.write('\t'.join(output) + '\n')


def write_primer(out_writer, searcher):
    '''
    This function take in a Dp_search object, call pick_primer_set()
    which yield as a generators forward and reverse primer pair
    for a tile.
    '''
    position_shift = searcher.region_start - searcher.max_tile +1
    region_start = searcher.region_start
    region_length = searcher.region_length
    tile_number = 0
    for f_primer, r_primer in pick_primer_set(searcher):
        tile_number += 1
        f_3_prime, f_length, f_dir = f_primer
        r_3_prime, r_length, r_dir = r_primer
        # using BED file, 0-based open-ended indexing on 
        # the sense strand of DNA
        f_start = f_3_prime - f_length + 1
        f_end = f_3_prime + 1
        r_start = r_3_prime
        r_end = r_3_prime + f_length

        f_seq = get_primer_seq(searcher.reference, f_primer)
        r_seq = get_primer_seq(searcher.reference, r_primer)

        f_score_data = map(str, searcher.primer_memo[f_primer])
        r_score_data = map(str, searcher.primer_memo[r_primer])

        basic_name = searcher.region_name + '_' \
                     + str(region_start) + '_' \
                     + str(region_length) + '_'

        f_name = basic_name + f_dir + str(tile_number)
        r_name = basic_name + r_dir + str(tile_number)

        f_output = [f_name] \
                   + map(str, [f_start, f_end, f_seq, f_length]) \
                   + f_score_data

        r_output = [r_name] \
                   + map(str, [r_start, r_end, r_seq, r_length]) \
                   + r_score_data

        map(out_writer.writerow, [f_output, r_output])
    logging.info("Finish writing primers of %s %s %s to file"
                  %(searcher.region_name, 
                    searcher.region_start, 
                    searcher.region_end))


def handle_bedfile(bedfile):
    file = open(bedfile)
    reader = csv.reader(file,delimiter='\t')
    bed_row = []
    bed_dictionary = {}
    for line in reader:
        chromo, start, end, name = line 
        start, end = map(int, [start, end])
        # we dont want underscore character in the bed specified naming
        # to interfere with our own naming convention
        bed_specified_name = name.split(',')[0].replace('_','')
        region_name = chromo + '_' + bed_specified_name 
        # last term is potentially [''], ie a list of single empty string
        bed_row.append((region_name, start, end))
        if chromo in bed_dictionary:
            bed_dictionary[chromo].append((start, end, name))
        else:
            bed_dictionary[chromo] = [(start, end, name)]
    return bed_row #bed_dictionary

if __name__ == '__main__':
    main()
