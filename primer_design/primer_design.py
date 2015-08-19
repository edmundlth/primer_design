"""
                Primer design
Authors: Bernard J. Pope, 
         Daniel J. Park,
         Tu Nguyen,
         Edmund Lau
Email: elau1@student.unimelb.edu.au

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
                                 run time, tiles size mean, variance,
                                 overlap mean, variance

    This program uses bottom up dynamic programming (see description below)
    to sift through the search space of all possible primer set
    which will completely tile all regions specified. We allow tiles
    to overlap and primers to vary in length up to the corresponding
    user specified constrains.


Conventions:
    - Indexing are 0-based (similar to BED format and python string)
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
"""
from argparse import ArgumentParser
import logging
import sys
import time

from Bio import SeqIO

from utils import handle_bedfile, rev_complement, write_aux
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
DEFAULT_SALTCORR = 5
DEFAULT_CONC = [50, 0, 0, 0, 0, 25, 25]
DEFAULT_GC_WEIGHT = 1.5
DEFAULT_AUXFILE = 'auxiliary_primer_out.tsv'


def parse_args():
    'Parse command line arguments for the program'
    parser = ArgumentParser(description='Hiplex primer design tool')
    parser.add_argument('--log', metavar="LOG_FILE", type=str,
                        default=DEFAULT_LOG_FILE,
                        help='Log file. Default to %s'%DEFAULT_LOG_FILE)
    parser.add_argument('--bed', metavar="BED_FILE", type=str,
                        required=True,
                        help=''' A BED file specifying all the coordinates
                        of the regions of interest ''')
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
                        Defaulted to 20''')
    parser.add_argument('--primer_length_var', metavar='PRIMER_LENGTH_VARIATION',
                        type=int, default=DEFAULT_PRIMER_LENGTH_VAR,
                        help='''An integer specifying the amount of variation
                        from the optimal primer length.
                        Defaulted to 8
                        Eg: optimal length of 20 and variation of 5 gives
                            primer length ranging from 15 to 25 inclusive''')
    parser.add_argument('--allowed_overlap', metavar="ALLOWED_TILE_OVERLAP",
                        type=int, default=DEFAULT_ALLOWED_OVERLAP,
                        help='''An integer specifying the allowed overlaping
                        between successive tiles.
                        Defaulted to 5''')
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
                        Defaulted to Tm_NN''')
    parser.add_argument('--target_tm', metavar="TARGET_TM", type=float,
                        default=DEFAULT_TARGET_TM,
                        help='''A floating point number specifying the
                        target melting point of the reaction mixture
                        in degree Celsius.
                        Defaulted to 64.0 degC''')
    parser.add_argument('--tm_underachieve', metavar="TM_UNDERACHIEVE_PENALTY",
                        type=float, default=DEFAULT_TM_UNDERACHIEVE,
                        help='''A floating point number specifying the
                        penalty weight given to an underachieving primer
                        in term of tm, ie has lower melting temperature than
                        target_tm. Lower melting temperature should be score
                        harsher than higher ones, there for the value should
                        be > 1.''')
    parser.add_argument('--saltcorr', metavar="SALTCORR", type=int,
                        default=DEFAULT_SALTCORR,
                        choices=range(0,8),
                        help=''' An integer from 0 to 7 inclusive indicating
                        the saltcorrection method to be used during
                        melting temperature prediction.
                        0 - no salt corrections.
                        To see individual method, refer to
                        Biopython Bio.SeqUtils.MeltingTemp module
                        Default to 5.''')
    parser.add_argument('--conc', metavar="CONCENTRATIONS", type=float,
                        nargs=7, default=DEFAULT_CONC,
                        help=''' 7 floating point numbers which give the
                        the concentration of the following chemical species:
                            - Na (in mM)
                            - K  (in mM)
                            - Mg (in mM)
                            - Tris-HCl (in mM)
                            - dNTPs (in mM)
                            - dnac1 (in nM) note the nM. This is the
                                            concentration of the more abundant
                                            DNA species (likely to be primer)
                            - dnac2 (in nM) Concentration of the less abundant
                                            DNA species.

                            Default to [50, 0, 0, 0, 0, 25, 25]
                            ''')
    parser.add_argument('--gc_weight', metavar="GC_WEIGHT", type=float,
                        default=DEFAULT_GC_WEIGHT,
                        help='''A floating point number >1 specifying the
                        weight given to G-C binding over A-T binding during
                        the scoring of hairpins.
                        Default to 1.5''')
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
    ''' Main() function'''
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

    # parse bed file and return a dictionary of {chrom: [(start, end), (start, end)]}
    bed_coords = handle_bedfile(user_inputs.bed)

    # write output file header
    outfile_header = ['name', 'start', 'end', 'sequence',
                      'length', 'tm', 'entropy', 'hairpin',
                      'gc', 'gc_clamp', 'run']
    with open(user_inputs.outfile, 'w') as outfile:
        outfile.write('\t'.join(outfile_header) + '\n')

        # initialise a auxiliary data dictionary that will collect data
        # from each run which will be written into an auxiliary file
        # a post process statistics is then done on it.
        auxiliary_data = {}

        # Now we loop through all specified regions,
        # search for optimal primer set (prepare the memos),
        # pick the optimal primer set,
        # write primers and their datas to output files
        for chromo in bed_coords:
            for region_coord in bed_coords[chromo]:
                before_time = time.time() # time each run

                coords = (chromo, region_coord[0], region_coord[1])
                searcher = Dp_search(user_inputs, coords, score_func)
                searcher.dp_search()
                #print(searcher.pos_memo)
                searcher.pick_primer_set()
                searcher.write_output(outfile)

                after_time = time.time()
                time_taken = after_time - before_time
                print("%s %s %s"%coords)
                print("length %i"%searcher.region_length)
                print("Search-pick-write time: %s\n"%time_taken)

                # assemble aux_data
                auxiliary_data[coords] = [list(searcher.aux_data['tiles']),
                                          list(searcher.aux_data['overlap']),
                                          searcher.region_length, time_taken,
                                          len(searcher.primer_memo)]

        # The primer design process is now completed
        # write auxiliary data into auxfile
        write_aux(auxiliary_data, user_inputs.auxfile)


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
        self.chrom, self.region_start, self.region_end = region_coords
        self.region_length = self.region_end - self.region_start

        self.reference = self._get_reference()

        # tiling_range is the maximum number of bases that will be
        # covered by any tiling pattern
        # this will be the number of postion needed to be memoized
        self.tiling_range = self.region_length + 2 * (self.max_tile -1)
        self.allowed_overlap = user_inputs.allowed_overlap
        self.primer_length = user_inputs.primer_length
        self.primer_length_var = user_inputs.primer_length_var

        # initiallise both memos
        # position memo is a list of 5-tuple recording 
        # (score, tile size, overlap, f_primer, r_primer)
        self.pos_memo = [(0, 0, 0, None, None)
                         for index in range(self.tiling_range)]
        self.primer_memo = {}

        # chosen primer_set,
        # non-empty only if pick_primer_set was successfully called
        self.primer_set = []

        # initialise an auxiliary data record
        self.aux_data = {'tiles':[], 'overlap':[]}



    def dp_search(self):
        """
        This method will try all legal tiling patterns, that is:
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
        end_search = self.region_end-1 + self.min_tile
        # notice 0 based indexing for the end of region is self.region_end -1
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = None
            for t_size in self.tile_sizes:
                # a tile is defined as (start pos, tile size)
                tile = (pos - t_size +1, t_size)
                # select the best primer with respect to primer length
                # and return the sum of their scores and primers themselves.
                tile_primer_pair = self.best_primers_in_tile(tile)
                tile_score, f_primer, r_primer = tile_primer_pair
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
            position_in_memo = pos - self.region_start + self.max_tile -1
            self.pos_memo[position_in_memo] = (best_score,
                                               best_tile_size,
                                               best_overlap,
                                               best_f,
                                               best_r)


        # redefine start and end searching position to handle
        # the tiles that go beyong region + min_tile -1
        start_search = end_search
        end_search = self.region_end-1 + self.max_tile
        # region_end is bed file coord, ie 0-based
        # hence the actual region ending position is region_end -1
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = None
            for t_size in self.tile_sizes:
                if pos - t_size < self.region_end: # if tile still overlap
                    tile = (pos - t_size +1, t_size)
                    tile_primer_pair = self.best_primers_in_tile(tile)
                    tile_score, f_primer, r_primer = tile_primer_pair
                    for overlap in range(self.allowed_overlap):
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
            position_in_memo = pos - self.region_start + self.max_tile -1
            self.pos_memo[position_in_memo] = (best_score,
                                               best_tile_size,
                                               best_overlap,
                                               best_f,
                                               best_r)
        print("DP search ended")

    def best_primers_in_tile(self, tile):
        '''
        Take a tile in the reference region and look for the best
        forward and reverse primers (best with respect to their length)
        '''
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
                f_sequence = self.get_seq(f_primer)
                f_scores = self.score_primer(f_sequence, direction='f')
                self.primer_memo[f_primer] = f_scores

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length, 'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = self.get_seq(r_primer)
                r_scores = self.score_primer(r_sequence, direction='r')
                self.primer_memo[r_primer] = r_scores


            if f_scores[0] > best_f_score:
                best_f_score = f_scores[0]
                best_f = f_primer
            if r_scores[0] > best_r_score:
                best_r_score = r_scores[0]
                best_r = r_primer
        total_score = best_f_score + best_r_score

        return (total_score, best_f, best_r)


    def _get_reference(self):
        """
        This is a private method used when initialising the class.
        (see __init__) It returns the whole sequence of chromosome
        of concern in this Dp_search object
        """
        seq_read = SeqIO.read('./fasta/%s.fa'%self.chrom, 'fasta')
        # This function is specific to the way we put our files
        # This should be generalised to handle other situations
        # for instance, the fasta files are with multiple chrom
        return seq_read.seq

    def get_seq(self, primer_specification):
        '''
        Return the sequence of primer from the reference based
        on its specification.

        primer_specification : (3' position, length, direction)
        '''
        three_prime_pos, length, direction = primer_specification
        if direction == 'f':
            return self.reference[three_prime_pos - length+1: three_prime_pos +1]
        elif direction == 'r':
            seq = self.reference[three_prime_pos: three_prime_pos + length]
             # return the reverse complement if the primer is a reverse primer
            return rev_complement(seq)
        else:
            print("Warning: please specify your primer direction 'f' or 'r'")


    def pick_primer_set(self):
        '''
        Look through starting from the end of pos_memo, pick the best scored
        "starting position" and then trace back using the tile_size and overlap
        chosen at each position to recover the optimal primer set.

        '''
        if not self.primer_memo:
            print('Warning: The search has not been run')
            return

        # look from the end of pos_memo until it hits
        # the end of the region and determine the best starting position
        best_reverse_start = None # this will be a negative integer
        best_start_score = 0
        for pos in range(-1, -self.max_tile -1, -1):
            pos_info = self.pos_memo[pos]
            pos_score = pos_info[0]
            if pos_score > best_start_score:
                best_start_score = pos_score
                best_reverse_start = pos

        region_size_remaining = self.region_length \
                                + self.max_tile \
                                + best_reverse_start
                                 # reverse_start is negative

        position = best_reverse_start
        # while we are still in the region
        # notice the definition of "region" is loosen to include
        # the overhang of the last tile
        while abs(position) < region_size_remaining:
            (pos_score, tile_size,
             overlap, f_primer, r_primer) = self.pos_memo[position]

            self.primer_set.append(f_primer)
            self.primer_set.append(r_primer)

            # Beware of tile_size, overlap = 0 case, an infinite loop result
            if tile_size == 0 and overlap == 0:
                print("There is error and this will go into an infinite loop")
                return

            position -= (tile_size - overlap)
            self.aux_data['tiles'].append(tile_size)
            self.aux_data['overlap'].append(overlap)
        print("Finished picking primers based on pos memo and primer memo")


    def write_output(self, outfile):
        for primer in self.primer_set:
            three_prime_pos, length, direction = primer

            # produce bed file format coordinates for primer name
            position_shift = self.region_start - self.max_tile +1
            start = three_prime_pos - length +1 + position_shift
            end = three_prime_pos +1 + position_shift
            primer_seq = self.get_seq(primer)
            score_data = map(str, self.primer_memo[primer])
            primer_name = self.chrom + '_' \
                          + str(start) + '_' \
                          + str(end) + '_' \
                          + direction + str(self.primer_set.index(primer) +1)

            output = map(str, [primer_name, start, end,
                               primer_seq, length]) + score_data
            row_output = '\t'.join(output)+'\n'
            outfile.write(row_output)
        outfile.flush()
        print("Finish writing primers of %s %s %s to file"
              %(self.chrom, self.region_start, self.region_end))


if __name__ == '__main__':
    main()
