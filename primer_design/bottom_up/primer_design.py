"""
Commandline arguments (user input):
    - BED file
    - fasta choice?
    - Heel sequence
    - tiles (from min to max inclusive)
    - optimal primer length
    - primer length variation
    - tile overlap allowed
    - scoring function used
    - target Tm
    - tm underachieve penalty weight
    - tm funciton used
    - Concentration of Na, K, Mg, Tris, dNTPs, DMSO, dnac1, dnc2, 
         * DMSO not added yet
    - saltcorr 
    - GCweight
    - Auxilarry output file for tile_size distributions and other data

"""
from argparse import ArgumentParser
import logging
import score
from Bio import SeqIO
import time
from utils import handle_bedfile, rev_complement
###############################################################################

DEFAULT_LOG_FILE =' primer_design.log'
def parse_args():
    'Parse command line arguments for the program'
    parser = ArgumentParser(description = 'Hiplex primer design tool')
    parser.add_argument('--log', metavar = "LOG_FILE", type = str,
                         default = DEFAULT_LOG_FILE, 
                         help = 'Log file. Default to %s'%DEFAULT_LOG_FILE)
    parser.add_argument('--bed', metavar = "BED_FILE", type = str,
                        required=True, 
                        help = ''' A BED file specifying all the coordinates
                        of the regions of interest ''')
    parser.add_argument('--sense_heel', metavar = 'SENSE_HEEL', type = str,
                        required=True,
                        help = '''A string of nucleic acid specifying the 
                        sense heel strand''')
    parser.add_argument('--antisense_heel', metavar = 'ANTISENSE_HEEL', type=str,
                        required=True,
                        help = '''A string of nucleic acid specifying the 
                        antisense heel strand''')
    parser.add_argument('--tiles', metavar="MAX_MIN_TILE_SIZE", type = int,
                        required=True, nargs = 2, 
                        help = ''' A pair of integers specifying the maximum
                        and the minimum tile sizes inclusive''')
    parser.add_argument('--primer_length', metavar = 'PRIMER_LENGTH', type = int,
                        default = 20, 
                        help = '''An integer specifying the optimal primer length
                        Defaulted to 20''')
    parser.add_argument('--primer_length_var', metavar = 'PRIMER_LENGTH_VARIATION",
                        type = int, default = 8,
                        help = '''An integer specifying the amount of variation 
                        from the optimal primer length. 
                        Defaulted to 8
                        Eg: optimal length of 20 and variation of 5 gives
                            primer length ranging from 15 to 25 inclusive''')
    parser.add_argument('--allowed_overlap', metavar ="ALLOWED_TILE_OVERLAP",
                        type=int, default = 5,
                        help = '''An integer specifying the allowed overlaping 
                        between successive tiles.
                        Defaulted to 5''')
    parser.add_argument('--score_func', metarvar = "SCORING_FUNCTION",
                        type = str, default = 'score_Lp', 
                        help = ''' A string specifying the name of the
                        scoring function to be used.
                        Options include:
                           score_Lp
                           score_linear
                           ''')
    parser.add_argument('--tm_func', metavar = 'TM_FUNCTION', type = str,
                       default = 'Tm_NN',
                       help = '''A string specifying the name of the 
                       melting temperature prediction algorithm used.
                       Defaulted to Tm_NN
                       Options include:
                           Tm_NN  (Nearest Neighbor thermodynamics)
                           Tm_GC  (Prediction using GC content of sequence)
                           Tm_Wallace (Wallace's 2-4 rule)
                           ''')
    parser.add_argument('--target_tm', metarvar="TARGET_TM", type = float,
                        default = 64.0, 
                        help = '''A floating point number specifying the 
                        target melting point of the reaction mixture
                        in degree Celsius.
                        Defaulted to 64.0 degC''')
    parser.add_argument('--tm_underachieve', metavar="TM_UNDERACHIEVE_PENALTY",
                        type = float, default = 1.0,
                        help = '''A floating point number specifying the 
                        penalty weight given to an underachieving primer
                        in term of tm, ie has lower melting temperature than
                        target_tm. Lower melting temperature should be score
                        harsher than higher ones, there for the value should
                        be > 1.''')
    parser.add_argument('--saltcorr', metavar ="SALTCORR", type = int,
                        default = 5, 
                        help = ''' An integer from 0 to 7 inclusive indicating
                        the saltcorrection method to be used during 
                        melting temperature prediction.
                        0 - no salt corrections. 
                        To see individual method, refer to
                        Biopython Bio.SeqUtils.MeltingTemp module
                        Default to 5.''')
    parser.add_argument('--conc', metavar = "CONCENTRATIONS", type = float,
                        nargs=7, default = [50, 0, 0, 0, 0, 25, 25],
                        help = ''' 7 floating point numbers which give the
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
    parser.add_argument('--gc_weight', metavar = "GC_WEIGHT", type = float,
                        default = 1.5, 
                        help = '''A floating point number >1 specifying the
                        weight given to G-C binding over A-T binding during
                        the scoring of hairpins.
                        Default to 1.5''')
    parser.add_argument('--auxfile', metavar = "AUXILIARY_FILE", type = str,
                        default = 'auxiliary_primer_out.txt',
                        help = '''A string specifying the name of an auxiliary
                        output files that will record some additional data
                        and statistics''')

    return parser.parse_args()

def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOG_FILE'''
    logging.basicConfig(filename = log,
                        level = logging.DEBUG,
                        filemode = 'w',
                        format = '%(asctime)s %(message)s',
                        datefmt = '%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


def main():
    user_inputs = parse_args()
    start_log(user_inputs.log)
    scoring = Score(user_inputs)
    score_func = scoring.score_func
    # handle scoring here
    bed_coords = handle_bedfile(args.bed)
    min_tile, max_tile = ags.tiles
    outfile = open(args.outfile, 'w')
    auxiliary_data = {} 
    # collect all aux data from each exon search
    # a post process statistics is then done on it.
    for chromo in bed_coords:
        for region_coord in bed_coords[chromo]:
            coords = (chromo, region_coord[0], region_coord[1])
            searcher = DP_exon_search(user_inputs, coords, score_func)
            searcher.dp_search()
            searcher.pick_primer_set()
            searcher.write_output()



class DP_exon_search(object):
    """
    DP_exon_search object initialise the 'search for best primer set'
    process. It takes in all command-line arguments 
       *may or may not be usefull
    , a scoring function (* see Score() class) and 
    exon_coords :: (chrom, start, end), eg (chr1, 0, 100).

    It initialise a class level pos_memo and primer_memo for 
    Dynamic programming purpose.
    """
    def __init__(self, user_inputs, exon_coords, score_func):
        self.score_primer = score_func

        self.min_tile, self.max_tile = user_inputs.tiles
        self.tiles = [min_tile + extend for extend in range(max_tile - min_tile + 1)]
        self.chrom, self.exon_start, self.exon_end = exon_coords
        self.exon_length = self.exon_end - self.exon_start
        self.background = self._get_background(exon_coords, max_tile)
        self.background_length = len(self.background)
        self.allowed_overlap = user_inputs.allowed_overlap
        self.primer_length = user_inputs.primer_length
        self.primer_length_var = user_inputs.primer_length_var
        # initiallise both memos
        assert self.background_length = self.exon_length + 2 * (max_tile -1)
        self.pos_memo = [ (0,0,None, None) for index in range(self.background_length)]
        self.primer_memo = {}

        # chosen primer_set, non-empty only if pick_primer_set was successfully called
        self.primer_set = []

        # initialise an auxiliary data record
        self.aux_data = {}



    def dp_search(self):
        """
        This method will try all legal tiling patterns, that is:
            The exon is a subset of the union of the tiles
            All tiles are of the sizes specified
            All tiles has at least 1bp intersection with the exon
            Tiles can overlap each other at most allowed_overlap bp.
            (specified by user_inputs, see __init__)

        DP memoization:
        The best scored tiling pattern for each positions are recorded 
        in the pos_memo and all scored primers are recored in the
        primer_memo
        """

        # the bottom up DP start at the first base of the exon
        # but ends PREMATURELY at the end of the exon + min_tile -1
        # so as to ensure that all tiles intersect with the region.
        # Though we can allow tile go beyond that case is 
        # handled in a separate for loop.
        # This separation should be handled with a separate private method
        # (because i am repeating myself

        exon_end_pos = self.background_length - max_tile
        start_search = self.max_tile - 1
        end_search = exon_end_pos + min_tile
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = None
            for t_size in self.tile_sizes:
                # a tile is defined as (start pos, tile size)
                tile = (pos - t_size +1, t_size)
                tile_primer_pair = self.best_primers_in_tile(tile)
                tile_score, f_primer, r_primer = tile_primer_pair
                for overlap in range(self.allowed_overlap):
                    prefix_pos = pos - t_size + overlap
                    prefix_score = self.pos_memo[prefix_pos][0]
                    total_score = prefic_score + suffix_score

                    if total_score > best_score:
                        best_score = total_score
                        best_f = f_primer
                        best_r = r_primer
                        best_overlap = overlap
                        best_tile_size = t_size
            self.pos_memo[pos] = (best_score, best_tile_size, best_overlap, best_f, best_r)

        # redefine start and end searching position to handle
        # the tiles that go beyong exon + min_tile -1
        start_search = end_search
        end_search = self.background_length
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = None
            for t_size in self.tile_sizes:
                if pos - t_size < exon_end_pos: # if tile still overlap
                    tile = (pos - t_size +1, t_size)
                    tile_primer_pair = self.best_primers_in_tile(tile)
                    tile_score, f_primer, r_primer = tile_primer_pair
                    for overlap in range(self.allowed_overlap):
                        prefix_pos = pos - t_size + overlap
                        prefix_score = self.pos_memo[prefix_pos][0]
                        total_score = prefix_score + suffix_score

                        if total_score > best_score:
                            best_score  = total_score
                            best_f = f_primer
                            best_r = r_primer
                            best_overlap = overlap
                            best_tile_size = t_size
            self.pos_memo[pos] = (best_score, best_tile_size, best_overlap, best_f, best_r)

    def best_primers_in_tile(self, tile):
        best_f = None
        best_r = None
        best_f_score = 0
        best_r_score = 0
        tile_start, t_size = tile
        tile_end = tile_start + tile_size -1
        var = self.primer_length_var
        primer_length = self.primer_length

        for vary in range(-var, var +1):
            this_primer_length = primer_length + vary
            # deal with forward primer
            f_primer = (tile_start -1, this_primer_length,'f')
            if f_primer in self.primer_memo:
                f_scores = self.primer_memo[f_primer]
            else:
                f_sequence = background[tile_start - this_primer_length: tile_start]
                f_scores = self.score_primer(f_sequence, direction = 'f')
                self.primer_memo[f_primer] = f_scores

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length,'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = background[tile_end +1 : tile_end + this_primer_length +1]
                r_scores = self.score_primer(r_sequence, direction = 'r')
                self.primer_memo[r_primer] = r_scores


            if f_scores[0] > best_f_score:
                best_f_score = f_scores[0]
                best_f = f_primer
            if r_scores[0] > best_r_score:
                best_r_score = r_scores[0]
                best_r = r_primer
        total_score = best_f_score + best_r_score

        return (total_score, f_primer, r_primer)


    def _get_background(self, exon_coords, max_tile):
        """
        This is a private method used when initialising the class.
        (see __init__) It takes in exon_coords of concern and
        the maximum tile size available and return the relevant region
        (called background)in the reference genome as a Bio.Seq object
            background = overhang + exon + overhang
                where overhang is the flanking sequence of size = max_tile-1
        The overhang allow the tiles to go beyond the exon but has
        at least 1 bp intersection. 
        """
        # there is a repeat in parsing user_inputs here, see __init__
        chrom, start, end = exon_coords
        region_start = start - max_tile +1
        region _end = end + max_tile -1
        seq_read = SeqIO.read('/fasta/%s.fa'%chrom, 'fasta')
        # This function is specific to the way we put our files
        # This should be generalised to handle other situations
        # for instance, the fasta files are with multiple chrom
        seq = seq_read.seq[region_start:region_end]
        return seq

    def get_seq(self, primer_specification):
        three_prime_pos, length, dir = primer_specification
        if dir == 'f':
            return self.background[three_prime_pos - length+1: three_prime_pos +1]
        elif dir == 'r':
            seq = self.background[three_prime_pos: three_prime_pos + length]
            return rev_complement(seq)
        else:
            print("Warning: please specify your primer direction 'f' or 'r'")


    def pick_primer_set(self):
        if not self.primer_memo:
            print('Warning: The search has not been run')
            return
        
        
        best_reverse_start = None
        best_start_score = 0
        for pos in range(-1, -self.max_tile,-1):
            pos_info = pos_memo[pos]
            pos_score = pos_info[0]
            if pos_score > best_start_score:
                best_start_score = pos_score
                best_reverse_start = pos
        region_size_remaining = self.exon_length \
                                + self.max_tile \
                                + best_reverse_start
                                 # reverse_start is negative

        position = best_reverse_start
        while abs(position) < region_size_remaining:
            pos_score, 
            tile_size, 
            overlap, 
            f_primer, r_primer = self.pos_memo[position]

            self.primer_set.append(f_primer)
            self.primer_set.append(r_primer)

            position -= (tile_size - overlap)
        return primer_set


    def write_output(self, outfile):
        for primer in primer_set:
            three_prime_pos, length, direction = primer
            five_primer_pos = three_primer_pos + length -1
            primer_seq = self.get_seq(primer)
            score_data = list(self.primer_memo[primer])
            primer_name = self.chrom + '_' \
                          + str(three_prime_pos) + '_' \
                          + str(five_primer_pos) + '_' \

            output = [primer_name, primer_seq] + score_data
            row_output = '\t'.join(output)+'\n'
            outfile.write(row_output)
        outfile.flush()







class Primer(object):
    def __init__(self, chrom, start, end, direction):
        assert direction == 'f' or direction == 'r'
        self.chrom = chrom
        self.direction = direction
        self.start_pos = start
        self.end_pos = end
        self.length = end - start
        self.raw_scores = None
        self.scores = None
        self.total_score = None
        # self.sequence = None
        # a decision to make here
        # to avoid storing the sequence we can 
        # Only calculate this when get_sequence is called
        # when the sequence is needed or,
        # to avoid recalculaton, the sequence is stored
        # as an attribute.

    def get_sequence(self):
        # return sequence
        # or if self.sequence == None:
        #    calculate
        #    self.sequence == sequence
        #    else:
        #       return self.sequence?? hmm...
        pass

    def get_3prime(self):
        pass







if __name__ == '__main__':
    main()
