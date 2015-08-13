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
    - outfile
    - Auxilarry output file for tile_size distributions and other data

"""
from argparse import ArgumentParser
import logging
import sys
import score
from Bio import SeqIO
import time
from utils import handle_bedfile, rev_complement, write_aux
###############################################################################

DEFAULT_LOG_FILE ='primer_design.log'
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
    parser.add_argument('--outfile', metavar = "OUTPUT_FILE", type = str,
                        default = 'primer_out.txt',
                        help = ''' A string specifying the name of the output
                        file''')
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
    parser.add_argument('--primer_length_var', metavar = 'PRIMER_LENGTH_VARIATION',
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
    parser.add_argument('--score_func', metavar = "SCORING_FUNCTION",
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
    parser.add_argument('--target_tm', metavar="TARGET_TM", type = float,
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
    print('These are your inputs:')
    print(user_inputs)
    
    start_log(user_inputs.log)
    
    scoring = score.Score(user_inputs)
    score_func = scoring.score_func

    # file handling
    bed_coords = handle_bedfile(user_inputs.bed)
    outfile_header = ['name', 'start', 'end', 'sequence', 'length',
                      'tm', 'entropy', 'hairpin', 'gc', 'gc_clamp', 'run']
    aux_outfile_header = ['tiles', 'overlap']
    outfile = open(user_inputs.outfile, 'w')
    outfile.write('\t'.join(outfile_header) + '\n')

    # collect all aux data from each exon search
    # a post process statistics is then done on it.
    auxiliary_data = {} # lists for tiles and overlap 


    for chromo in bed_coords:
        for region_coord in bed_coords[chromo]:
            before_time = time.time()

            coords = (chromo, region_coord[0], region_coord[1])
            searcher = DP_exon_search(user_inputs, coords, score_func)
            searcher.dp_search()
            #print(searcher.pos_memo)
            searcher.pick_primer_set()
            searcher.write_output(outfile)

            after_time = time.time()
            time_taken = after_time - before_time
            print("%s %s %s"%coords)
            print("length %i"%searcher.exon_length)
            print("Search-pick-write time: %s\n"%time_taken)

            # assemble aux_data
            auxiliary_data[coords] = [list(searcher.aux_data['tiles']),
                                    list(searcher.aux_data['overlap']),
                                    searcher.exon_length, time_taken,
                                    len(searcher.primer_memo)]
            tiles_chosen = searcher.aux_data['tiles']
            overlaps = searcher.aux_data['overlap']

    write_aux(auxiliary_data, user_inputs.auxfile)
    outfile.close()


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
        self.tile_sizes = [self.min_tile + extend 
                      for extend in range(self.max_tile - self.min_tile + 1)]
        self.chrom, self.exon_start, self.exon_end = exon_coords
        self.exon_length = self.exon_end - self.exon_start

        self.background = self._get_background()

        self.tiling_range = self.exon_length + 2 * (self.max_tile -1)
        self.allowed_overlap = user_inputs.allowed_overlap
        self.primer_length = user_inputs.primer_length
        self.primer_length_var = user_inputs.primer_length_var
        # initiallise both memos
        self.pos_memo = [(0,0,0,None, None) 
                         for index in range(self.tiling_range)]
        self.primer_memo = {}

        # chosen primer_set, 
        # non-empty only if pick_primer_set was successfully called
        self.primer_set = []

        # initialise an auxiliary data record
        self.aux_data = {'tiles':[],'overlap':[]}



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
        # Though we can allow tile go beyond, that case is 
        # handled in a separate for loop.
        # This separation should be handled with a separate private method
        # (because i am repeating myself in codes)

        start_search = self.exon_start
        end_search = self.exon_end-1 + self.min_tile
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
                    # Access the score of the previous tiling in pos_memo
                    # if this tile is chosen and take the best one, 
                    # allowing for overlap of this and previous tile
                    # To access previous pos:
                    # 1) bring pos back to 0 if pos is at self.exon_start
                    # 2) add self.max_tile -1 to go to current pos_memo entry
                    # 3) shift back by t_size to get previous tiles score
                    # 4) but allow for the overlaping.
                    prefix_pos = pos - self.exon_start \
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
            position_in_memo = pos - self.exon_start + self.max_tile -1
            self.pos_memo[position_in_memo] = ( best_score, 
                                                best_tile_size, 
                                                best_overlap, 
                                                best_f, 
                                                best_r)


        # redefine start and end searching position to handle
        # the tiles that go beyong exon + min_tile -1
        start_search = end_search
        end_search = self.exon_end-1 + self.max_tile 
        # exon_end is bed file coord,
        # hence the actual exon_ending position is exon_end -1
        for pos in range(start_search, end_search):
            best_score = 0
            best_f = None
            best_r = None
            best_overlap = 0
            best_tile_size = None
            for t_size in self.tile_sizes:
                if pos - t_size < self.exon_end: # if tile still overlap
                    tile = (pos - t_size +1, t_size)
                    tile_primer_pair = self.best_primers_in_tile(tile)
                    tile_score, f_primer, r_primer = tile_primer_pair
                    for overlap in range(self.allowed_overlap):
                        prefix_pos = pos - self.exon_start \
                                     + self.max_tile -1 \
                                      - t_size + overlap
                        prefix_score = self.pos_memo[prefix_pos][0]
                        total_score = prefix_score + tile_score

                        if total_score > best_score:
                            best_score  = total_score
                            best_f = f_primer
                            best_r = r_primer
                            best_overlap = overlap
                            best_tile_size = t_size
            position_in_memo = pos - self.exon_start + self.max_tile -1
            self.pos_memo[position_in_memo] = (best_score, 
                                              best_tile_size, 
                                              best_overlap, 
                                              best_f, 
                                              best_r)
        print("DP search ended")

    def best_primers_in_tile(self, tile):
        best_f = None
        best_r = None
        best_f_score = 0
        best_r_score = 0
        tile_start, t_size = tile
        tile_end = tile_start + t_size -1
        var = self.primer_length_var
        primer_length = self.primer_length


        for vary in range(-var, var +1):
            this_primer_length = primer_length + vary
            # deal with forward primer
            f_primer = (tile_start -1, this_primer_length,'f')
            if f_primer in self.primer_memo:
                f_scores = self.primer_memo[f_primer]
            else:
                f_sequence = self.get_seq(f_primer)
                f_scores = self.score_primer(f_sequence, direction = 'f')
                self.primer_memo[f_primer] = f_scores

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length,'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = self.get_seq(r_primer)
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


    def _get_background(self):
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
        seq_read = SeqIO.read('./fasta/%s.fa'%self.chrom, 'fasta')
        # This function is specific to the way we put our files
        # This should be generalised to handle other situations
        # for instance, the fasta files are with multiple chrom
        return seq_read.seq

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
        for pos in range(-1, -self.max_tile -1,-1):
            pos_info = self.pos_memo[pos]
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
            position_shift = self.exon_start - self.max_tile +1
            start = three_prime_pos - length +1 + position_shift
            end = three_prime_pos +1 + position_shift
            primer_seq = self.get_seq(primer)
            score_data = map(str,self.primer_memo[primer])
            primer_name = self.chrom + '_' \
                          + str(start) + '_' \
                          + str(end) + '_' \
                          + direction + str(self.primer_set.index(primer) +1) 

            output = map(str,[primer_name, start, end,
                      primer_seq, length]) + score_data
            row_output = '\t'.join(output)+'\n'
            outfile.write(row_output)
        outfile.flush()
        print("Finish writing primers of %s %s %s to file"
                %(self.chrom, self.exon_start, self.exon_end))







#class Primer(object):
#    def __init__(self, chrom, start, end, direction):
#        assert direction == 'f' or direction == 'r'
#         self.chrom = chrom
#         self.direction = direction
#         self.start_pos = start
#         self.end_pos = end
#         self.length = end - start
#         self.raw_scores = None
#         self.scores = None
#         self.total_score = None
        # self.sequence = None
        # a decision to make here
        # to avoid storing the sequence we can 
        # Only calculate this when get_sequence is called
        # when the sequence is needed or,
        # to avoid recalculaton, the sequence is stored
        # as an attribute.

#    def get_sequence(self):
#        # return sequence
#        # or if self.sequence == None:
#        #    calculate
#        #    self.sequence == sequence
#        #    else:
#        #       return self.sequence?? hmm...
#        pass
#
#    def get_3prime(self):
#        pass







if __name__ == '__main__':
    main()
