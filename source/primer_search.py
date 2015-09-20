import logging
import sys
import os
import time
import csv

from Bio import SeqIO

from utils import rev_complement, mean, std_deviation


class Dp_search(object):
    """
    Dp_search creates an searcher given a region of interest,
    a reference sequence where the region is embedded in,
    a scoring function that maps (purely) strings of 'ATGC'
    to scores, and user command-line input to perform 
    a dyanmic programming search and optimisation process.
    The process terminates and create 2 memo where the optimum
    solution can be found by feeding the object to pick_primer_set()
    as argument.
    """
    def __init__(self, user_inputs, region_and_ref_seq, score_func):
        """
        user_inputs :: user commandline inputs as Namespace() object
        region_and_ref_seq :: ((chrom_name, start, end), ref_seq)
        score_func :: function object
        """
        self.score_primer = score_func
        self.min_tile, self.max_tile = user_inputs.tiles
        self.tile_sizes = [self.min_tile + extend
                           for extend in
                           range(self.max_tile - self.min_tile + 1)]
        self.region_coords, self.reference = region_and_ref_seq
        (self.region_name, 
         self.region_start, 
         self.region_end) = self.region_coords
        self.chrom = self.region_name.split('_')[0]
        self.region_length = self.region_coords[2] - self.region_coords[1]
        # tiling_range is the maximum number of bases that will be
        # covered by any tiling pattern
        # this will be the number of postion needed to be memoized
        self.tiling_range = self.region_length + 2 * (self.max_tile -1)
        self.primer_length = user_inputs.primer_length
        self.primer_length_var = user_inputs.primer_length_var
        self.allowed_overlap = user_inputs.allowed_overlap
        # position memo is a list of 5-tuple recording 
        # (score, tile size, overlap, tile_count, f_primer, r_primer)
        self.pos_memo = [(0, 0, 0, 0, None, None)
                         for index in range(self.tiling_range)]
        self.primer_memo = {}
        self.aux_data = {'tiles':[], 'overlap':[]}

        logging.info('Initialised Dp_search() for %s with length %i'
                     %(str(self.region_coords), self.region_length))
        logging.info('Tilling range = %s'%self.tiling_range)
        before = time.time()
        self.dp_search()
        self.time_taken = time.time() - before
        logging.info('dp_search() ended with %i entries in primer_memo'
                     %(len(self.primer_memo)))


    def dp_search(self):
        """
        dp_search:
        the bottom up DP start at the first base of the region
        since legal tile must intersect at least 1 base
        and it ends when even the maximum tile can only intersect by 1 base
        Each tile is checked to ensure they intersect with then region
        Let 
        '*' region base
        '-' flanking base
        '>' forward primer base
        '<' reverse primer base
        '=' a tile base
        Let 
        max_tile = 5, 
        max_primer_len = 4, 
        region length = 8
        S = start_search_index = max_primer + (max_tile -1) -1
        E = end_search_index =  max_primer + 2*(max_tile-1) + region_len -1 
        R = region_end_index = end_search_index - (max_tile -1)
        p = pos = current position in the Dp loop
        s = tile start index = p - tile size +1
                0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5   # pos_memo indexing
        0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3  # ref_seq indexing 
        > > > > - - - - * * * * * * * * - - - - < < < <
                        S       s = = R p     E
        """
        max_primer_len = self.primer_length + self.primer_length_var
        start_search_index = max_primer_len + (self.max_tile -1)
        end_search_index = (max_primer_len 
                            + 2 * (self.max_tile -1) 
                            + self.region_length -1) #+1 for exclusive end
        for pos in range(start_search_index, end_search_index +1):
            # Initialise the best_result to their null values
            # The final best_result is the best choice 
            # at this pos given the best choice
            # from all previous position
            # best_result 
            #   = (best_score, 
            #      best_tile, 
            #      best_overlap, 
            #      tile_count for this solution starting at this pos, 
            #      best_forward, 
            #      best_reverse)
            best_result = (0, 0, 0, 0, None, None)
            for t_size in self.tile_sizes:
                region_end_index = end_search_index - (self.max_tile -1)
                tile_start_index = pos - t_size +1
                if tile_start_index <= region_end_index: #if still overlap
                    tile = (tile_start_index, t_size)
                    # select best f and r primers of this tile
                    tile_score, f_primer, r_primer = self.best_primers_in_tile(tile)
                    for overlap in range(self.allowed_overlap):
                        # access the score of all the tiles before
                        # given the choice of current tile
                        #              p    
                        #     **********
                        #           -===    # t_size = 4, overlap = 1
                        #           x       # x = prefix pos index
                        prefix_pos = tile_start_index -1 + overlap
                        prefix_pos_in_memo = prefix_pos - max_primer_len
                        prefix_score = self.pos_memo[prefix_pos_in_memo][0]
                        prefix_tile_count = self.pos_memo[prefix_pos_in_memo][3]
                        total_score = prefix_score + tile_score 
                        tile_count = prefix_tile_count +1
                        if total_score > best_result[0]:
                            best_result = (total_score,
                                           t_size,
                                           overlap,
                                           tile_count,
                                           f_primer,
                                           r_primer)
            position_in_memo = pos - max_primer_len
            self.pos_memo[position_in_memo] = best_result




    def best_primers_in_tile(self, tile):
        '''
        Take a tile in the reference region and look for the best
        forward and reverse primers (best with respect to their length)
        '''
        logging.info('Start choosing best primers in tile %s'%str(tile))
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

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length, 'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = get_primer_seq(self.reference, r_primer)
                r_scores = self.score_primer(r_sequence, direction='r')
                self.primer_memo[r_primer] = r_scores

            if f_scores[0] > best_f_score:
                best_f_score = f_scores[0]
                best_f = f_primer
            if r_scores[0] > best_r_score:
                best_r_score = r_scores[0]
                best_r = r_primer
        total_score = best_f_score + best_r_score
        best_result = (total_score, best_f, best_r)
        logging.info('Finished choosing best primers in tile')
        logging.info('Chosen primer at tile %s Score:%s Forward:%s Reverse:%s'
                     %(str(tile), best_result[0], best_result[1], best_result[2])
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


def pick_primer_set(searcher, score_correction_exponent):
    '''
    Look through starting from the end of pos_memo, pick the best scored
    "starting position" and then trace back using the tile_size and overlap
    chosen at each position to recover the optimal primer set.

    '''
    if not filter(lambda x: x != (0, 0, 0, 0, None, None), searcher.pos_memo):
        logging.info('Warning the search has not been done')
        return
    # look from the end of pos_memo until it hits
    # the end of the region and determine the best starting position
    best_reverse_start = None # this will be a negative integer
    best_start_score = 0
    for pos in range(-1, -searcher.max_tile -1, -1):
        pos_info = searcher.pos_memo[pos]
        pos_score, pos_tile_count = pos_info[0], pos_info[3]
        adjusted_score = adjust_score(pos_score, pos_tile_count, score_correction_exponent)
        if adjusted_score > best_start_score:
            best_start_score = adjusted_score
            best_reverse_start = pos
            sys.stderr.write('pos:%s, tile: %s, score: %s, avg: %s\n'
                    %(pos, pos_tile_count, pos_score, adjusted_score))
    sys.stderr.write('region_len:%s, best_start_score: %s, best_pos: %s\n\n'
            %(searcher.region_length, best_start_score, best_reverse_start))

    region_size_remaining = searcher.region_length \
                            + searcher.max_tile -1\
                            #+ best_reverse_start
                             # reverse_start is negative

    position = best_reverse_start
    # while we are still in the region
    # notice the definition of "region" is loosen to include
    # the overhang of the last tile
    while abs(position) <= region_size_remaining:
        (pos_score, tile_size,
         overlap, tile_count, f_primer, r_primer) = searcher.pos_memo[position]
        assert searcher.pos_memo[len(searcher.pos_memo) + position] == searcher.pos_memo[position]
        yield f_primer, r_primer

        # Beware of tile_size, overlap = 0 case, an infinite loop result
        if tile_size == 0 and overlap == 0:
            logging.info("There is error and this will go into an infinite loop")
            raise RuntimeError('''Attempt to shift by 0 position, 
                    this will go into infinite loop''')

        position -= (tile_size - overlap)
        searcher.aux_data['tiles'].append(tile_size)
        searcher.aux_data['overlap'].append(overlap)



def adjust_score(score, tile_count, exponent = 1):
    """
    Adjust total score of a tiling pattern by the number 
    of tiles it has so that the algorithm would not be 
    biased towards having more tiles.
    """
    return score / (float(tile_count) ** exponent)

def write_aux(aux_writer, searcher):
    """
    This functions takes a Dp_search object write
    the auxiliary information to a file 
    via a specified csv.writer object
    """
    tiles = searcher.aux_data['tiles']
    overlap = searcher.aux_data['overlap']
    num_primers = len(searcher.primer_memo)
    output = [searcher.chrom, searcher.region_start,
              searcher.region_end, searcher.region_length,
              searcher.time_taken, mean(tiles),
              len(tiles), std_deviation(tiles),
              mean(overlap), num_primers]
    aux_writer.writerow(output)


def write_primer(out_writer, searcher, score_correction_exponent):
    '''
    This function take in a Dp_search object, and a csv.writer
    object and write all the primers generated by 
    pick_primer_set() (a generator that yields a pair
    of forward and reverse primers together with their data) 
    to file.
    '''
    region_start = searcher.region_start
    region_length = searcher.region_length
    tile_number = 0
    for f_primer, r_primer in pick_primer_set(searcher, score_correction_exponent):
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


def regions_ref_seqs_generator(bedfile, user_inputs):
    max_tile = user_inputs.tiles[1]
    max_primer_len = user_inputs.primer_length + user_inputs.primer_length_var
    fa_path = user_inputs.fa
    regions_and_ref_seqs = get_regions_and_ref_seq(obtain_regions(bedfile), 
                                                   max_tile, 
                                                   max_primer_len, 
                                                   fa_path)
    for region_ref_seq in regions_and_ref_seqs:
        yield region_ref_seq

def obtain_regions(bedfile):
    """
    bedfile specifies the path to a particular BED file
    The function parse the BED file and return a dictionary
    with chromosome as keys and a list of regions in the particular
    chromosome as values.
    region is of the form:
        (name, start position, end position)
    name :: <chromosome>_<first csv value of 4 column in BED file>
    start and end position are 0-based index of the region.
    eg:
    regions_dict = 
    {'chr1':[('chr1_something', 1234, 1890),('chr1_somethingelse',2345,2890)],
     'chrX':[('chrX_thing', 5678,5999)]}
    """
    with open(bedfile) as file:
        reader = csv.reader(file, delimiter='\t')
        region_dict = {}
        # Obtain a list of regions from bedfile
        for line in reader:
            chrom, start, end, name = line[:4]
            start, end = map(int, [start, end])
            # Use first comma separated value in the 
            # name field and the chromosome as
            # the region name
            bed_specified_name = name.split(',')[0].replace('_','-')
            # The name is potentially empty
            # replace '_' by '-' so that all '_' occurance is 
            # controlled by the program itself instead of the user.
            region_name = chrom + '_' + bed_specified_name
            if chrom in region_dict:
                region_dict[chrom]. append((region_name, start, end))
            else:
                region_dict[chrom] = [(region_name, start, end)]
    return region_dict

def get_regions_and_ref_seq(region_dict, max_tile, max_primer_len, fa_path):
    """
    Let S = starting index of region from the start of 
    the chromosome and let '-' denote bases outside of region and
    '*' denote bases in region.
    We require S >= max_tile + max_primer_len -1
    The edge case is shown below:
                               S
           0 1 2 3 4 5 6 7 8 9 10
         | - - - - - - - - - - * * * * * * * * * * *
                     ===========    # max_tile = 6, need 1 base intersect
           -------->                # max_primer_len = 5
    Similar condition is imposed for the other end of the chromosome.

    Also, we require that the distant between 2 regions, L to satisfy
        L >= 2 * (max_tile -1) + max_primer_len +1
    so that the primers between the regions will never overlap 
    completely. (DO WE WANT NOT OVERLAP AT ALL??)
    The edge case is shown below:
                  <-- 
            ======    ======        
                   -->
        *****--------------******
    """
    regions_and_reference = []
    for chrom, region_list in region_dict.items():
        chrom_file_path = os.path.join(fa_path, chrom + '.fa')
        chrom_sequence = SeqIO.read(chrom_file_path, 'fasta').seq
        chrom_len = len(chrom_sequence)
        region_list.sort(key = lambda x:x[1])
        num_region = len(region_list)

        # Give warning if region is too close to the chromosome boundary
        # Only check the first and last region in a given chromosome
        first_region = region_list[0]
        last_region = region_list[-1]
        # Notice that first and last region might be the same if
        # there is only one region in this particular chromosome.
        min_margin = max_tile + max_primer_len -1
        if first_region[1] < min_margin:
            logging.info("""WARNING: %s is too close to left boundary.\n
                    Require at least %i margin\nbut actual margin is:%i
                    """%(first_region, min_margin, first_region[1]))
        if chrom_len - last_region[2] < min_margin:
            logging.info("""WARNING: %s is too close to right boundary
                    at %i.\n
                    Require at least %i margin\nbut actual margin is:%i\n
                    """%(last_region,
                        chrom_len-1,
                        min_margin,
                        last_region[1]))
        # Merge regions which are too close together
        region_index = 0
        while region_index < num_region:
            num_merged = 1 # Merging A with itself count as 1
            current_region = region_list[region_index]
            while region_index + num_merged < num_region:
                next_region = region_list[region_index+num_merged]
                distance = next_region[1] - current_region[2]
                # Give warning for repeated or overlapping regions
                if distance < 0:
                    logging.info("""WARNING: repeated or overlapping region:
                    %s and %s"""%(str(current_region), str(next_region)))
                min_separation = 2 * (max_tile -1) + max_primer_len +1
                if distance < min_separation:
                    new_name = _merge_name(current_region[0], next_region[0])
                    # new name joins the bedspecified name of both regions
                    # if the names are the same, only one of them is used
                    new_start = current_region[1]
                    new_end = next_region[2]
                    logging.info("MERGED: %s and %s\ndistance: %i"
                                 %(current_region, next_region, distance))
                    current_region = (new_name, new_start, new_end)
                    logging.info("NEW REGION: %s"%str(current_region))
                    num_merged +=1
                else:
                    break

                # obtain reference by slicing chrom sequence
                # Reference = union of 
                #       1) Region
                #       2) Max_tile_size -1 + max_primer_len on both side
                # reference is obtained here to avoid reading seq again later.
                #    ----->=====       =====<----- 
                #              *********
            ref_start = current_region[1] - max_tile - max_primer_len +1
            ref_end = current_region[2] + max_tile + max_primer_len
            ref_seq = str(chrom_sequence[ref_start:ref_end])
            regions_and_reference.append((current_region, ref_seq))
            region_index += num_merged
    return regions_and_reference

                            
def _merge_name(name1, name2):
    """
    Name1 and name2 are of the form
        <chromosome>_<name from bedfile>
    This function will return
        <chromomsome>_<name from bedfile 1>_<name from bedfile 2>
    if the names are the same, only one of the will be used, ie
    we return
        <chromosome>_<name from bedfile>
    """
    name1 = name1.split('_')
    name2 = name2.split('_')
    chrom = name1[0]
    if not (name1[0] == name2[0] == chrom):
        logging.info('WARNING: inconsistent chromosome naming')
    new_name = chrom + '_' + '_'.join(set(name1[1:] + name2[1:]))
    return new_name



#def handle_bedfile(bedfile):
#    file = open(bedfile)
#    reader = csv.reader(file,delimiter='\t')
#    bed_row = []
#    for line in reader:
#        chromo, start, end, name = line 
#        start, end = map(int, [start, end])
#        # we dont want underscore character in the bed specified naming
#        # to interfere with our own naming convention
#        bed_specified_name = name.split(',')[0].replace('_','')
#        region_name = chromo + '_' + bed_specified_name 
#        # last term is potentially [''], ie a list of single empty string
#        bed_row.append((region_name, start, end))
#    return bed_row

