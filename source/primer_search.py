import logging
import sys
import os
import time

from utils import rev_complement


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
        region_end_index = end_search_index - (self.max_tile -1)
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
            best_result = (-10000, 0, 0, 0, None, None)
            for t_size in self.tile_sizes:
                tile_start_index = pos - t_size +1
                if tile_start_index <= region_end_index: #if still overlap
                    tile = (tile_start_index, t_size)
                    # select best f and r primers of this tile
                    tile_score, f_primer, r_primer = self.best_primers_in_tile(tile)
                    for overlap in range(self.allowed_overlap +1):
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
        #logging.info('Start choosing best primers in tile %s'%str(tile))
        best_f_score = -10000
        best_r_score = -10000
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
                f_scores = self.score_primer(f_sequence)
                self.primer_memo[f_primer] = f_scores

            # deal with reverse primer
            r_primer = (tile_end +1, this_primer_length, 'r')
            if r_primer in self.primer_memo:
                r_scores = self.primer_memo[r_primer]
            else:
                r_sequence = get_primer_seq(self.reference, r_primer)
                r_scores = self.score_primer(r_sequence)
                self.primer_memo[r_primer] = r_scores

            if f_scores[0] > best_f_score:
                best_f_score = f_scores[0]
                best_f = f_primer
            if r_scores[0] > best_r_score:
                best_r_score = r_scores[0]
                best_r = r_primer
        total_score = best_f_score + best_r_score
        best_result = (total_score, best_f, best_r)
        #logging.info('Finished choosing best primers in tile')
        #logging.info('Chosen primer at tile %s Score:%s Forward:%s Reverse:%s'
        #             %(str(tile), best_result[0], best_result[1], best_result[2]))
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
        return seq #rev_complement(seq)


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
    best_start_score = -10000
    for pos in range(-1, -searcher.max_tile -1, -1):
        pos_info = searcher.pos_memo[pos]
        pos_score, pos_tile_count = pos_info[0], pos_info[3]
        adjusted_score = adjust_score(pos_score, pos_tile_count, score_correction_exponent)
        if adjusted_score > best_start_score:
            best_start_score = adjusted_score
            best_reverse_start = pos
            logging.info('pos:%s, tile: %s, score: %s, avg: %s'
                    %(pos, pos_tile_count, pos_score, adjusted_score))
    logging.info('region_len:%s, best_start_score: %s, best_pos: %s'
            %(searcher.region_length, best_start_score, best_reverse_start))

    region_size_remaining = searcher.region_length \
                            + searcher.max_tile -1\
                            #+ best_reverse_start
                             # reverse_start is negative

    position = best_reverse_start
    empirical_tile_count = 0
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
        empirical_tile_count +=1 
        searcher.aux_data['tiles'].append(tile_size)
        searcher.aux_data['overlap'].append(overlap)
    logging.info("Empirical Tile count: %i"%empirical_tile_count)



def adjust_score(score, tile_count, exponent = 1):
    """
    Adjust total score of a tiling pattern by the number 
    of tiles it has so that the algorithm would not be 
    biased towards having more tiles.
    """
    return score / (float(tile_count) ** exponent)

