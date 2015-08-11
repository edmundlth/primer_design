


        
def dp_search(background, region_coord, user_inputs):
    global POS_MEMO

    region_start, region_end = region_coord
    region_len = region_end - region_start
    min_tile, max_tile = user_inputs.tiles
    background_len = region_len + 2 * (max_tile -1)
    tile_sizes = [min_tile + extend 
                  for extend in range(max_tile - min_tile +1)]
    allowed_overlap = user_inputs.overlap

    best_score = 0
    best_f = None
    best_r = None
    best_overlap = 0
    for pos in range(region_start, region_end +1):
        for t_size in tile_sizes:
            tile = (pos - , t_size)
            suffix = best_primers_in_tile(background,
                                          tile,
                                          user_inputs)
            suffix_score, f_primer, r_primer = suffix
            for overlap in range(allowed_overlap):
                prefix_pos = pos - region_start + max_tile - t_size + overlap
                prefix_score = POS_MEMO[prefix_pos][0]
                total_score = prefix_score + suffix_score

                if total_score > best_score:
                    best_score = total_score
                    best_f = f_primer
                    best_r = r_primer
                    best_overlap = overlap
        POS_MEMO[pos - region_start + max_tile] = (best_score, best_f, best_r)
    return



def best_primers_in_tile(background, tile, user_inputs):
    global PRIMER_MEMO, SCORE_PRIMER
    primer_length = user_input.primer_length
    length_var = user_input.length_var
    tile_start, tile_size = tile
    tile_end = tile_start + tile_size -1

    best_f = None
    best_r = None
    best_f_score = 0
    best_r_score = 0
    for vary in range(-length_var, length_var +1):
        this_primer_length = primer_length + vary
        
        # deal with forward primer
        f_primer = (tile_start -1, this_primer_length)
        if f_primer in PRIMER_MEMO:
            f_scores = PRIMER_MEMO[f_primer]
        else:
            f_sequence = background[tile_start - this_primer_length: tile_start]
            f_scores = SCORE_PRIMER(f_sequence)
            PRIMER_MEMO[f_primer] = f_scores

        # deal with reverse primer
        r_primer = (tile_end +1, this_primer_length)
        if r_primer in PRIMER_MEMO:
            r_scores = PRIMER_MEMO[r_primer]
        else:
            r_sequence = background[tile_end +1 : tile_end + this_primer_length +1]
            r_scores = SCORE_PRIMER(r_sequence)
            PRIMER_MEMO[r_primer] = r_scores


        if f_scores[0] > best_f_score:
            best_f_score = f_scores[0]
            best_f = f_primer
        if r_scores[0] > best_r_score:
            best_r_score = r_scores[0]
            best_r = r_primer
    total_score = best_f_score + best_r_score

    return (total_score, f_primer, r_primer)
