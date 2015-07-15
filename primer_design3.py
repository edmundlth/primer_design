"""

Reimplementation of primer_design1.py to rethink correctness
and improve readibitlity

"""

from utils import score, visualise
from time import time
from random import randint

import sys
sys.setrecursionlimit(10000)

PRIMER_LENGTH = 20

MEMO = {}
MEMO_CALLS = 0
CALLS = 0

def design(template, tile_sizes, pos, template_length):
    best_combination = []
    best_score = -10000000
    for size in tile_sizes:
        for overlap in range(1, size +1):
            prefix_pos = pos - (size - overlap)
            suffix_pos = pos + overlap
            suffix_length = template_length - overlap
            
            first_tile = (prefix_pos, size)
            first_tile_score = score(template, first_tile, PRIMER_LENGTH)

            best_suffix_design = suffix_design(template,
                                               tile_sizes,
                                               suffix_pos,
                                               suffix_length)
            
            new_score = first_tile_score + best_suffix_design[1]
            if new_score > best_score:
                best_combination = [prefix_pos , size] +\
                                   best_suffix_design[0][1:]
                best_score = new_score
    print("Best primer score = %f" %best_score)
    return best_combination
                                               

def suffix_design(template, tile_sizes, suffix_pos, length):
    global CALLS, MEMO_CALLS
    CALLS +=1
    
    if length <= 0:
        return ([],0)
    # Memoization
    if (suffix_pos, length) in MEMO:
        MEMO_CALLS +=1
        return MEMO[(suffix_pos,length)]

    best_tiling = [suffix_pos]
    best_score = -10000
    for size in tile_sizes:
        prefix_score = score(template, (suffix_pos, size), PRIMER_LENGTH)
        if size >= length:
            new_tiling = [suffix_pos,size]
            new_score = prefix_score
        else:
            # recursion: solve same problem for the suffix of suffix
            rest = suffix_design(template,
                                 tile_sizes,
                                 suffix_pos + size,
                                 length - size,)
            new_tiling = [suffix_pos,size] + rest[0][1:]
            new_score = prefix_score + rest[1]
        if new_score > best_score:
            best_tiling = new_tiling
            best_score = new_score
    result = (best_tiling, best_score)
    MEMO[ (suffix_pos, length) ] = result
    return result















if __name__ == "__main__":
    bases = ['A','T','G','C']
    
    template1 = 'AAATGCACGAAAAATCGTGGCGTTTTATTTGTGCTAGTCGTGCGTGAAAATTCGTCCC'
    template_length1 = len(template1) - 11
    tile_sizes1 = [8,9,10]
    
    template2 = "AATGCGTCAGTTGAC"
    template_length2 = len(template2) -6
    tile_sizes2 = [3,4,5]
    
    tile_sizes3 = [100,101,102,103,104,105,106,107,108,109,110]
    template_length3 = 500
    template3 = ''
    for b in [bases[randint(0,3)] for i in range(template_length3)]:
        template3 +=b
    template_length3 -= max(tile_sizes3)

    template = template3
    tile_sizes = tile_sizes3
    template_length = template_length3
    pos = max(tile_sizes)
    
    before = time()
    tiling = design(template, tile_sizes, pos, template_length)
    after = time()

##    visualise(template,tiling, pos, pos + template_length, PRIMER_LENGTH )
    print(tiling)
    print("CALLS = %i "%CALLS)
    print("MEMO call = %i"%MEMO_CALLS)
    print("MEMO length = number of new CALLS = %i" %len(MEMO))
    print("Time_taken = %f" %(after - before))

