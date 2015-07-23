"""
Utilities functions that will be called by the primer_design.py program
Author : Edmund Lau (elau1@student.unimelb.edu.au)
------------------------------------------------------------------------
Description:
Provide the utilitis which enable
- visualisation
- file handling
------------------------------------------------------------------------


"""
from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_Wallace, Tm_GC
import math
import itertools
import random







def visualise(template,tiling,left,right,primer_lengths):
    """
    This function generate a visual of the position of the primers base
    on the tiling scheme inputed.
    Primer pair of a tile is simply the sequences of length = primer_length
    which flank the tile.
    Forward primers are shown above the sequence at the right position
    Reverse primers are shown below the sequence at the right position
    The region of interest is in upper case, lower case otherwise
    """
    
    template = template[:left].lower()+\
               template[left:right+1].upper()+\
               template[right+1:].lower()
    upper = ' '* (tiling[0] - primer_lengths[0][0]) +\
            template[tiling[0]- primer_lengths[0][0] : tiling[0]]
    lower = ' ' * (tiling[0]+ primer_lengths[0][1] + 1)
    template = template[:tiling[0]] + '|' + template[tiling[0]:]
    
    last_tile = tiling[0]
    for i in range(1,len(primer_lengths)):
        inc = tiling[i]
        f_primer_length = primer_lengths[i][0]
        r_primer_length = primer_lengths[i][1]
        last_tile += inc +1
        lower += ' '* (inc - r_primer_length+1) + \
                 template[last_tile: last_tile+r_primer_length]
        upper += ' '* (inc -f_primer_length +1) + \
                 template[last_tile - f_primer_length: last_tile]
        template = template[:last_tile] + '|'+template[last_tile:]
    upper = upper[:-f_primer_length]

    chunk_size = 60
    for i in range(0,len(template),chunk_size):
        print('   ',upper[i: i+chunk_size],'   ')
        print("5'>",template[i: i+chunk_size],">3'")
        print('   ',lower[i: i+chunk_size],'   ')

        


