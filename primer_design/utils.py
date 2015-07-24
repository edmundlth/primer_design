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


from score import rev_complement

def visualise_tile(template, tile, f_r_lens):
    template = str(template)
    start = tile[0]
    end = tile[0] + tile[1]
    f_len, r_len = f_r_lens
    relevant_template = template[start - f_len:end + r_len]
    f_primer = relevant_template[:f_len]
    f_bond = '|' * f_len
    r_primer = ' ' * (end - start + f_len)\
                 + relevant_template[- r_len:]
    r_bond = ' ' * (end - start + f_len) + '|' * r_len
    #print(f_primer)
    #print(f_bond)
    #print(relevant_template)
    #print(r_bond)
    #print(r_primer)
    return (f_primer, f_bond, relevant_template, r_bond, r_primer)



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


def bed_coords(file_name):
    bed_file = open(file_name)
    start_len_coords = []
    for line in bed_file:
        coords = line.split('\t')[:3]
        start = int(coords[1])
        length = int(coords[2]) - start
        start_len_coords.append((coords[0],start,length))
    bed_file.close()
    return start_len_coords

