import sys
import os
import logging
import csv
from Bio import SeqIO

def weighted_num_complement(seq1,seq2,gc_weight = 1.5):
    """Return the number of complementary base of 
    the first min(len_seq1,len_seq2) base pair of seq1 and 2.
    AT binding is counted as 1, GC is counted as gc_weight
    eg: weighted_num_complement('atgc',tcagatt', gc_weight=2)
           A T G C
           |     |
           T C A G T T
        gets a score of 1 + 2 = 3
    """
    weighted_num_complement = 0
    for b1,b2 in zip(seq1,seq2):
        if is_complement(b1,b2):
            b1 = b1.upper()
            if b1 == 'G' or b1 == 'C':
                weighted_num_complement += gc_weight
            else:
                weighted_num_complement +=1                                                 
    return weighted_num_complement
    

def is_complement(base1,base2):
    """
    Return true if base2 is a Watson-Crick pair of base1
    false otherwise
    """
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    if base1 in complement_dic:
        return complement_dic[base1] == base2
    else:
        return False


complement = {'A':'T','G':'C','C':'G','T':'A'}
def rev_complement(seq):
    """
    Return the reverse complement of seq
    """
    seq = seq.upper()
    return ''.join(complement[base] for base in seq)[::-1]

def mean(sequence):
    """
    Calculate the average of sequence.
    If sequence is empty, return 0.0
    """
    if sequence:
        return sum(sequence)/float(len(sequence))
    else:
        return 0.0

def std_deviation(sequence):
    """
    Calculate the standard deviation of a sequence of number
    if sequence is empty, return 0.0
    """
    if sequence:
        x_bar = mean(sequence)
        square_deviation = sum([(x - x_bar)**2 for x in sequence])
        variance = square_deviation / float(len(sequence))
        return variance ** 0.5
    else:
        return 0.0


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
            print("""WARNING: %s is too close to left boundary.\n
                    Require at least %i margin\nbut actual margin is:%i
                    """%(first_region, min_margin, first_region[1]))
            logging.info("""WARNING: %s is too close to left boundary.\n
                    Require at least %i margin\nbut actual margin is:%i
                    """%(first_region, min_margin, first_region[1]))
        if chrom_len - last_region[2] < min_margin:
            print("""WARNING: %s is too close to right boundary
                    at %i.\n
                    Require at least %i margin\nbut actual margin is:%i\n
                    """%(last_region,
                        chrom_len-1,
                        min_margin,
                        last_region[1]))
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
                    print("""WARNING: repeated or overlapping region:
                    %s and %s"""%(str(current_region), str(next_region)))
                    logging.info("""WARNING: repeated or overlapping region:
                    %s and %s"""%(str(current_region), str(next_region)))
                min_separation = 2 * (max_tile -1) + max_primer_len +1
                if distance < min_separation:
                    new_name = _merge_name(current_region[0], next_region[0])
                    # new name joins the bedspecified name of both regions
                    # if the names are the same, only one of them is used
                    new_start = current_region[1]
                    new_end = next_region[2]
                    print("MERGED: %s and %s\ndistance: %i"
                                 %(current_region, next_region, distance))
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
