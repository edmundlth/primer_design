import logging

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
