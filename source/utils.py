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


complement_bases = {'A':'T','G':'C','C':'G','T':'A'}
def rev_complement(seq):
    """
    Return the reverse complement of seq
    """
    return complement(seq)[::-1]

def complement(seq):
    """
    return the complement of seq
    """
    return ''.join(complement_bases[b] for b in seq.upper())

############################ Statistics ########################
def quantile(val, sorted_data):
    """Return the quantile of "val" in the sorted data"""
    return sorted_data.index(val)/float(len(sorted_data))

def percentile_score(val, sorted_data):
    """Return percentile of val in the sorted data"""
    return quantile(val, sorted_data) * 100.0 - 50.0

def mean(data):
    """Calculate the mean of the given data -- a list
    of numbers. If data is empty, 0 is returned"""
    if data:
        return sum(data) / float(len(data))
    return 0.0

def variance(data):
    """Calculate the variance of data. Return 0.0 
    if data is empty"""
    if data:
        mu = mean(data)
        return sum((x - mu)**2 for x in data) / float(len(data))
    return 0.0

def std(data):
    """Return the standard deviation of data"""
    return variance(data) ** 0.5

def normalised_distance(val, mean_std):
    """Return the number (rational number allowed) of 
    standard deviation of val from the mean.
    mean_std == (mean, std)"""
    mu, std = mean_std
    return (val - mu)/ float(std)

