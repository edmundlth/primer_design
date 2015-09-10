import sys

def weighted_num_complement(seq1,seq2,gc_weight = 1.5):
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
    complement_dic = {'A':'T','G':'C','T':'A','C':'G'}
    if base1 in complement_dic:
        return complement_dic[base1] == base2
    else:
        return False


complement = {'A':'T','G':'C','C':'G','T':'A'}
def rev_complement(seq):
    seq = seq.upper()
    return ''.join(complement[base] for base in seq)[::-1]

def mean(sequence):
    if sequence:
        return sum(sequence)/float(len(sequence))
    else:
        return 0.0

def std_deviation(sequence):
    if sequence:
        x_bar = mean(sequence)
        variance = sum([(x - x_bar)**2 for x in sequence])/ float(len(sequence))
        return variance ** 0.5
    else:
        return 0.0

