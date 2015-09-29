"""
Objective of this module is to provide the tool to generate 
statistics, graphs and annotated files for the output of
primer_design tool for Hi-Plex. 

The main function should take in (as command-line arguments)
the result files for primer_designa and the name of the desired
output files to store the program.
"""

import os
import sys
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import itertools
import csv
import swalign
from time import time
import random
############################ parse_arg() and main() #######################
DEFAULT_OUTFILE = 'primer_analysis.txt'
DEFAULT_DIRECTORY = '.'
DEFAULT_SW_SCORES = [2, -1, -1, -1, 0.0]

def parse_args():
    parser = ArgumentParser(description='''Analysis tools for primer
                                          design output''')
    parser.add_argument('--dir', metavar='DIR', type=str,
                        default=DEFAULT_DIRECTORY,
                        help='''A string specifying the path to the
                        directory containing the primer design result
                        files. Default to current directory''')
    parser.add_argument('--primer_file', metavar='PRIMER_FILE',
                        type=str, required=True,
                        help='''The file name of the result file
                        containing all the primer's raw data''')
    parser.add_argument('--aux_file', metavar='AUXFILE',
                        type=str, required=True,
                        help='''The file name of the file containing
                        axiliary data from primer design''')
    parser.add_argument('--outfile', metavar="OUTFILE",
                        type=str, 
                        default=DEFAULT_OUTFILE,
                        help='''The name of the output file of this
                        program. Default to %s'''%DEFAULT_OUTFILE)
    parser.add_argument('--pairplot', metavar='COLUMN_NAME',nargs='*',
                        type=str,
                        default=[''],
                        help='''string value to specify the columns in 
                        primer data frame to be plotted in a pairplot''')
    parser.add_argument('--scores', metavar='SCORE', nargs=5, type=float,
                        default=DEFAULT_SW_SCORES,
                        help='''List of floats that specify,
                        (match score, mismatch score, gap_penalty,
                        gap_extension_penalty, gap_extension_decay)''')
    return parser.parse_args()


def main():
    # parse user inputs and initialise swaligner object
    user_inputs = parse_args()
    (match, mismatch, gap_penalty,
    gap_extension_penalty, gap_extension_decay) = user_inputs.scores
    sw_aligner = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(match,mismatch),
                                        gap_penalty=gap_penalty,
                                        gap_extension_penalty=gap_extension_penalty,
                                        gap_extension_decay=gap_extension_decay,
                                        prefer_gap_runs=True,
                                        verbose=False,
                                        globalalign=False,
                                        wildcard=None,
                                        full_query=False)
    # produce data frames from the provided files
    primer_file = os.path.join(user_inputs.dir, user_inputs.primer_file)
    aux_file = os.path.join(user_inputs.dir, user_inputs.aux_file)
    data = Primer_Data_Frames(primer_file, aux_file)

    # annotate data frames with dimer scores
    # including swalign.score, 3' end scores, match, mismatch
    annotate_dimer(data.primer_df, sw_aligner)
    sys.stderr.write(str(data.primer_df))
    with open(user_inputs.outfile,'w') as outfile:
        outfile.write(str(data.primer_df.describe()) + '\n')
        outfile.write(str(data.aux_df.describe())+'\n')
    if user_inputs.pairplot != ['']:
        fig_name = user_inputs.outfile.split('.')[0] + '.png'
        sns.pairplot(data.primer_df[user_inputs.pairplot])
        plt.savefig(fig_name, format='png')

#################### Smith-Waterman Alignment #########################
def align_to_pool(primer, pool, sw_aligner):
    """
    Align a single query primer to a pool of primer
    and return the best scored alignment 
    (a swalign.align object). The query primer
    will be reverse complemented before alignment.
    """
    best_align = None
    best_score = 0
    query = rev_complement(primer)
    for reference in pool:
        align = sw_aligner.align(reference, query)
        if align.score > best_score:
            best_score = align.score
            best_align = align
    return best_align

def annotate_dimer(primer_df, sw_aligner):
    """
    Given the primer data frame (containing a column call sequence),
    produce another 3 columns in the data frame (mutating it):
        align_score: best sw_aligner.score among the pool
        3'_end: produced using the cigar string
        3'_end_normalised
        match_mismatch: number of (mis)matches in the best alignment
    """
    pool = list(primer_df['sequence'])
    length = len(pool)
    align_scores = []
    end_matches = []
    match_mismatch = []
    dimer_scores = []
    for query in pool:
        best_align = align_to_pool(query, pool, sw_aligner)
        a_score = best_align.score
        end = _get_end_matches(best_align.cigar)
        mm = (best_align.matches, best_align.mismatches)
        d_score = dimer_scoring(a_score, end, mm, query)

        align_scores.append(a_score)
        end_matches.append(end)
        match_mismatch.append(mm)
        dimer_scores.append(d_score)
    primer_df['align_score'] =  align_scores
    primer_df['end_matches'] = end_matches
    primer_df['match_mismatch'] = match_mismatch
    primer_df['dimer_score'] = dimer_scores

# !!! still in development, the precise scoring function is not determined
def dimer_scoring(align_score, end_matches, match_mismatch, sequence):
    """
    Taking in alignment data and the sequence itself,
    produce a number between 0-100 as a score for primer dimer
    """
    length = float(len(sequence))
    gc_content = _gc_content(sequence)
    return 100 * (align_score/(2 * length) 
                  + end_matches/length 
                  + gc_content
                  + match_mismatch[0]/length
                  - match_mismatch[1]/length)


def _gc_content(sequence):
    sequence = sequence.upper()
    length = float(len(sequence))
    gc_count = sequence.count('C') + sequence.count('G')
    return gc_count / length
def _get_end_matches(cigar):
    """
    Given cigar input, return the number
    of initial matches.
    eg.
    given [(3, 'M'), (2, 'I'), (3, 'M')]
    return 3
    """
    if cigar and cigar[0][1] == 'M':
        return cigar[0][0]
    else:
        sys.stderr.write('WARNING: Cigar is empty')
        return 0




def annotate_dimer_random_pool(primer_df, sw_aligner, align_proportion=0.5):
    """
    Mutate the given panda.DataFrame object so that it
    will have a new column called 'dimer' which has the
    data about the dimer score of the sequences in 
    the 'sequence' column aligned to all the sequences
    in the pool.
    The score is the best sw_aligner.score from a random
    pool (align_proportion of the total pool)
    """
    pool = list(primer_df['sequence'])
    length = len(pool)
    dimer_score = []
    for query in pool:
        random.shuffle(pool)
        random_draw = pool[:int(length * align_proportion)]
        dimer_score.append(align_to_pool(query, random_draw, sw_aligner).score)
    primer_df['dimer'] = dimer_score


def gen_primers_fasta(primer_file, primer_fa=None):
    """
    A function that takes the result file containing
    all the primer sequence and raw data and write all
    primer sequence to a single fasta file together with
    their name as sequence header.
    eg. 
    primer_file = trial_out.tsv:
    name          start  end   seq       Tm   ...
    chr16_PALB2_f1 23546  23552  ATGCTTG  30 ...
    chr16_PALB@_r1 23646  23652  GTGTAGC  40 ....

    will be written to trial_out_primers.fa as:
    >chr16_PALB2_f1
    ATGCTTG
    >chr16_PALB2_r1
    GTGTAGC
    """
    if not primer_fa:
        primer_fa = os.path.basename(primer_file).split('.')[0] + '_primers.fa'
    with open(primer_file) as data_file:
        with open(primer_fa,'w') as outfile:
            reader = csv.reader(data_file, delimiter='\t')
            next(reader)
            primers = [(line[0],line[3]) for line in reader]
            for name, seq in primers:
                outfile.write('>' + name + '\n' \
                              + seq + '\n')


def rev_complement(seq):
    """
    Return the reverse compliment of a sequence
    in capital letters. The sequence have to be
    composed of the alphabets "ATGC" or their
    lower case.
    """
    bases = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
    return ''.join([bases[b] for b in seq.upper()])[::-1]






#################### Primer Data Frames ##############################

class Primer_Data_Frames(object):
    def __init__(self, primer_file=None, aux_file=None):
        self.primer_df = pd.read_csv(primer_file, sep='\t')
        self.aux_df = pd.read_csv(aux_file, sep='\t')

        #self.primer_numeric_df = get_numeric_df(self.primer_df)
        #self.aux_numeric_df = get_numeric_df(self.aux_df)
        #self.kmer_table = None
        #self.k = 0

#    def annotate_dimer(self):
#        if self.kmer_table == None:
#            print('Warning: kmer counting have not been done')
#        else:
#            self.primer_df['dimer'] = pd.Series(
#                                                map(self.dimer_score, 
#                                                    self.primer_df['sequence'])
#                                               )
#
#
#    def dimer_score(self, seq, k,gc_weight=1.5, gc_cutoff=0.3):
#        if seq != '':
#            match_to = complement(seq[-k:])
#            rev_match_to = match_to[::-1]
#        if match_to in self.kmer_table:
#            match = self.kmer_table[match_to]
#        else:
#            match = 0
#
#        if rev_match_to in self.kmer_table:
#            rev_match = self.kmer_table[rev_match_to]
#        else:
#            rev_match = 0
#
#        return (match + rev_match) * gc_factor(seq, gc_weight, gc_cutoff)
#
#def make_kmer_table(k, bedfilename, fasta_dir, outfilename=None):
#    return generate_kmer_
#    self.kmer_table = generate_kmer_table(k, bedfilename, fasta_dir, outfilename)
#
#
#def gc_factor(seq,gc_weight=1.5, gc_cutoff=0.3):
#    seq = seq.upper()
#    gc_fraction = (seq.count('G') + seq.count('C'))/float(len(seq))
#    if gc_fraction > gc_cutoff:
#        return 1 + gc_fraction * gc_weight
#    else:
#        return 0
#
#def complement(seq):
#    bases = {'A':'T','G':'C','T':'A','C':'G'}
#    return ''.join([bases[b] for b in seq.upper()])
#
#def get_numeric_df(df):
#    '''Return a new data frame which consist of the all the 
#    columns of the original data frame which has numeric type
#    '''
#    droplist = []
#    for col in df.columns:
#        column_type = df[col].dtype
#        if column_type not in [np.dtype('float64'), np.dtype('int64')]:
#            droplist.append(col)
#    return df.drop(droplist,axis=1)

##################### kmer counting ####################################
#
#def generate_kmer_table(k, bedfilename, fasta_dir, outfilename=None):
#    if not outfilename:
#        outfilename = bedfilename.split('.')[0] + '_region.fa'
#    generate_ref_file(bedfilename, outfilename, fasta_dir)
#    ref_seqs = get_ref_seqs(outfilename)
#    kmer_table = kmer_count(ref_seqs, k)
#    return kmer_table 
#
#def kmer_count(ref_seqs, k):
#    count_dict = {}
#    for seq in ref_seqs:
#        length = len(seq)
#        for index in range(length - k +1):
#            kmer = seq[index:index + k]
#            if kmer in count_dict:
#                count_dict[kmer] +=1
#            else:
#                count_dict[kmer] = 1
#    return count_dict
#
#
#def get_ref_seqs(ref_seqs_filename):
#    ref_seqs = [str(reader.seq).upper() for reader in SeqIO.parse(ref_seqs_filename,'fasta')]
#    return ref_seqs
#
#def generate_ref_file(bedfilename, outfilename,fasta_dir, slack = 150):
#    coords = handle_bed(bedfilename)
#    coords = sorted(coords)
#    with open(outfilename,'w') as outfile:
#        for chrom_group in itertools.groupby(coords, lambda x:x[0]):
#            chrom = chrom_group[0]
#            fasta_filename = os.path.join(fasta_dir, '%s.fa'%chrom)
#            with open(fasta_filename) as f_file:
#                ref_seq = SeqIO.read(f_file,'fasta').seq
#                length_ref = len(ref_seq)
#                for coord in chrom_group[1]:
#                    start = max(0, coord[1] - slack)
#                    end = min(length_ref, coord[2] + slack)
#                    seq = str(ref_seq[start:end])
#                    seq_name = '_'.join(map(str,coord))
#                    outfile.write('>' + seq_name + '\n')
#                    outfile.write(seq + '\n')
#
#
#def handle_bed(bedfilename, delimiter='\t'):
#    with open(bedfilename) as bedfile:
#        reader = csv.reader(bedfile, delimiter=delimiter)
#        coords = [tuple([line[0]] + map(int,line[1:3])) for line in reader]
#        return coords

if __name__ == '__main__':
    main()
