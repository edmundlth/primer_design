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
NO_PRIMER_FILE = 'None'
NO_AUX_FILE = 'None'
NO_FA = 'None'
DEFAULT_SW_SCORES = [2, -1, -1, -1, 0.0]

def parse_args():
    parser = ArgumentParser(description='''Analysis tools for primer
                                          design output''')
    parser.add_argument('--fa', metavar='FA_FILE', type=str,
                        default=NO_FA,
                        help='''A string specifying the fasta file
                        containing all the primer sequence''')
    parser.add_argument('--primer_file', metavar='PRIMER_FILE',
                        type=str, default=NO_PRIMER_FILE,
                        help='''The file name of the result file
                        containing all the primer's raw data''')
    parser.add_argument('--aux_file', metavar='AUXFILE',
                        type=str, default=NO_AUX_FILE,
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
    if user_inputs.primer_file != NO_PRIMER_FILE:
        # produce data frames from the provided files
        primer_file = user_inputs.primer_file
#        primer_file = os.path.join(user_inputs.dir, user_inputs.primer_file)
        primer_df = pd.read_csv(primer_file, sep='\t')
        
        # annotate data frames with dimer scores
        # including swalign.score, 3' end scores, match, mismatch
        annotate_dimer(primer_df, sw_aligner)
        sys.stderr.write(str(primer_df.describe()))
        with open(user_inputs.outfile, 'w') as outfile:
            primer_df.to_csv(outfile, sep='\t')
        if user_inputs.pairplot != ['']:
            fig_name = user_inputs.outfile.split('.')[0] + '.png'
            sns.pairplot(primer_df[user_inputs.pairplot])
            plt.savefig(fig_name, format='png')
    if user_inputs.aux_file != NO_AUX_FILE:
        aux_file = user_inputs.aux_file
#        aux_file = os.path.join(user_inputs.dir, user_inputs.aux_file)
        aux_df = pd.read_csv(aux_file, sep='\t')
        sys.stderr.write(str(aux_df.describe()))


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
    matches = []
    mismatches = []
    dimer_scores = []
    for query in pool:
        best_align = align_to_pool(query, pool, sw_aligner)
        a_score = best_align.score
        end = _get_end_matches(best_align.cigar)
        mm = (best_align.matches, best_align.mismatches)
        d_score = dimer_scoring(a_score, end, mm, query)

        align_scores.append(a_score)
        end_matches.append(end)
        matches.append(mm[0])
        mismatches.append(mm[1])
        dimer_scores.append(d_score)
    primer_df['align_score'] =  align_scores
    primer_df['end_matches'] = end_matches
    primer_df['matches'] = matches
    primer_df['mismatches'] = mismatches
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


################### Auxiliary data processing ######################

def wastage_annotate(aux_df):
    pass



#################### Primer Data Frames ##############################
#
#class Primer_Data_Frames(object):
#    def __init__(self, primer_file=None, aux_file=None):
#        self.primer_df = pd.read_csv(primer_file, sep='\t')
#        self.aux_df = pd.read_csv(aux_file, sep='\t')


if __name__ == '__main__':
    main()
