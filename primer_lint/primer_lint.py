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
DEFAULT_OUTFILE = 'primer_analysis.tsv'
NO_AUX_FILE = 'None'
NO_FA = 'None'
NO_BOXPLOT = ['']
NO_PAIRPLOT = ['']
NO_COLORING = ''
DEFAULT_MATCH_SCORE = 2
DEFAULT_MISMATCH_SCORE = -1
DEFAULT_GAP_PENALTY = -1
DEFAULT_GAP_EXTENSION_PENALTY = -1
DEFAULT_GAP_EXTENSION_DECAY = 0.0
DEFAULT_END_MATCH_SCORE = 2



################### Danny's Dimer scoring ############################

from primer_dimer11_new import Score_gaps

#!!! Danny's algorithm scoring parameters (not included in the module
#!!! commandline parameters yet. Need manual modifications.
SCORE_KEY = {'lastx':20, 'mm1':8, 'mm2':7, 'mm3':6,
             'gc':4, 'at':2, 'mm':4, 'w':'yes', 'wx':1.0,
             'gapping':'yes', 'tile':3, 'tf':1, 'loopl':4,
             'loopu':10, 'stop':7, 'agg_thresh':0}
d_score_dimer = lambda seq1, seq2: Score_gaps('', 
                                              seq1.upper(), 
                                              '', 
                                              seq2.upper(), 
                                              SCORE_KEY)[0]

def one_to_all_score(primer, pool):
    """
    produce a dimer score (higher is bad) for primer against the
    whole pool. using Danny's dimer prediction algorithm.
    """
    best_score = 0
    query = rev_complement(primer)[:20]
    for reference in pool:
        score = d_score_dimer(query, reference)
        if score > best_score:
            best_score = score
    return best_score


#######################################################################
def parse_args():
    parser = ArgumentParser(description='''Analysis tools for primer
                                          design output''')
    parser.add_argument('--fa', metavar='FA_FILE', type=str,
                        default=NO_FA,
                        help='''A string specifying the fasta file
                        containing all the primer sequence''')
    parser.add_argument('--primer_file', metavar='PRIMER_FILE',
                        type=str, required=True,
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
                        default=NO_PAIRPLOT,
                        help='''string values to specify the columns in 
                        primer data frame to be plotted in a pairplot''')
    parser.add_argument('--color_with', metavar='CLASS', type=str,
                        default=NO_COLORING,
                        help='''The name of the column of column score
                        with to classify (binary) the data points in
                        the pairplot to above average and below avg''')
    parser.add_argument('--boxplot', metavar='COLUMN_NAME', nargs='*',
                        type=str,
                        default=NO_BOXPLOT,
                        help='''string values to specify the columns of
                        primer data frame which are to be boxplotted.''')
    parser.add_argument('--match', metavar='SCORE', type=int,
                        default=DEFAULT_MATCH_SCORE,
                        help='''Match score for Smith Waterman alignment''')
    parser.add_argument('--mismatch', metavar='SCORE', type=int,
                        default=DEFAULT_MISMATCH_SCORE,
                        help='''Mismatch score for Smith Waterman alignment''')
    parser.add_argument('--gap_penalty', metavar='SCORE', type=int,
                        default=DEFAULT_GAP_PENALTY,
                        help='''Gap penalty for Smith Waterman alignment''')
    parser.add_argument('--gap_extension_penalty', metavar='SCORE', type=int,
                        default=DEFAULT_GAP_EXTENSION_PENALTY,
                        help='''Gap extension penalty for Smith Waterman alignment''')
    parser.add_argument('--gap_extension_decay', metavar='SCORE', type=float,
                        default=DEFAULT_GAP_EXTENSION_DECAY,
                        help='''Gap extension decay for Smith Waterman alignment''')
    parser.add_argument('--end_match', metavar='SCORE', type=int,
                        default=DEFAULT_END_MATCH_SCORE,
                        help='''The scoring weight for 3' end matches''')
    parser.add_argument('-d', action='store_true',
                        help='''Use Danny's dimer prediciton algorithm if specified''')
    return parser.parse_args()


def main():
    # parse user inputs and initialise swaligner object
    user_inputs = parse_args()
    sw_aligner = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(user_inputs.match,
                                                                        user_inputs.mismatch),
                                        gap_penalty=user_inputs.gap_penalty,
                                        gap_extension_penalty=user_inputs.gap_extension_penalty,
                                        gap_extension_decay=user_inputs.gap_extension_decay,
                                        prefer_gap_runs=True,
                                        verbose=False,
                                        globalalign=False,
                                        wildcard=None,
                                        full_query=False)
    # produce data frames from the provided files
    primer_file = user_inputs.primer_file
    primer_df = pd.read_csv(primer_file, sep='\t')
    
    # annotate data frames with dimer scores
    if user_inputs.d: #!!! use Danny's algorithm
        pool = primer_df['sequence']
        primer_df['dimer_score'] = map(lambda primer: one_to_all_score(primer, pool), pool)
    else:
        annotate_dimer(primer_df, sw_aligner, user_inputs.end_match)
    # annotate each primer with it's dimer score rank
    # higher rank, better dimer, worse primer.
    primer_df['dimer_rank'] = rank(primer_df['dimer_score'])

    # Annotate the off target rank as well if the the off target
    # count is provided in column which has name
    # either 'off' or 'off target'
    # !!! this handles fergus and kplex mapping data files.
    try:
        primer_df['off_rank'] = rank(primer_df['off'])
    except KeyError:
        try:
            primer_df['off_rank'] = rank(primer_df['off target'])
        except KeyError:
           pass 


    sys.stderr.write(str(primer_df.describe()))
    with open(user_inputs.outfile, 'w') as outfile:
        primer_df.to_csv(outfile, sep='\t')


    fig_name = user_inputs.outfile.split('.')[0]
    if user_inputs.pairplot != NO_PAIRPLOT:
        if user_inputs.color_with == NO_COLORING:
            coloring = None
        else: # if user wants coloring, produce the binary class first
            coloring = user_inputs.color_with
            hue_column = 'above_avg_' + coloring
            primer_df[hue_column] = above_threshold(list(primer_df[coloring]))
        sns.pairplot(data=primer_df, vars=user_inputs.pairplot, hue=hue_column)
        plt.savefig(fig_name + '.png', format='png')
    if user_inputs.boxplot != NO_BOXPLOT:
        for column_name in user_inputs.boxplot:
            plt.figure()
            sns.boxplot(primer_df[column_name])
            plt.savefig(fig_name+'_box_%s.png'%column_name, format='png')

    if user_inputs.aux_file != NO_AUX_FILE:
        aux_file = user_inputs.aux_file
        aux_df = pd.read_csv(aux_file, sep='\t')
        wastage_annotate(aux_df)
        sys.stderr.write(str(aux_df.describe()))
        total = sum(aux_df['region_length'])
        tiled = sum([count * mean 
                     for count, mean in 
                     zip(aux_df['tile_count'], aux_df['mean_tile'])])
        wasted = tiled - total
        assert wasted >= 0
        sys.stderr.write('\n\nTotal region length: %i\n'%total)
        sys.stderr.write('Total tiled length: %.1f\n'%tiled)
        sys.stderr.write('Wasted tiling length: %.1f\n'%wasted)





#################### Smith-Waterman Alignment #########################
def align_to_pool(primer, pool, sw_aligner, end_match_weight):
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
        end_matches = _get_end_matches(align.cigar)
        this_score = align.score + end_match_weight * end_matches
        if this_score > best_score:
            best_score = this_score
            best_align = align
    return best_align

def annotate_dimer(primer_df, sw_aligner, end_match_weight):
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
        best_align = align_to_pool(query, pool, sw_aligner, end_match_weight)
        a_score = best_align.score
        end = _get_end_matches(best_align.cigar)
        mm = (best_align.matches, best_align.mismatches)
        d_score = dimer_scoring(a_score, end, mm, query, end_match_weight)

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
def dimer_scoring(align_score, end_matches, match_mismatch, sequence, end_match_weight):
    """
    Taking in alignment data and the sequence itself,
    produce a number between 0-100 as a score for primer dimer
    """
    length = float(len(sequence))
    gc_content = _gc_content(sequence)
    return 100 * (align_score/(2 * length) 
                  + end_matches * end_match_weight/length 
                  + gc_content
                  + match_mismatch[0]/length
                  - match_mismatch[1]/length)


def rank(score_list):
    """
    Produces the 0-base ranks of a list of scores.
    The higher the rank the higher the score.
    eg:
    >>> rank([1, 2, 4, 3])
    [3, 2, 0, 1]
    """
    indexing = sorted(zip(score_list,range(len(score_list))), reverse = False)
    result = [0 for i in range(len(score_list))]
    for index in range(len(score_list)):
        score, original_index = indexing[index]
        result[original_index] = index
    return result

def above_threshold(score_list, threshold=None):
    """
    Produces a list of 1's and 0's representing
    the truth value at the ith position in
    the score_list whether score_list[i] > threshold.
    1 -> True
    0 -> False
    If threshold = None, the mean is used as threshold
    """
    if threshold == None:
        threshold = mean(score_list)
    return [int(score > threshold) for score in score_list]

def mean(scores):
    """
    Return the mean of scores
    """
    if scores:
        return sum(scores)/float(len(scores))
    else:
        return 0

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



#!!! Unused
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
    aux_df['covered_length'] = map(lambda x: x[1] * x[0], 
                                   zip(aux_df['tile_count'], 
                                       aux_df['mean_tile']))
    aux_df['wastage'] = map(lambda x: x[1] - x[0],
                            zip(aux_df['region_length'],
                                aux_df['covered_length']))


if __name__ == '__main__':
    main()
