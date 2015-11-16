"""
(This module is still under active development)
Description:
    This module contains the class "Score" which implements the scoring
    method for measuring the goodness of a particular primer sequence
    given the users' commandline parameters. 
    The method will be used within the dynamic programming loop in
    the primer_design.py module, thus the scoring is strictly 
    a function of primer (nucleic acid sequence) to a real number, 
    for instance,
    primer positions and current chosen pool of primers would not affect
    the score of a primer.
    
    Nov 2015:
    This module now apply statistical method to produce a scoring
    function.
    Statistical scoring method 1:
        Assume normal variate for all variables (individual criteria
        such as tm, gc, etc .. and also the combined weighted score).
        It would score all possible primers from the regions provided
        in the bedfile and calculate the mean and variance.
        The score would be "number of std deviation from the mean"

    Statistical scoring method 2:
        Built empirical distribution of all scoring criteria by
        scoring all possible primers.
        The score of a primer would be given by its percentile
        in the empirical distribution.

"""


import logging
from math import log
from time import time

from Bio.SeqUtils.MeltingTemp import (Tm_NN, 
                                      Tm_GC, 
                                      Tm_Wallace, 
                                      DNA_NN1, 
                                      DNA_NN2, 
                                      DNA_NN3, 
                                      DNA_NN4)
from utils import weighted_num_complement, mean, std, normalised_distance 
from specificity import Specificity

# Nearest neighbor thermodynamic data for Tm prediction
# from a Biopython module
# Refer to 
# http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html
THERMO_TABLES = [DNA_NN1, DNA_NN2, DNA_NN3, DNA_NN4]

# Functions that map sequence to their relevant raw data.
def raw_entropy(seq):
    """Calculate Shanon's entropy for seq using natural logarithm"""
    seq = seq.upper()
    length = float(len(seq))
    prob = [seq.count(base)/length for base in set(seq)]
    return - sum( p * log(p) for p in prob)

def raw_max_hairpin(seq, direction='f'):
    """
    Hairpin prediction algorithm. 
    Currently returning just the maximum number of bond that is
    sterically possible for seq together with its attached heel sequence.
    Only considering 1 loop of size 3 occuring at where the sequence
    folds onto itself.
    GC bond is weighted more than AT bond based on user_inputs.gc_weight
    """
    if direction == 'f':
        seq = self.sense_heel + seq
    elif direction == 'r':
        seq = self.antisense_heel + seq
    else:
        logging.info('Warning: primer_direction not specified for hairpin prediction')
        raise ValueError
    length = len(seq)
    score = 0
    loop_size = 3
    start_index = loop_size + 2 
    # Start at 5:
    #           * * * * * * * * *    (Original sequence and index)
    #           1 2 3 4 5 6 7 8 9

    #             3 2 1 
    #             * * *              (First possible loop)
    #         4 *   | |
    #             * * * * * 
    #             5 6 7 8 9 
    # Thus top = 2 1 
    #      bottom = 6 7 8 9
    for index in range(start_index,length):
        # only consider loop of size 3
        # smaller loop are assummed to be sterically impossible
        # and case where top and bottom with
        # only 1 bp is not considered
        top = seq[:index-loop_size][::-1]
        bottom = seq[index:]
        new_score = weighted_num_complement(top,bottom,
                                            gc_weight = self.gc_weight)
        if new_score > score:
            score = new_score
    return score
        

def raw_gc(seq):
    """ 
    Return the gc content as a fraction of gc_count/length of seq
    """
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/float(len(seq))

def raw_gc_clamp(seq):
    """
    Return 1 if sequence has a gc_clamp else return 0
    """
    end_with = seq[-1].upper()
    return float(end_with == 'G' or end_with == 'C')

def raw_run(seq):
    """
    Return the maximum number of run we occuring in 
    the sequence
    eg: raw_run(AAAGGGGGATCTATTCT) == raw_run(GGGGG) == 5
    """
    seq = seq.upper()
    max_run = 1
    run = 1
    current_base = seq[0]
    for base in seq[1:]:
        if base == current_base:
            run +=1
        else:
            current_base = base
            if run > max_run:
                max_run = run
            run = 1
    if run > max_run:
        max_run = run
    return max_run

#################### Utilities ##################################

def get_all_primers(regions_and_ref_seqs, primer_length, primer_length_var):
    """
    regions_and_ref_seqs: [(regions coords, ref seq)]
    Return a dictionary of all possbile primers given
    the reference sequences
    """
    all_primers = {}
    num_distinct_primers = 0
    for region, ref_seq in regions_and_ref_seqs:
        for pos in range(len(ref_seq)):
            for var in range(-primer_length_var, primer_length_var +1):
                primer = ref_seq[pos: pos + primer_length + var].upper()
                if len(primer) == primer_length + var:
                    if primer not in all_primers:
                        num_distinct_primers += 1
                        all_primers[primer] = []
    logging.info("Number of all possible primers: %i"%num_distinct_primers)
    return all_primers




####################


class Score(object):
    """
    This class handles the user inputs and produces method for primer
    scoring. 
    It takes in user commandline input from running primer_design,
    extract relevant information about users' preference, specification
    of reaction condition, and then initialise a scoring function
    based on them. 
    """
    def __init__(self, user_inputs, regions_and_ref_seqs):
        """
        Initialise Score class base on user_inputs.
        After initialisation, score_func and tm_func would be
        defined and can be used in the DP_search process.
        """
        self.score_func = None
        # Initialise relevant informations for scoring 
        self.tm_func = None
        self.user_inputs = user_inputs
        # Look up data tables for scoring
        self.combined_data = {'specificity':[], 'tm':[], 'entropy':[], 'gc':[]}
        self.primer_data = get_all_primers(regions_and_ref_seqs,
                                           user_inputs.primer_length,
                                           user_inputs.primer_length_var)
        #self.primer_data = {}
        # Parse user_inputs and assign them to relevant attribute
        self._define_functions(user_inputs, regions_and_ref_seqs)
        logging.info('Finished initialising Score object')

    def _define_functions(self, user_inputs, regions_and_ref_seqs):
        """
        This private method is used to initialised the Score class.
        It parses user_inputs, decide which tm_prediciton function to
        use and then provide the relevant parameters from user_inputs
        to the tm_prediciton functions.
        It also parse user choice of scoring and then build (define)
        a function that scores primers base on statistical methods.
        """
        # Define tm prediction function to be used.
        tm_func_choice = user_inputs.tm_func
        if tm_func_choice == "Tm_NN":
            self.tm_func = lambda seq : Tm_NN(seq, 
                                              Na = user_inputs.Na,
                                              K = user_inputs.K,
                                              Mg = user_inputs.Mg,
                                              Tris = user_inputs.Tris,
                                              dNTPs = user_inputs.dNTPs,
                                              nn_table=THERMO_TABLES[user_inputs.NN_table-1],
                                              dnac1 = user_inputs.dnac1,
                                              dnac2 = user_inputs.dnac2,
                                              saltcorr = user_inputs.saltcorr)
        elif tm_func_choice == "Tm_GC":
            self.tm_func = lambda seq : Tm_GC(seq,
                                              Na = user_inputs.Na,
                                              K = user_inputs.K,
                                              Mg = user_inputs.Mg,
                                              Tris = user_inputs.Tris,
                                              dNTPs = user_inputs.dNTPs,
                                              saltcorr = user_inputs.saltcorr)

        elif tm_func_choice == "Tm_Wallace":
            self.tm_func = Tm_Wallace

        # Define scoring function by building a lookup table for it
        # first we score all possible primers
        
        self._collect_all_raw_scores()
        #self._score_all_primers(regions_and_ref_seqs,
        #                        user_inputs.primer_length,
        #                        user_inputs.primer_length_var)
        if user_inputs.score_func == 'normal':
            self.score_func = self._create_normal_score_func()
        elif user_inputs.score_func == 'empirical':
            self.score_func = self._create_empirical_score_func()


    def _collect_all_raw_scores(self):
        before = time()
        self.primer_data = Specificity(self.user_inputs.fa, self.primer_data).primers
        for primer in self.primer_data:
            tm, entropy, gc = (self.tm_func(primer),
                               raw_entropy(primer),
                               raw_gc(primer))
            specificity = self.primer_data[primer][0]
            self.primer_data[primer] += [tm, entropy, gc]
            self.combined_data['tm'].append(tm)
            self.combined_data['entropy'].append(entropy)
            self.combined_data['gc'].append(gc)
            self.combined_data['specificity'].append(specificity)
        logging.info("Finished collecting all raw data in %s"
                     %(time() - before))






#    def _score_all_primers(self, regions_ref_seqs, primer_length, primer_length_var):
#        """
#        After __init__ had finished parsing all user_inputs from commandline,
#        it will call this method to obtain all possible primers from the 
#        reference sequences provided by regions_ref_seqs and simultaneously
#        compute their relevant properties and store them in a dictionary
#        that map Sequence to tuples of their data.
#        """
#        before = time()
#        num_distinct_primers = 0
#        for region, ref_seq in regions_ref_seqs:
#            for pos in range(len(ref_seq)):
#                for var in range(-primer_length_var, primer_length_var +1):
#                    primer = ref_seq[pos: pos + primer_length + var].upper()
#                    if len(primer) == primer_length + var:
#                        if primer not in self.primer_data:
#                            num_distinct_primers += 1
#                            tm, entropy, gc = (self.tm_func(primer), 
#                                               raw_entropy(primer), 
#                                               raw_gc(primer))
#                            self.primer_data[primer] = [tm, entropy, gc]
#                            self.combined_data['tm'].append(tm)
#                            self.combined_data['entropy'].append(entropy)
#                            self.combined_data['gc'].append(gc)
#        if not num_distinct_primers > 0:
#            logging.info("Warning: there is no primers to be scored")
#            raise RuntimeError
#        elif not all([len(lst) == num_distinct_primers for lst in self.combined_data.values()]):
#            logging.info("Warning: The number of data point doesn't equal num primers")
#            raise RuntimeError
#        logging.info("Finished scoring all %i primers in %s seconds"
#                     %(num_distinct_primers, time() - before))
#
    

    def _create_normal_score_func(self):
        """
        Return a function object which assume normality in
        distribution of primer data.
        """
        before = time()
        mean_std = {}
        for feature, data in self.combined_data.items():
            mean_std[feature] = (mean(data), std(data))
        logging.info("%s"%(str(mean_std)))
        for primer in self.primer_data:
            specificity, tm, entropy, gc = self.primer_data[primer]
            scores = [normalised_distance(specificity, mean_std['specificity']),
                      normalised_distance(tm, mean_std['tm']),
                      normalised_distance(entropy, mean_std['entropy']),
                      normalised_distance(gc, mean_std['gc'])]
            primer_score = sum(scores)/ float(len(scores))
            self.primer_data[primer] = [primer_score, specificity, tm, entropy, gc]
        scoring_function= lambda seq: self.primer_data[seq.upper()]
        logging.info("Finished creating scoring function in %s"%(time() - before))
        return scoring_function

    def _build_quantile_lookup(self):
        """
        This private method is called to build a lookup table
        for each scoring feature we are interested in. The lookup
        table is of the form:
            {feature: {scoring data : percentile} }
        This data structure is intended to speed up the process
        of scoring primers by looking up its percentile score.
        """
        # First we build the empirical distribution of each 
        #!!! features of the primers (currently they're just
        #!!! tm, entropy and gc). This amount to building up
        #!!! a lookup list which is a function that map a 
        #!!! value to the quantile it is in.
        quantile_lookup = {'specificity':{}, 'tm':{}, 'entropy':{}, 'gc':{}}
        for feature, data in self.combined_data.items():
            data.sort()
            num_data = float(len(data))
            current_val = data[0]
            quantile_lookup[feature][current_val] = 0.0 - 50.0
            # There's always at least one data guaranteed by _score_all_primers
            for index in range(1, len(data)):
                next_val = data[index]
                if next_val != current_val:
                    quantile_lookup[feature][next_val] = index * 100.0/ num_data -50.0
                    current_val = next_val
        return quantile_lookup

    def _create_empirical_score_func(self):
        """
        Returns a function object that maps primers to
        their empirical percentile score.
        """
        before = time()
        quantile_lookup = self._build_quantile_lookup()
        missing_data_count = 0
        for primer in self.primer_data:
            try:
                specificity, tm, entropy, gc = self.primer_data[primer]
                scores = [quantile_lookup['specificity'][specificity],
                          quantile_lookup['tm'][tm],
                          quantile_lookup['entropy'][entropy],
                          quantile_lookup['gc'][gc]]
                primer_score = sum(scores) / float(len(scores))
                self.primer_data[primer] = [primer_score, specificity, tm, entropy, gc]
            except:
                missing_data_count += 1
        scoring_function = lambda seq: self.primer_data[seq.upper()]
        logging.info("Finished creating scoring function in %s"%(time() - before))
        return scoring_function
