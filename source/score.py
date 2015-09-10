"""
(This module is still under active development)
Description:
    This module contains the class Score which implements the scoring
    method for measuring the goodness of a particular primer sequence
    given the users' commandline parameters. 
    The method will be used within the dynamic programming loop in
    the primer_design.py module, thus the scoring is strictly 
    a function of primer (nucleic acid sequence) to a real number, 
    for instance,
    primer positions and current chosen pool of primers would not affect
    the score of a primer.

"""


from utils import weighted_num_complement
from Bio.SeqUtils.MeltingTemp import (Tm_NN, 
                                      Tm_GC, 
                                      Tm_Wallace, 
                                      DNA_NN1, 
                                      DNA_NN2, 
                                      DNA_NN3, 
                                      DNA_NN4)
from math import log
import logging

# Nearest neighbor thermodynamic data for Tm prediction
# from a Biopython module
# Refer to 
# http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html
THERMO_TABLES = [DNA_NN1, DNA_NN2, DNA_NN3, DNA_NN4]

class Score(object):
    """
    This class handles the user inputs and produces method for primer
    scoring. 
    It takes in user commandline input from running primer_design,
    extract relevant information about users' preference, specification
    of reaction condition, and then initialise a scoring function
    based on them. 
    """
    def __init__(self, user_inputs):
        logging.info('Initialising Score object')
        self.score_func = None

        # Initialise relevant informations for scoring 
        self.tm_func = None
        self.score_weights = None
        self.target_tm = None
        self.tm_underachieve_weight = None
        self.NN_table = None
        self.Na = None
        self.Mg = None
        self.K = None
        self.dNTPs = None
        self.Tris = None
        self.saltcorr = None
        self.sense_heel = None
        self.antisense_heel = None
        self.gc_weight = None
        self.dnac1 = None
        self.dnac2 = None
        # Parse user_inputs and assign them to relevant attribute
        self._handle_user_inputs(user_inputs)
        logging.info('Finished initialising Score object')


    def _handle_user_inputs(self, user_inputs):
        """
        This private method is called during the class initiation
        to assign all relevant user input to corresponding 
        class attributes.
        """
        self.target_tm = user_inputs.target_tm
        self.tm_underachieve_weight = user_inputs.tm_underachieve
        assert self.tm_underachieve >= 1.0
        self.NN_table = user_inputs.NN_table
        self.sense_heel = user_inputs.sense_heel.upper()
        self.antisense_heel = user_inputs.antisense_heel.upper()

        (self.Na, self.Mg, self.K,self.dNTPs,
         self.Tris,self.dnac1,self.dnac2) = user_inputs.conc

        self.saltcorr = user_inputs.saltcorr
        self.gc_weight = user_inputs.gc_weight
        assert self.gc_weight >= 1.0

        score_choice = user_inputs.score_func
        self.score_weights = user_inputs.score_weights
#        if score_choice == 'score_Lp':
        self.score_func = (
                lambda seq, direction: self.score_Lp(seq, direction, p = score_choice)
                          )
#        elif score_choice == 'score_linear':
#            self.score_func = self.score_linear
#            self.linear_weights = user_inputs.linear_weights

        tm_func_choice = user_inputs.tm_func
        if tm_func_choice == "Tm_NN":
            self.tm_func = lambda seq : Tm_NN(seq, Na = self.Na,
                                                   K = self.K,
                                                   Mg = self.Mg,
                                                   Tris = self.Tris,
                                                   dNTPs = self.dNTPs,
                                                   nn_table=THERMO_TABLES[self.NN_table-1],
                                                   dnac1 = self.dnac1,
                                                   dnac2 = self.dnac2,
                                                   saltcorr = self.saltcorr)
        elif tm_func_choice == "Tm_GC":
            self.tm_func = lambda seq : Tm_GC(seq,Na = self.Na,
                                                   K = self.K,
                                                   Mg = self.Mg,
                                                   Tris = self.Tris,
                                                   dNTPs = self.dNTPs,
                                                   saltcorr = self.saltcorr)

        elif tm_func_choice == "Tm_Wallace":
            self.tm_func = Tm_Wallace
        else:
            print("Warning: Tm prediction function not chosen, default to Tm_NN")


    # Functions that map sequence to their relevant raw data.
    def raw_tm(self, seq):
        """ Use the chosen Tm function which is aware of user specified
        reaction conditions such as reactants concentration, 
        salt correction method and nearest neighbour data, 
        to predict the melting temperature of seq"""
        return self.tm_func(seq)
    
    def raw_entropy(self, seq):
        """Calculate Shanon's entropy for seq using natural logarithm"""
        seq = seq.upper()
        length = float(len(seq))
        prob = [seq.count(base)/length for base in set(seq)]
        return - sum( p * log(p) for p in prob)

    def raw_max_hairpin(self, seq, direction):
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
            print('Warning: primer_direction not specified for hairpin prediction')
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
            

    def raw_gc(self, seq):
        """ 
        Return the gc content as a fraction of gc_count/length of seq
        """
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        return gc_count/float(len(seq))

    def raw_gc_clamp(self,seq):
        """
        Return 1 if sequence has a gc_clamp else return 0
        """
        end_with = seq[-1].upper()
        return float(end_with == 'G' or end_with == 'C')

    def raw_run(self, seq):
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


    # Normalisation of the raw data

    def normalise_tm(self, tm):
        """
        Return a value between 0 - 100 as a normalised
        score for a given tm value.
        The further tm is from the target_tm specified by
        user, the lower the score.
        However, tm < target_tm is scored more harshly
        based on user specified tm_underachieve_weight
        which is a value >= 1.0.
        """
        target = self.target_tm
        tm_difference = tm - target
        lower = 20.0
        upper = 80.0
        range = (upper - target) - (lower - target)
        normalised = tm_difference / range
        if normalised < 0:
            return 100*(1 + self.tm_underachieve_weight * normalised)
        elif normalised > 0:
            return 100 * (1 -  normalised)
        else:
            return 100

    def normalise_entropy(self, entropy):
        """
        Normalise entropy to a value between 0 - 100
        entropy is maximum when all four of ATGC
        has equal probability of occuring, ie
        pA = pT = pG = pC = 1/4, hence
        max_entropy = -log(1/4) = log(4) ~= 1.386294
        """
        return entropy * 100 / 1.386294

    def normalise_gc(self, gc_fraction):
        """
        Normalise gc_fraction to a value between
        0-100
        """
        return gc_fraction * 100

    def normalise_run(self, run, length):
        """
        Normalise run as a proportion of 
        length to a value between 0 - 100
        """
        return 100 * (1 - float(run) / length)

    def normalise_hairpin(self, hairpin, length, gc_fraction):
        """
        Normalise hairpin to a value between 0 - 100 according
        to the fraction hairpin/ max_hairpin, where
        max_hairpin is the hairpin value when half of the
        sequence bind to the other half.
        # maximum_hairpin = AT is counted as 1 while GC are 
        # weighted according to user specifed gc weight
        """
        maximum_hairpin = 0.5 * length * ( gc_fraction * self.gc_weight + 
                                        (1 - gc_fraction))
        return 100 * (1 - hairpin / maximum_hairpin)

    def normalise_gc_clamp(self, gc_clamp):
        '''
        Return 100 is there is a gc_clamp,
        zero otherwise'''
        return float(bool(gc_clamp)) * 100



    # Scalarisation of gathered data from the sequence to get
    # a single number as the score for the sequence
    def score_Lp(self, seq, direction, p = 2):
        """
        Return the score of the nucleic acid sequence: seq,
        together will all the raw values of its properties.
        Score is calculated using the Lp-norm:
           Score = sum( wi * ni),
              where wi = weight of the i-th property of seq
                    ni = normalised score of the i-th property
        """
        p = float(p)
        length = len(seq)
        tm = self.raw_tm(seq)
        entropy = self.raw_entropy(seq)
        hairpin = self.raw_max_hairpin(seq,direction)
        gc_fraction = self.raw_gc(seq)
        gc_clamp = self. raw_gc_clamp(seq)
        run = self.raw_run(seq)

        
        normalised_scores = [self.normalise_tm(tm),
                            self.normalise_entropy(entropy),
                            self.normalise_hairpin(hairpin, length, gc_fraction),
                            self.normalise_gc(gc_fraction),
                            self.normalise_gc_clamp(gc_clamp),
                            self.normalise_run(run, length)]


        # score = (sum of of the p power of normalised scores) ** (1/p)
        score = sum( 
                    map( lambda x: x**p,
                         [normalised * weight
                          for normalised, weight in 
                          zip(normalised_scores, self.score_weights)]
                        )
                    ) ** (1/p)
        return (score, tm, entropy, hairpin, gc_fraction, gc_clamp, run)
