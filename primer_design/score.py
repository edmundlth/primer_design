from utils import stability, weighted_num_complement
from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_GC, Tm_Wallace
from math import log


class Score(object):
    def __init__(self, user_inputs):
        self.score_func = None
        self.tm_func = None
        self.target_tm = None
        self.tm_underachieve_weight = None
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
        self._handle_user_inputs(user_inputs)


    def _handle_user_inputs(self, user_inputs):
        """
        This private method is called during the class initiation
        to assign all relevant user input to corresponding 
        class attributes.
        """
        self.target_tm = user_inputs.target_tm
        self.tm_underachieve_weight = user_inputs.tm_underachieve
        self.sense_heel = user_inputs.sense_heel.upper()
        self.antisense_heel = user_inputs.antisense_heel.upper()

        (self.Na, self.Mg, self.K,self.dNTPs,
         self.Tris,self.dnac1,self.dnac2) = user_inputs.conc

        self.saltcorr = user_inputs.saltcorr
        self.gc_weight = user_inputs.gc_weight

        score_choice = user_inputs.score_func
        if score_choice == 'score_Lp':
            self.score_func = self.score_Lp
        elif score_choice == 'score_linear':
            self.score_func = self.score_linear


        tm_func_choice = user_inputs.tm_func
        if tm_func_choice == "Tm_NN":
            self.tm_func = lambda seq : Tm_NN(seq, Na = self.Na,
                                                   K = self.K,
                                                   Mg = self.Mg,
                                                   Tris = self.Tris,
                                                   dNTPs = self.dNTPs,
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
        return self.tm_func(seq)
    
    def raw_entropy(self, seq):
        seq = seq.upper()
        length = float(len(seq))
        prob = [seq.count(base)/length for base in set(seq)]
        return - sum( p * log(p) for p in prob)

    def raw_max_hairpin(self, seq, direction):
        if direction == 'f':
            seq = self.sense_heel + seq
        elif direction == 'r':
            seq = self.antisense_heel + seq
        else:
            print('Warning: primer_direction not specified for hairpin prediction')
            raise ValueError
        length = len(seq)
        score = 0
        for i in range(5,length):
            # only consider loop of size 3
            # smaller loop are assummed to be sterically impossible
            # and case where top and bottom with
            # only 1 bp is not considered
            top = seq[:i-3][::-1]
            bottom = seq[i:]
            new_score = weighted_num_complement(top,bottom,
                                                gc_weight = self.gc_weight)
            if new_score > score:
                score = new_score
        return score
            

    def raw_gc(self, seq):
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        return gc_count/float(len(seq))

    def raw_gc_clamp(self,seq):
        end_with = seq[-1].upper()
        if end_with == 'G' or end_with == 'C':
            return 1
        else:
            return 0

    def raw_run(self, seq):
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
        target = self.target_tm # looking up attribute expensive?
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
        # max entropy is log(4) = 1.386294
        return entropy * 100 / 1.386294

    def normalise_gc(self, gc_fraction):
        return gc_fraction * 100

    def normalise_run(self, run, length):
        return 100 * (1 - float(run) / length)

    def normalise_hairpin(self, hairpin, length, gc_fraction):
        # maximum_hairpin = AT is counted as 1 while GC are 
        # weighted according to user specifed gc weight
        maximum_hairpin = length * ( gc_fraction * self.gc_weight + 
                                        (1 - gc_fraction))
        return 100 * (1 - hairpin / maximum_hairpin)

    def normalise_gc_clamp(self, gc_clamp):
        if gc_clamp:
            return 100
        else:
            return 0



    # Scalarisation of gathered data from the sequence to get
    # a single number as the score for the sequence
    def score_Lp(self, seq, direction, p = 2):
        p = float(p)
        length = len(seq)
        tm = self.raw_tm(seq)
        entropy = self.raw_entropy(seq)
        hairpin = self.raw_max_hairpin(seq,direction)
        gc_fraction = self.raw_gc(seq)
        gc_clamp = self. raw_gc_clamp(seq)
        run = self.raw_run(seq)

        # score = (sum of of the p power of normalised scores) ** (1/p)
        score = sum( 
                    map( lambda x: x**p, 
                            [self.normalise_tm(tm),
                            self.normalise_entropy(entropy),
                            self.normalise_hairpin(hairpin, length, gc_fraction),
                            self.normalise_gc(gc_fraction),
                            self.normalise_gc_clamp(gc_clamp),
                            self.normalise_run(run, length)]
                        )
                    ) ** (1/p)
        return (score, tm, entropy, hairpin, gc_fraction, gc_clamp, run)

    def score_linear(self, seq):
        pass








