"""
This module implements Smith-Waterman alignment algorithm
in a primer-dimer detection context.
(although the module would be in a form suitable for 
general sequence alignment too.

The idea is to have a much more costumisable scoring during
alignment to capture the behaviour intended for dimerisation.
The approach is to use the bottom up dynamic programming idea inherent in
the SW alignment and make the objective functions more general.
"""

import sys
from argparse import ArgumentParser


DEFAULT_MATCH = 2.0
DEFAULT_MISMATCH = -2.0
DEFAULT_GAP = -3.0
DEFAULT_GAP_EXT = 0.3
DEFAULT_MISMATCH_EXPONENT = 1.0 # ie no decay
DEFAULT_GAP_EXPONENT = 1.0

MIN_LOOP_SIZE = 3

def parse_args():
    parser = ArgumentParser(description="""Smith Waterman alignment
                                           tool with costumisable gap
                                           scoring for Hiplex primer
                                           design""")
    parser.add_argument('-q', '--query', metavar="QUERY", type=str,
                        required=True, 
                        help='''The query string that will be matched
                        to the reference string specified in --ref''')
    parser.add_argument('-r', '--ref', metavar="REFERENCE", type=str,
                        required=True,
                        help='''The reference string''')
    parser.add_argument('--match', metavar='SCORE', type=float,
                        default=DEFAULT_MATCH,
                        help='''A floating point number specifying
                        the score given to a matching pair''')
    parser.add_argument('--mismatch', metavar='PENALTY', type=float,
                        default=DEFAULT_MISMATCH,
                        help='''A floating point number specifying
                        the penalty (hence a negative number) given
                        to a mismatched pair.''')
    parser.add_argument('--gap', metavar='PENALTY', type=float,
                        default=DEFAULT_GAP,
                        help='''A floating point number specifying
                        the penalty for starting / opening a gap''')
    parser.add_argument('--gap_ext', metavar='GRADIENT', type=float,
                        default=DEFAULT_GAP_EXT,
                        help='''A floating point number specifying
                        the gradient g in the equation:
                           gap penalty = max{0,
                                      gap_open_penalty + 
                                      g * gap_length}
                        if g is a positive number the, the penalty
                        is less harsh for longer gap than short gaps.''')
    parser.add_argument('--costume_score', action='store_true',
                        help='''If specified, the costumised gap scoring
                        scheme will be used instead''')
    parser.add_argument('--mismatch_exponent', metavar='EXPONENT', type=float,
                        default=DEFAULT_MISMATCH_EXPONENT,
                        help='''A floating point number specifying the
                        exponent in the equation:
                            mismatch penalty = max penalty / (distance from 3' end) ** exponent
                        where max penalty = mismatch penalty specified in --mismatch.
                        This will be applyied if --costume_score is specified''')
    parser.add_argument('--gap_exponent', metavar='EXPONENT', type=float,
                        default=DEFAULT_GAP_EXPONENT,
                        help='''A floating point number specifying the
                        exponent in the equation:
                            gap penalty = max penalty / (gap length) ** exponent
                        if gap length < min loop size = 3, the penalty is -inf,
                        if gap length >= min loop size, the above equation is used''')
    parser.add_argument('--verbose', action='store_true',
                        help='''If specified, the decision matrix and alignment string
                        will be printed out''')
    return parser.parse_args()

    



class SWalign(object):
    """
    Smith Waterman alignment object.
    """
    def __init__(self, match=2, 
                 mismatch=-2, gap=-3, gap_ext=0.3,
                 mismatch_exponent=1.0,
                 gap_exponent=1.0,
                 costume_score=False,
                 verbose=False):
        """
        Ultimately the decision matrix would be of the form

        eg: matching query = 'AB' to reference = 'ABC'
           -    A      B      C
        - [(0,T) (0,T) (0,T) (0,T)]
        A [(0,T) (2,M) (1,D) (2,M)]
        B [(0,T) (1,I) (4,M) (1,D)]
        """
        self.match_score = match
        self.mismatch_score = mismatch
        self.gap_penalty = gap
        self.gap_extension_decay = gap_ext
        self.costume_score = costume_score
        self.mismatch_exponent = mismatch_exponent
        self.gap_exponent = gap_exponent
        self.verbose = verbose

        self.num_row = None
        self.num_col = None
        self.decision_matrix = None
        self.query = None
        self.reference = None
        #self.alignment_is_done = False

    def align(self, query, reference):
        """
        Initialise a decision matrix, fill it in, obtain
        optimal alignment and return a score.
        """
        self.query = query
        self.reference = reference
        self.num_row = len(query) +1
        self.num_col = len(reference) +1
        self.decision_matrix = [[(0, 'T') for col in range(self.num_col)]
                                for row in range(self.num_row)]
        for row in range(1, self.num_row):
            for col in range(1, self.num_col):
                self.decision_matrix[row][col] = self.best_decision(row, col)
        alignment = self.get_max_pos_and_score()
        if self.verbose:
            self.print_matrix()
            self.print_alignment()
        return alignment[2]

    def best_decision(self, row, col):
        """
        Should be a function of only the decision matrix
        and the row, col index to fill in.
        It return a tuple: (best_score, decision)
        best_decision = max{(match score, M), (insert score, I), (delete score, D)}
        """
        decisions = [self.match(row, col), 
                     self.insertion(row, col), 
                     self.deletion(row, col), 
                     (0, 'T')] # T = terminate
        return max(decisions, key=lambda x:x[0])


    def match(self, row, col):
        """
        Examine the current state of the decision matrix
        and return the score if the character at
        (row, col) is to be a match or mismatch
        """
        q_base = self.query[row -1]
        r_base = self.reference[col -1]
        prefix_score = self.decision_matrix[row -1][col -1][0]
        if q_base == r_base:
            extension_score = self.match_score
        else:
            if self.costume_score:
                extension_score = self.mismatch_score / (row ** self.mismatch_exponent)
            else:
                extension_score = self.mismatch_score
        return (extension_score + prefix_score, 'M')

    def insertion(self, row, col):
        """
        Return the score 
        """
        all_scores = [self.penalise_gap(length) + \
                      self.decision_matrix[row][col-length][0] 
                      for length in range(1, col + 1)]
        score = max(all_scores)
        return (score, 'I')

    def deletion(self, row, col):
        score = max([self.penalise_gap(length) + \
                     self.decision_matrix[row-length][col][0]
                     for length in range(1, row + 1)])
        return (score, 'D')

    def penalise_gap(self, gap_length):
        if self.costume_score:
            if gap_length < MIN_LOOP_SIZE:
                gap_penalty = -1000
            else:
                gap_penalty = self.gap_penalty / gap_length ** self.gap_exponent
        else:
            gap_penalty = min(0, self.gap_penalty + self.gap_extension_decay * gap_length)
        return gap_penalty


    def print_matrix(self):
        """
        Print the decision_matrix in tabular form
        """
        for row in self.decision_matrix:
            print(row)
        print('\n')

    #!!! Notice that there might be multiple position with the max score, ie multiple solutions
    #!!! We are biased towards choosing the bottom right.
    def get_max_pos_and_score(self):         
        best_row = 0
        best_col = 0
        best_score = 0
        for row in range(self.num_row):
            for col in range(self.num_col):
                this_score = self.decision_matrix[row][col][0]
                if this_score >= best_score:
                    best_row = row
                    best_col = col
                    best_score = this_score
        return (best_row, best_col, best_score)
        #else:
        #    sys.stderr.write("WARNING: Attempt to get best score before alignment")

    def print_alignment(self):
        """
        Base on the completed decision matrix, obtain the
        best starting position and then trace the path
        towards the top left using the decisions records.
        The alignement representation is build from right to left
        """
        row, col, best_score = self.get_max_pos_and_score()
        current_decision = self.decision_matrix[row][col][1]
        num_insertion = 0
        num_deletion = 0
        ref_string = ''
        query_string = ''
        while current_decision != 'T':
            current_decision = self.decision_matrix[row][col][1]
            if current_decision == 'M':
                ref_string = self.reference[col-1] + ref_string
                query_string = self.query[row-1] + query_string
                row -= 1
                col -= 1
            elif current_decision == 'D':
                ref_string = '-'+ ref_string
                query_string = self.query[row-1] + query_string
                row -= 1
                num_deletion += 1
            elif current_decision == 'I':
                ref_string = self.reference[col-1] + ref_string
                query_string = '-' + query_string
                col -= 1
                num_insertion += 1
            elif current_decision == 'T':
                print('Alignment score: %s'%best_score)
                print('Reference position: %s'%str((col, col + len(ref_string))))
                print('Query position: %s'%str((row, row + len(query_string))))
                print(ref_string)
                print(query_string)



def main():
    """
    Main(). Parse user inputs, initialise
    alignment object, and print relevant results.
    """
    user_inputs = parse_args()
    alignment = SWalign(match=user_inputs.match,
                        mismatch=user_inputs.mismatch,
                        gap=user_inputs.gap,
                        gap_ext=user_inputs.gap_ext,
                        mismatch_exponent=user_inputs.mismatch_exponent,
                        gap_exponent=user_inputs.gap_exponent,
                        costume_score=user_inputs.costume_score,
                        verbose=user_inputs.verbose)
    print(alignment.align(user_inputs.query, user_inputs.ref))


if __name__ == '__main__':
    main()
