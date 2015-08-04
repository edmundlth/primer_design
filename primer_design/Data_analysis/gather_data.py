
from score import (
                   score_primer_Lp,
                   score_primer_Lp_3,
                   score_primer_linear,
                   score_primer_linear_3,
                   score_primer_sum,
                   score_primer_sum_3,
                  )
from utils import (
                    raw_tm,
                    raw_hairpin,
                    raw_entropy,
                    raw_dimer,
                    raw_self_dimer,
                    raw_gc,
                    raw_gc_clamp,
                    raw_run,
                    k_entropy,
                  )
from score_utils import (
                         entropy_score,
                         hairpin_score,
                         dimer_score,
                         self_dimer_score,
                         run_score,
                         gc_score,
                         gc_clamp_score,
                         tm_score,
                        )
FERGUS_OUT_HEADER =  ['dimer', 'entropy_3prime_half', 
                      'entropy_3prime_10mer', 'k_entropy', 
                      'L2_3_score', 'L2_score', 'L3_3_score', 
                      'linear_3_score', 'sum_3_score','sum_score']
def gather_all_fergus():
    infile = open('fergus_primer_L2.tsv')
    header = next(infile).rstrip().split('\t')
    primer_set = []
    for line in infile:
        line = line.rstrip().split('\t')
        primer = line[1]
        primer_set.append(primer)
    max_primer_length = len(max(primer_set, key = len))
    infile = open('fergus_primer_L2.tsv')
    next(infile)
    outfile = open('fergus_primers_data.tsv','w')
    out_header = FERGUS_OUT_HEADER
    outfile.write('\t'.join(out_header) + '\n')
    for line in infile:
        line = line.rstrip().split('\t')
        sequence = line[0]
        dimer = raw_all_to_all_dimer(primer, primer_set)


        L2_3_score = score_primer_Lp_3(sequence,TARGET_TM, p = 2)
        L2_score = score_primer_Lp(sequence,TARGET_TM, p =2)
        L3_3_score = score_primer_Lp_3(sequence, TARGET_TM,p = 3)
        linear_3_score = score_primer_linear_3(sequence,TARGET_TM)
        linear_score = score_primer_linear(sequence,TARGET_TM)
        sum_3_score = score_primer_sum(sequence,TARGET_TM)
        sum_score = score_primer_sum_3(sequence,TARGET_TM)


        #3' end effects using 3' = primer[l//2:]
        primer_3prime = sequence[len(sequence)//2:]
        entropy_3prime_half = raw_entropy(primer_3prime)

        #3' end effects using 10-mer
        primer_3primer = sequence[-10:]
        entropy_3prime_10mer = raw_entropy(primer_3prime)

        # all k_mer entropy score using k_entropy
        k_entrop = k_entropy(sequence)

        output = line +  map(str, [dimer, entropy_3prime_half, 
                                  entropy_3prime_10mer, k_entropy, 
                                  L2_3_score, L2_score, L3_3_score, 
                                  linear_3_score, sum_3_score, sum_score])


        outfile.write('\t'.join(output)+'\n')
        outfile.flush()
    outfile.close()
