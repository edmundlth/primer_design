"""
This module aim to build a tool that interpret bwa output
and use it for specificity (or even dimer) checking.

The idea is to incorporate the information generated into
the initialisation of Score class. 
"""


import pysam
import os
import sys
import logging
import subprocess

from utils import rev_complement, complement

DEFAULT_PRIMER_FASTA = './all_primers.fa'

class Specificity(object):
    def __init__(self, reference_file, all_primers):
        """
        reference_file: a fasta file that is assummed to
        be bwa-indexed.
        all_primers: a dictionary with primers sequences as keys
        """
        self.primer_fasta = DEFAULT_PRIMER_FASTA
        self.file_dir = os.path.dirname(self.primer_fasta)
        self.file_name = os.path.basename(self.primer_fasta).split('.')[0]

        self.primers = all_primers
        self.tags_dict = {}
        generate_fasta(all_primers, self.primer_fasta)
        generate_alignment_files(reference_file, self.file_dir, self.file_name)
        self.make_specificity_lookup()

    def make_specificity_lookup(self):
        bam_file = os.path.join(self.file_dir, self.file_name + '.bam')
        alignments = pysam.AlignmentFile(bam_file)
        count = 0
        for align_entry in alignments.fetch():
            seq = align_entry.seq
            # Flag #16 implies the sequence is mapped to the reverse strand
            if align_entry.flag == 16:
                seq = rev_complement(seq)
            tags = {"X0":None, "X1":None, "XA":None,
                    "XM":None, "XO":None, "XG":None,
                    "AS":None, "NM":None}
            for tag in tags:
                if align_entry.has_tag(tag):
                    tags[tag] = align_entry.get_tag(tag)
            self.primers[seq] = [specificity_score(tags)]




def specificity_score(alignment_data):
    #!! Specificity score is negated number of hits
    #!! should also include information about the match quality 
    #!! and the binding energy 
    if alignment_data["X0"] != None:
        return -1 * alignment_data["X0"]
    else:
        sys.stderr.write("Warning: Num best hit not available\n")
        raise RuntimeError


#!! a bignumber is supplied to bwa samse -n option so that
#!! it retains all "alternative hits" of the sequence.
BIGNUMBER = 10 
def generate_alignment_files(reference_file, file_dir, file_name):
    """
    Generate required BAM and BAM-Index file needed for
    pysam to have random access to alignment data.
    The procedure:
        1) Use "bwa aln reference_file primer_fasta" to align 
           all the primers in the "primer_fasta" file to 
           the reference sequence from "reference_file".
           This produces the .sai file (suffix array index)
        2) Use "bwa samse reference_file <..>.sai primer_fasta" 
           to generate the corresponding .sam file
        3) use "samtools view -b <..>.sam" to convert the .sam
           to .bam
        4) use "samtools sort <..>.bam" to sort the .bam file
           so that we can use "samtools index sorted_<..>.bam"
           to generate .bam.bai indexed bamfile for pysam
    """
    fa_file = os.path.join(file_dir, file_name + '.fa')
    sai_file = os.path.join(file_dir, file_name + '.sai')
    sam_file = os.path.join(file_dir, file_name + '.sam')
    bam_file = os.path.join(file_dir, file_name + '.bam')
    bam_index = os.path.join(file_dir, file_name)

    with open(sai_file, 'w') as sai:
        aln_process = subprocess.call(['bwa', 'aln', 
                                        reference_file, 
                                        fa_file],
                                        stdout=sai)
    with open(sam_file, 'w') as sam:
        subprocess.call(['bwa', 'samse', '-n %i'%BIGNUMBER, 
                         reference_file, sai_file, fa_file],
                         stdout=sam)
    with open(bam_file, 'w') as bam:
        subprocess.call(['samtools', 'view', '-b', sam_file], stdout=bam)

    subprocess.call(['samtools', 'sort', bam_file, bam_index])
    subprocess.call(['samtools', 'index', bam_file])




def generate_fasta(sequences, fasta_name=DEFAULT_PRIMER_FASTA):
    with open(fasta_name, 'w') as outfile:
        for seq in sequences:
            outfile.write('>\n' + seq + '\n')


def robust_command_call(command_string):
    command = command_string.split(' ')
    exit_code = subprocess.call(command)
    if exit_code != 0:
        logging.info("Exit Code: %i"%exit_code)














def _test():
    reference_file = "../primer_design/fasta/mock_ref.fa"
    all_primers = {'ggggTTTT'.upper():{},
                   'aaaaaaaa'.upper():{}}
    test_specificity = Specificity(reference_file, all_primers)
    print(test_specificity.primers)


if __name__ == '__main__':
    _test()
