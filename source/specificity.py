"""
This module aim to build a tool that interpret bwa output
and use it for specificity (or even dimer) checking.

The idea is to incorporate the information generated into
the initialisation of Score class. 

The main function in the module is "make_specificity_lookup"
which takes a reference sequence fasta file 
(single or multi sequences) an iterable object containing all
the sequences to be aligned to the reference and return
a dictionary {sequence:[specificity score]}.
"""


import pysam
import os
import sys
import shutil
import logging
import subprocess
import tempfile

from utils import rev_complement
import error

#!! a bignumber is supplied to bwa samse -n option so that
#!! it retains all "alternative hits" of the sequence.
DEFAULT_NUM_THREADS = 2
BIGNUMBER = 10
def make_specificity_lookup(reference_file, all_primers):
    specificity_lookup = {}
    alignments = get_sam_alignments(reference_file, all_primers)
    for align_entry in alignments:
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
        assert seq in all_primers
        specificity_lookup[seq] = [specificity_score(tags)]
    return specificity_lookup

def specificity_score(alignment_data):
    #!! Specificity score is negated number of hits
    #!! should also include information about the match quality 
    #!! and the binding energy 
    if alignment_data["X0"] != None:
        num_optimal_hits = alignment_data["X0"]
    else:
        sys.stderr.write("Warning: Num best hit not available\n")
        raise RuntimeError
    if alignment_data["X1"] != None:
        num_suboptimal_hits = alignment_data["X1"]
    else:
        num_suboptimal_hits = 0
    return -1 * (num_optimal_hits + num_suboptimal_hits)


#############################################################

def get_sam_alignments(reference_fa_name, all_primers):
    """
    Using commandline tools to generate the necessary 
    files for pysam to read and fetch all the 
    alignment informations.
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

    All the files are generated in a temporary directory which 
    will be deleted once all the alignments had been obtained.
    """
    tempdir = tempfile.mkdtemp()
    try:
        primer_fa_name = os.path.join(tempdir, 'primers.fa')
        sai_filename = os.path.join(tempdir, 'samindex.sai')
        sam_filename = os.path.join(tempdir, 'samfile.sam')
        bam_filename = os.path.join(tempdir, 'bamfile.bam')
        bam_index_prefix = os.path.join(tempdir, 'bamfile')

        generate_fasta(all_primers, primer_fa_name) 
        with open(sai_filename, 'w') as sai_file:
            #!! Should ask user for number of threads??
            subprocess.call(['bwa', 'aln', 
                             '-t %i'%DEFAULT_NUM_THREADS,
                             reference_fa_name, 
                             primer_fa_name],
                             stdout=sai_file)
        with open(sam_filename, 'w') as sam_file:
            subprocess.call(['bwa', 'samse', '-n %i'%BIGNUMBER, 
                             reference_fa_name, 
                             sai_filename, 
                             primer_fa_name],
                             stdout=sam_file)
        with open(bam_filename, 'w') as bam_file:
            subprocess.call(['samtools', 'view', 
                             '-b', sam_filename], 
                            stdout=bam_file)

        subprocess.call(['samtools', 'sort', 
                         bam_filename, 
                         bam_index_prefix])
        subprocess.call(['samtools', 'index', 
                         bam_filename])
        alignments = [align for align in 
                      pysam.AlignmentFile(bam_filename).fetch()]
    finally:
        shutil.rmtree(tempdir)
        assert not os.path.exists(tempdir)
    return alignments


def generate_fasta(sequences, fasta_name):
    with open(fasta_name, 'w') as outfile:
        for seq in sequences:
            outfile.write('>\n' + seq + '\n')


def _error_and_exit(exit_code, process_name):
    processes = {"aln":"BWA ALN FAILED",
                 "samse":"BWA SAMSE FAIL"}
    if process_name == "aln":
        if exit_code != 0:
            sys.exit(error.ERRORS["BWA ALN FAILED"])
    elif process_name == "samse":
        sys.exit(error.ERRORS["BWA SAMSE FAILED"])
    elif process_name == "samview":
        pass
    pass



