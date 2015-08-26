Hi-Plex Primer Design Tool
==========================

Downloading
-----------
We suggest using python virtual environment 'virtualenv' (ref: https://virtualenv.pypa.io/en/latest/) to install the program by following the steps below:
1. If you haven't install virtualenv, use pip (Python's package manager) to install it:  
  `pip install -U virtualenv`  
2. Change to a working directory of your choice (eg. 'hiplex_dir')  
  `mkdir hiplex_dir`  
  `cd hiplex_dir`  
3. Create the virtual environment with the name of your choice (eg. 'hiplex_env'):  
  `virtualenv hiplex_env`  
4. Activate the virtual environment:  
  `source hiplex_env/bin/activate`  
5. Clone the git repository containing the package:  
  `git clone https://github.com/edmundlth/primer_design`  
6. Now we are ready to install the package into this virtual environment by:  
  `pip install -U primer_design`  
7. Use `primer_design -h` to see command-line arguments available and their descriptions.  

Example usage
-------------
A directory called 'example' is included in the packaged.
As an example usage of the program:  
1. (say we are currently in hiplex_dir from above)  
  `cd primer_design/example`  
You will find the following files here:  
  'example.bed' -- contain 2 BED-file formated coordinates  
  'example.fa' -- contain a sequence of DNA in fasta format  
2. run the program in this directory using the following command-line argument:  
  `primer_design --log example.log --bed example.bed --fa . --outfile example.tsv --sense_heel ctctctatgggcagtcggtgatt --antisense_heel ctgcgtgtctccgactcag --tiles 50 55 --primer_length 10 --primer_length_var 8 --allowed_overlap 5 --score_func score_Lp --tm_func Tm_NN --target_tm 50.0 --tm_underachieve 1.0 --NN_table 3 --saltcorr 5 --conc 50 0 0 0 0 25 25 --gc_weight 1.5 --auxfile example_aux.tsv`  

Program Description
--------------------
Hi-Plex is a highly multiplexed PCR system. This software is the primer design software
for the system. Using the input parameters, for each region coordinate specifed in the BED-file,
the program make a exhausive search of all possible primer set satisfying the
constrains that completely tile the region.


Note:
The program require that the BED-file naming convention to be consistent with
the fasta files names containing the reference sequences, i.e.,
the name of each fasta file should be the same of the names provided at the
the first column of the BED-file
