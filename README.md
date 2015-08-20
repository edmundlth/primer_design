Hi-Plex Primer Design Tool
==========================

Downloading
-----------
The program can be installed using pip (Python's package manager) as follow:  

    $pip install https://github.com/edmundlth/primer_design  


Example usage
-------------
A test case is included in the Python package downloaded.  
The BED-file for the test case (example.bed) is located in the same directory
as primer\_design.py - the main module - together with a directory named 'fasta'
containing a file named 'example.fa' containing the reference sequence needed.  

Use the following command in the directory containing primer\_design.py:  

$python primer\_design.py --log 'example.log' --bed 'example.bed' --fa './example\_fasta' --outfile 'example.tsv' --sense\_heel ctctctatgggcagtcggtgatt --antisense\_heel ctgcgtgtctccgactcag --tiles 50 55 --primer\_length 10 --primer\_length\_var 8 --allowed\_overlap 5 --score\_func 'score\_Lp' --tm\_func 'Tm\_NN' --target\_tm 50.0 --tm\_underachieve 1.0 --saltcorr 5 --conc 50 0 0 0 0 25 25 --gc\_weight 1.5 --auxfile 'example\_aux.tsv'

Program Description
--------------------
Hi-Plex is a highly multiplexed PCR system. This software is the primer design software
for the system.  
  
The program is a command-line tool which takes in the following arguments,  
(use command-line --help for more description):  

* --bed     **Specifying the path to the BED-file which contains the coordinates**  
* --fa     **Specify the path to the directory containing the fasta files of reference sequences**  
* --log     **Specify the logging file name**
* --outfile     **Specify the name of the file recording primer information**
* --auxfile     **Specify the name of the file recording auxiliary data**
* --sense\_heel     **Nucleic acid string specifying the sense heel sequence**
* --antisense\_heel     **Nucleic acid string specifying the antisense heel sequence**  
* --tiles     **2 number specifying the minimum and maximum tile size allowed (inclusive)**  
* --primer\_length     **A number specifying the optimal primer length**  
* --primer\_length\_var     **A number specifying the allowed variation in primer length**
* --allowed\_overlap     **A number specifying the maximum allowed overlaping between tiles**
* --score\_func     **Specify the scoring function to be used**
* --tm\_func     **Specify the Melting temperature prediction function to be used**
* --target\_tm     **Specify the ideal melting temperature**  
* --tm\_underachieve     **Specify the penalty weight for primer with tm < target\_tm**
* --saltcorr     **Specify the salt correction  method to be used**
* --conc     **Specify the concentration of various ions**
* --gc\_weight     **Specify the importance of a G-C binding relative to A-T**

  
Note:  
The program require that the BED-file naming convention to be consistent with
the fasta files names containing the reference sequences, i.e.,  
the name of each fasta file should be the same of the names provided at the  
the first column of the BED-file

Using the input parameters, for each region coordinate specifed in the BED-file,  
the program make a exhausive search of all possible primer set satisfying the  
constrains that completely tile the region.


