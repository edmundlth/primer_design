Hi-Plex Primer Design Tool
==========================

Program Description:
--------------------
Hi-Plex is a highly multiplexed PCR system. This software is the primer design software
for the system. 
The program can be installed using pip (Python's package manager) as follow:
    $pip install https://github.com/edmundlth/primer_design

The program is a command-line tool which takes in the following arguments,
(use command-line --help for more description):

* --bed *Specifying the path to the BED-file which contains the coordinates*  
* --fa *Specify the path to the directory containing the fasta files of reference sequences*  
* --log *Specify the logging file name*
* --outfile *Specify the name of the file recording primer information*
* --auxfile *Specify the name of the file recording auxiliary data*
* --sense\_heel *Nucleic acid string specifying the sense heel sequence*
* --antisense\_heel *Nucleic acid string specifying the antisense heel sequence*  
* --tiles *2 number specifying the minimum and maximum tile size allowed (inclusive)*  
* --primer\_length *A number specifying the optimal primer length*  
* --primer\_length\_var *A number specifying the allowed variation in primer length*
* --allowed\_overlap *A number specifying the maximum allowed overlaping between tiles*
* --score\_func *Specify the scoring function to be used*
* --tm\_func *Specify the Melting temperature prediction function to be used*
* --target\_tm *Specify the ideal melting temperature*  
* --tm\_underachieve *Specify the penalty weight for primer with tm < target\_tm*
* --saltcorr *specify the salt correction  method to be used*
* --conc *specify the concentration of various ions*
* --gc\_weight *Specify the importance of a G-C binding relative to A-T*


Note: The program require that the BED-file naming convention to be consistent with
      the fasta files names containing the reference sequences, i.e.,  
      the name of each fasta file should be the same of the names provided at the  
      the first column of the BED-file

Using the input parameters, for each region coordinate specifed in the BED-file,  
the program make a exhausive search of all possible primer set satisfying the  
constrains that completely tile the region.


