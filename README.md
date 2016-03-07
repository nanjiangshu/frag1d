# Frag1D

## Description
Frag1D is a program for predicting one-dimensional (1D) structures of proteins
from their amino acid sequences. The program Frag1D predicts three types of 1D
structures, i.e. three-state secondary structure (H: Helix, S: Sheet, R: Random 
Coil), eight-state Shape Strings (R, S, U, V, A, K, G and T, see the definition 
in the paper) and three-state Shape Strings (R, S, U, V -> S; A, K -> H; G, T -> T). 
The program is witten in c/c++ and shell scripts.  Currently, Frag1D can be 
run on Linux and Windows (with Cygwin)

Frag1D is copyrighted (c) to Nanjiang Shu and Tuping Zhou, Structural
Chemistry, Stockholm University, Sweden and is free for academic use.

##Reference:
Tuping Zhou<sup>1</sup>, Nanjiang Shu<sup>1</sup> and Sven Hovmoller. A Novel Method for
Accurate One-dimensional Protein Structure Prediction based on Fragment
Matching, Bioinformatics, 2010;26(4):470-477. (<sup>1</sup>Co-first author)

##Contact:
Nanjiang Shu, Science for Life Laboratory, Stockholm
Stockholm University, Sweden
####    Email: 

nanjiang.shu@scilifelab.se

nanjiang.shu@mmk.su.se

##Web-server
The web server of Frag1D is available at 
http://frag1d.bioshu.se


##Installation:

Download the package by

    $ git clone https://github.com/nanjiangshu/frag1d

Enter the frag1d directory

    $ cd frag1d

Fetch the data set by 
    
    $ git lfs fetch
    $ git lfs checkout

If you don't have git-lfs installed, please installed follow the instructions
at https://git-lfs.github.com/

After that, run

    $ make 
    $ make install

Make sure that the NCBI nr database formatted for PSI-BLAST is installed. The
environmental variable BLASTDB points to the directory storing nr blast
database needed by PSI-BLAST

    $ export BLASTDB=path-storing-blast-nr-database


##Usage
```
usage: frag1d.sh  [-cpu INT] [-outpath DIR] [-blastdb STR] [-db STR]
                  [-pssm FILE] [-pssmfilelist LISTFILE]
                  [-showexample] [-verbose] [-version] [-not-clean] [-h]
                  [-outname STR]
                   FILE

Note: the supplied sequence should be in FASTA format

Options:
 -cpu           INT  Set the number of cores to be used to run blastpgp
                     (default: 1)
 -outpath       DIR  Output the result to the specified path, (default:./)
 -outname       STR  Output name, (default: query)
 -blastdb      FILE  Database for psi-blast, (default: nr)
 -db            STR  Database for Frag1D, (default: nr30)
 -pssm         FILE  Supply pssm file in PSI-BLAST -Q flag output format,
                     If supplied, PSI-BLAST will not run
 -pssmfilelist FILE  Supply a file with a list of pssm files for batch mode
                     prediction
 -not-clean          If supplied, the intermediate files will be kept
 -verbose            Print verbose information 
 -version            Print version 
 -showexample        Print How to use with examples
 -h|--help           Print this help message and exit

```
##Examples:
In the subfolder `test`

* Carry out the prediction by supplying a single sequence file in FASTA format

        $ ../frag1d.sh test1.aa

* Carry out the prediction by supplying a single pssm file

        $ ../frag1d.sh --pssm test2.pssm

* Carry out the prediction by a list file with a number of sequence files

        $ ../frag1d.sh --seqfilelist test3.seqfilelist


##Others:
The blastpgp used in this version of FRAG1D is 2.2.25.
If you want to use a different version of PSIBLAST, please copy the program
`blastpgp` to $FRAG1D/bin and the corresponding matrix file to $FRAG1D/data

    $ cp blastpgp $FRAG1D/bin
    $ cp BLOSUM62 $FRAG1D/data

##Running time:

It takes on average about 6 minutes to predict a protein sequence with 300
amino acids, when running on a single core with 2GHZ cpu and 1GB RAM. Most
time is taken by PSI-BLAST for building the sequence profile. When the
sequence profile is obtained, it takes about half a minute to get the
prediction result.


##Result example
An example output result file can be found in `test/testseq.predfrag1d`

```
#One-dimensional structure prediction by Frag1D (2009), for entry testseq
# Explanation:
# Num      Residue index in the sequence.
# AA       One-letter amino acid code.
# Sec      Predicted three state secondary structure,
#          (H, S and R).
# ConfSec  Confidence of the predicted secondary structure.
# S8       Predicted 8 state Shape String,
#          (R, S, U, V, A, K, G and T).
# ConfS8   Confidence of the predicted 8 state Shape String.
# S3       Predicted 3 state Shape String
#          (R, S, U, V -> S; A, K -> H; G, T - > T). 
# ConfS3   Confidence of the predicted 3 state Shape String.
# Num AA Sec ConfSec  S8 ConfS8  S3 ConfS3
    1  E   R  0.874   A  0.667   H  0.726
    2  E   R  0.979   K  0.429   H  0.521
    3  G   R  0.943   T  0.727   T  0.778
    4  L   R  0.821   R  0.565   S  0.788
    5  D   R  0.807   R  0.552   S  0.834
    6  F   R  0.850   R  0.543   S  0.914
    7  P   R  0.797   R  0.488   S  0.698
    8  E   R  0.774   S  0.478   S  0.620
    ...
```


