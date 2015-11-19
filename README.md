#MetaDisorder

Disordered proteins are highly abundant in regulatory processes such as transcription and cell-signaling. Different methods
have been developed to predict protein disorder often focusing on different types of disordered regions.
Here, we present MD, a novel META-Disorder prediction method that molds various sources of information predominantly obtained 
from orthogonal prediction methods. MD significantly outperformed its constituents, and compared favorably to other top 
prediction methods. In sustained cross-validation, MD not only outperforms its origins, but it also compares favorably 
to other state-of-the-art prediction methods in a variety of tests that we applied. MD is capable of predicting disordered 
regions of all "flavors", and identifying new ones that are not captured by other predictors.

##How to install the package
=============================================
1)sudo apt-get install python-software-properties
2)sudo apt-add-repository "deb http://rostlab.org/debian/ stable main contrib non-free"
3)sudo apt-get update (ignore GPG error)
4)sudo apt-get install rostlab-debian-keyring (without verification)
5)sudo apt-get update
6)sudo apt-get install metadisorder 

SYNOPSIS
metadisorder [OPTION]

## Method Description

* Authors: Markus Schmidberger and Guy Yachdav
* Development year: 2010
* Languages: perl
* Publications: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0004433
* Description (ML ? )
* Training / Test Data
* ...


Output format
Self-annotating, see example outputs in /usr/share/metadisorder/example.

##OPTIONS

        chk 
         Path to psiblast checkpoint file for fasta input. Required.

       debug
         Debugging mode: [0|1]. Default: 0.

       fasta
         Path to fasta input

        -h, --help

       hssp 
         Path to hssp file for fasta input. Required.

       norsnet
         Path to norsnet prediction for fasta input

       prof
         Path to prof prediction for fasta input

       profbval_raw
         Path to profbval prediction for fasta input. Use mode 5 of profbval.

        profcon
          Path to profcon prediction for fasta input

        out 
         Path to output file

        out_mode
          Output format: [0|1] DESCRIBE FORMATS HERE

        --version
         Print version

        workdir
         Work directory, optional. If not defined a temporary directory is used.


##EXAMPLES
  
  Obtaining input files: see <https://rostlab.org/owiki/index.php/How_to_generate_an_HSSP_file_from_alignment>

  Example without profcon input:

     fasta=/usr/share/metadisorder/example/tmdfast.fasta hssp=/usr/share/metadisorder/example/tmdfast.hsspPsiFil prof=/usr/share/metadisorder/example/tmdfast.rdbProf profbval_raw=/usr/share/metadisorder/example/tmdfast.profbval
    norsnet=/usr/share/metadisorder/example/tmdfast.norsnet chk=/usr/share/metadisorder/example/tmdfast.chk out=tmdfast.noprofcon_mdisorder out_mode=1

  Example with profcon input (please note: profcon is really slow and is known not to improve predictions     significantly):

    metadisorder fasta=/usr/share/metadisorder/example/tmdfast.fasta hssp=/usr/share/metadisorder/example/tmdfast.hsspPsiFil prof=/usr/share/metadisorder/example/tmdfast.rdbProf profbval_raw=/usr/share/metadisorder/example/tmdfast.profbval norsnet=/usr/share/metadisorder/example/tmdfast.norsnet chk=/usr/share/metadisorder/example/tmdfast.chk profcon=/usr/share/metadisorder/example/tmdfast.profcon out=tmdfast.profcon_mdisorder out_mode=1

##ENVIRONMENT
  METADISORDERCONF
    Location of metadisorderrc configuration file to use overriding other configuration files

  FILES
  /usr/share/metadisorder/metadisorderrc.default
    Default configuration file. See this file for a description of the parameters.

  /etc/metadisorderrc
    System configuration file overriding values in /usr/share/metadisorder/metadisorderrc.default

  ~/.metadisorderrc
    User configuration file overriding values in /etc/metadisorderrc

   $METADISORDERCONF
     If this environment variable is set ~/.metadisorderrc is disregarded and the value of the variable is read        for configuration options overriding /etc/metadisorderrc

##RESTRICTIONS
Right now all input files must be given on the command line as you see in the examples. Autmatical generation of input files is not supported at present.  Let us know if you need this feature.

##REFERENCES
Schlessinger, A., Punta, M., Yachdav, G., Kajan, L., and Rost, B.
(2009). Improved disorder prediction by combination of orthogonal
approaches. PLoS ONE, 4(2), e4433.




