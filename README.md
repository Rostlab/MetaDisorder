#MetaDisorder
MD, a novel META-Disorder prediction method that molds various sources of information predominantly obtained from orthogonal prediction methods which focus on different types of disordered regions.

##How to Install the Package

```shell
sudo apt-get install python-software-properties
sudo apt-add-repository "deb http://rostlab.org/debian/ stable main contrib non-free"
sudo apt-get update # ignore GPG error
sudo apt-get install rostlab-debian-keyring # without verification
sudo apt-get update
sudo apt-get install metadisorder
```

## How To Run, Basics

* **Usage:** metadisorder [OPTION]
    * In case of getting the error: "Can't locate Config/IniFiles.pm in @INC (you may need to install the Config::IniFiles module)", you can resolve it by executing the following command: ```sudo cpan install Config::IniFiles```

* **Obtaining Input Files:**  https://rostlab.org/owiki/index.php/How_to_generate_an_HSSP_file_from_alignment

* **Example of How to Run without Profcon Input:**

   ```
metadisorder fasta=/usr/share/metadisorder/example/tmdfast.fasta hssp=/usr/share/metadisorder/example/tmdfast.hsspPsiFil prof=/usr/share/metadisorder/example/tmdfast.rdbProf profbval_raw=/usr/share/metadisorder/example/tmdfast.profbval
    norsnet=/usr/share/metadisorder/example/tmdfast.norsnet chk=/usr/share/metadisorder/example/tmdfast.chk out=tmdfast.noprofcon_mdisorder out_mode=1
```
* **Example of How to Run with Profcon Input:** (please note: profcon is really slow and is known not to improve predictions significantly):

  ```
metadisorder fasta=/usr/share/metadisorder/example/tmdfast.fasta hssp=/usr/share/metadisorder/example/tmdfast.hsspPsiFil prof=/usr/share/metadisorder/example/tmdfast.rdbProf profbval_raw=/usr/share/metadisorder/example/tmdfast.profbval norsnet=/usr/share/metadisorder/example/tmdfast.norsnet chk=/usr/share/metadisorder/example/tmdfast.chk profcon=/usr/share/metadisorder/example/tmdfast.profcon out=tmdfast.profcon_mdisorder out_mode=1
  ```

* **Output** is specified by out='outputFileName'

* **Expected Results:** The output file is self-annotating. It contains a table with the following information:

    * Number - *residue number*
    * Residue - *amino-acid type*
    * NORSnet - *raw score by NORSnet (prediction of unstructured loops)*
    * NORS2st - *two-state prediction by NORSnet; D=disordered*
    * PROFbval - *raw score by PROFbval (prediction of residue flexibility from sequence)*
    * Bval2st - *two-state prediction by PROFbval*
    * Ucon - *raw score by Ucon (prediction of protein disorder using predicted internal contacts)*
    * Ucon2st - *two-state prediction by Ucon*
    * MD - *raw score by MD (prediction of protein disorder using orthogonal sources)*
    * MD_rel - *reliability of the prediction by MD; values range from 0-9. 9=strong prediction*
    * MD2st - *two-state prediction by MD*

You can also see example outputs in **/usr/share/metadisorder/example**.

## Method Description

* **Authors:** Avner Schlessinger, Marco Punta, Guy Yachdav, Laszlo Kajan, and Burkhard Rost
* **Publications:** http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0004433
* **Year:** 2009
* **Languages:** Perl
* **Description:** Disordered proteins are highly abundant in regulatory processes such as transcription and cell-signaling. Different methods have been developed to predict protein disorder often focusing on different types of disordered regions. Here, we present MD, a novel META-Disorder prediction method that molds various sources of information predominantly obtained  from orthogonal prediction methods. MD significantly outperformed its constituents, and compared favorably to other top prediction methods. In sustained cross-validation, MD not only outperforms its origins, but it also compares favorably to other state-of-the-art prediction methods in a variety of tests that we applied. MD is capable of predicting disordered regions of all "flavors", and identifying new ones that are not captured by other predictors.

##Link to Elixir
https://bio.tools/tool/mytum.de/MetaDisorder/1

##How To Run Extended

###Options

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

##Installation from Repository
[Warning] the source code includes many dependency and local directory, some of which might require installation/alteration. For the convenience, clone it to the "home directory".
- Blast, individual compilation of Psipred and the executable of Disopred are required for the package. 
```shell
cd 
git clone https://github.com/Rostlab/MetaDisorder
cd MetaDisorder
perl runMD.pl fasta=tmdfast.fasta
```
###Environment

* **METADISORDERCONF** - location of metadisorderrc configuration file to use overriding other configuration files
* **FILES**
    * **/usr/share/metadisorder/metadisorderrc.default** - default configuration file. See this file for a description of the parameters.
    * **/etc/metadisorderrc** -system configuration file overriding values in **/usr/share/metadisorder/metadisorderrc.default**
    * **~/.metadisorderrc** - user configuration file overriding values in **/etc/metadisorderrc**
    * **$METADISORDERCONF** - if this environment variable is set **~/.metadisorderrc** is disregarded and the value of the variable is read for configuration options overriding **/etc/metadisorderrc**.

##Training
MetaDisorder uses DisProt 3.4, which contains:
460 IDPs (intrinsically disordered proteins) 
1103 disordered regions, encompassing 35 functional categories (all based on published experimental data)

60 proteins with >780 residues were discarded as these could not be handled by all of the methods tested. From the remaining set, 17 more proteins crashed when applying at least one of the predictors in this study, and were also discarded.

##Evaluation
![Tag](https://cloud.githubusercontent.com/assets/13695363/13203314/aa2cafe0-d8b5-11e5-9670-cabc664b4289.png)
![Tag](https://cloud.githubusercontent.com/assets/13695363/13203315/aa2de428-d8b5-11e5-9f6f-f02e205cb078.png)

##Restrictions
Right now all input files must be given on the command line as you see in the examples. Automatical generation of input files is not supported at present.

##References
Schlessinger, A., Punta, M., Yachdav, G., Kajan, L., and Rost, B.
(2009). Improved disorder prediction by combination of orthogonal
approaches. PLoS ONE, 4(2), e4433.
