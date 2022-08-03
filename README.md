# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database and an already curated model.

![Alt text](./Flowchart_PlantGEMs.drawio.svg)

# Contents :
 - [Folders structure](#folders-structure-common-to-all-modules-and-main-) 
 - [How to use](#how-to-use)
   - [main.py](#mainpy-)
   - [blasting.py](#blastingpy-only-)
   - [mpwting.py](#mpwtingpy-only-)
   - [merging.py](#mergingpy-only-)
 - [Files description](#files-description-)
 - [NEWS](#news-)
 
NB : Each script can be launched separately. However, the pipeline can be entirely launched with _main.py_.<br />
NB' : Every module can be used in CLI or in a Python script. Don't forget to use a virtual environment meeting the requirements.txt file specifications.

## Folders' structure to follow (common to all modules and _main_) :

```text
main_directory/
  ├── files/
  │    ├── metacyc.json (merging.py)
  │    ├── model.sbml (blasting.py)
  │    ├── model.faa (blasting.py)
  │    ├── species_1.gff (blasting.py + mpwting.py)
  │    ├── species_1.faa (blasting.py + EggNOG Mapper)
  │    ├── species_1.fna (mpwting.py)
  │    ├── species_1.tsv (=EggNOG annotation on species.faa file) (mpwting.py)
  │    ├── species_2.gff
  │    ├── species_2.faa
  │    ├── species_2.fna
  │    ├── species_2.tsv
  │    └── ...
  └── main.ini
```

## Output of the entire pipeline :
```text
main_directory/
  ├── files/
  │    ├── metacyc.json (merging.py)
  │    ├── metacyc_ids.tsv (created by merging.py)
  │    ├── model.sbml (blasting.py)
  │    ├── model.faa (blasting.py)
  │    ├── species_1.gff (blasting.py + mpwting.py)
  │    ├── species_1.faa (blasting.py + EggNOG Mapper)
  │    ├── species_1.fna (mpwting.py)
  │    ├── species_1.tsv (=EggNOG annotation on species.faa file) (mpwting.py)
  │    ├── species_2.gff
  │    ├── species_2.faa
  │    ├── species_2.fna
  │    ├── species_2.tsv
  │    └── ...
  ├── blast/
  │    ├── species_1/
  │    │    ├── species_1_blast_draft.json
  │    │    ├── protein_gene_correspondence.tsv
  │    │    └── objects_history/
  │    │         ├── blasted.pkl
  │    │         ├── drafted.pkl
  │    │         └── genes_selected.pkl
  │    ├── species_2/
  │    └── ...
  ├── mpwt/ 
  │    └── (see mpwt at https://github.com/AuReMe/mpwt)
  ├── merge/
  │    ├── species_1/
  │    │    ├── species_1_blast_draft.json
  │    │    ├── species_1_merged.json (Currently PlantGEMs' final result)
  │    │    ├── enzrxns.dat
  │    │    ├── proteins.dat
  │    │    └── reactions.dat
  │    ├── species_2/
  │    └── ...
  └── main.ini
```

## HOW TO USE


## **PlantGEMs.py :**

### - The whole pipeline :

You will need all the files described in "Folders structure" above and the correct folders structure. If not, 
PlantGEMs will ask you for every file one by one, which is definitely not convenient. 

_main.ini_ is mandatory for each and every use, the process will exit as soon as it does not see it in the _main_directory_. It needs to be filled with file with the name, genetic element and taxon information for each species.

__Example of use :__
```bash
~PlantGEMs/PlantGEMs/src/$ python PlantGEMs.py run path/to/main/directory/
```
__Help displayed with the associated argument :__
```txt
(venv) [asa@asa-precision PlantGEMs]$ python PlantGEMs/src/PlantGEMs.py -h
usage: PlantGEMs.py [-h] [-le] [-rr RERUN] [-i [0-100]] [-d [0-100]] [-ev [0-1]] [-c [0-100]] [-bs [0-1000]] [-m] {run,blasting,mpwting,merging} main_directory

positional arguments:
  {run,blasting,mpwting,merging}
                        Choice of the module you'd like to use.
                        run = whole pipeline
                        B or blasting = only blasting module
                        P or mpwting = only mpwt from AuReMe
                        M or merging = only merging module
  main_directory        The path to the main directory where the 'files/' directory is stored

options:
  -h, --help            show this help message and exit
  -le, --log_erase      Erase the existing log file to create a brand new one
  -rr RERUN, --rerun RERUN
                        Use this option if you want to rerun the blast selection on an existing blasted.pkl object. The species' name is expected here
  -i [0-100], --identity [0-100]
                        The blast's identity percentage tolerated. Default=50
  -d [0-100], --difference [0-100]
                        The tolerated length difference between the two aligned sequences. Default=30
  -ev [0-1], --e_val [0-1]
                        The blast's e-value threshold value. Default=e-100
  -c [0-100], --coverage [0-100]
                        The minimum sequence coverage tolerated. Default=20
  -bs [0-1000], --bit_score [0-1000]
                        The blast's bit-score threshold value. Default=300
  -m, --migrate         Take the files previously created by blasting or mpwting modules and store them in a new merge folder, ready to be merged.

```

### - blasting only :
You will need : (see also _Folders structure_)
* The model's sbml file (.sbml)
* The model's proteome fasta (.faa)
* The subject's proteome fasta (.faa)

You can store the files in a _files/_ directory or specify a path for each file needed, except for the main.ini which has to be in the **main_directory**.
<br />

**NB** : Please note that the "name" of the model is its _**id**_ in the SBML file.

__Example of use :__
```bash
~PlantGEMs/PlantGEMs/src/$ python PlantGEMs.py blasting path/to/main/directory/
```

### - mpwting only :
This module needs you to create a "files/" directory in a directory of your choice (can be the same as for _blasting.py_
above) and put in there the following files, for each organism you want to reconstruct :
* all the .fna (e.g.: _grape.fna_, _kiwi.fna_...)
* all the .gff (e.g.: _grape.gff_, _kiwi.gff_...)
* all the .tsv from EggNOG-mapper annotation (you have to do that manually here : http://eggnog-mapper.embl.de/
on the proteomic fasta (.faa) file) (e.g.: _grape.tsv_).

__Example of use :__
```bash
~PlantGEMs/PlantGEMs/src/$ python PlantGEMs.py mpwting path/to/main/directory
```

**NB** : every dependency needed is normally listed in the _requirements.txt_, but you will also need **Pathway-Tools** to be installed. Please see mpwt's GitHub page for more information : https://github.com/AuReMe/mpwt.


### - merging only :
This module merges all the models (sbml/json) and the .dat flat files from Pathway Tools that are present in the merge folder.
This folder is either created automatically if the pipeline is launched with _main.py_, or you need to make it and fill it with the files you want.
See the structure (above) to make it. You can merge as many sbml and json models you want but only one of each flat files (reactions/enzrxns/proteins.dat) from Pathway Tools.

__Example of use :__
```bash
~PlantGEMs/PlantGEMs/src/$ python PlantGEMs.py merging path/to/main/directory
```

## Files description :

### -- Python files --

- ``blasting.py`` -- Creation of plant draft from a template model using Blast.

- ``graphing.py`` -- Creation of UpsetPlots (to be continued...).

- ``PlantGEMs.py`` -- Main file to launch all the workflow with a single command line.

- ``merging.py`` -- Merge metabolic networks, one from the Metacyc database using Pathway Tools (cf. mpwting.py) and the others from homemade reconstructions (cf. blasting.py) or from already curated models or other draft software.

- ``module.py`` -- File for the parent class of all the modules, contains useful methods that can be inherited in all module's classes.

- ``mpwting.py`` -- Preparation of files to run Pathway Tools (http://bioinformatics.ai.sri.com/ptools/) automatically and create a draft based on Metacyc with the mpwt library (see https://github.com/AuReMe/mpwt).

- ``utils.py`` -- Utility file to avoid code redundancy.

[//]: # (- ``menecoing.py`` -- Performs a gap filling of the model with Meneco &#40;Work In Progress&#41;.)
[//]: # (- ``graph.py`` -- Utility file to create different graphs and statistical analysis on the networks.)

<br />

___main.ini_ structure (example with 3 fruits) :__
```text
[1]                     # order doesn't matter as long as it begins with 1 and that all subsequents numbers are present (e.g.: 2,1,3,5,4...).
ORGANISM_NAME = kiwi    # name of every file related to the species/organism (e.g.: kiwi.faa, kiwi.gff...).
ELEMENT_TYPE = :CONTIG  # structure of the .fna file (depending on the organism's DNA assembly, can be NONE, :CHRSM or :CONTIG).
NCBI_TAXON_ID = 3625    # well... the NCBI taxon id ! Find it on Taxonomy DB for the exact species you work on.
[2]
ORGANISM_NAME = grape
ELEMENT_TYPE = :CHRSM
NCBI_TAXON_ID = 29760
[3]
ORGANISM_NAME = cucumber
ELEMENT_TYPE = :CHRSM
NCBI_TAXON_ID = 3659
```
## NEWS :

Some improvements are to come : 
- A complete pipeline (blast draft (done) & mpwt draft (re-WIP) + merging (done) + compartmentalization (WIP) + gap-filling + global analysis of the reconstructed networks).
  - The rework of the MPWT module with easier use
  - The compartmentalization of the reactions.
  - The Meneco-based gap-filling (https://github.com/bioasp/meneco) will be done subsequently.
  
Further ideas :
- Multiple models reconstruction in the blast module.
- Integration of the graphical interface developed by 4 Master's students.
- A PyPI release.
- A docker version for more compatibility and easy deployment.
