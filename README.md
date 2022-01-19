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
  ├── mpwt/ (see mpwt at https://github.com/AuReMe/mpwt)
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


## **main.py :**

You will need all the files described in "Folders structure" above and the correct folder structure. If not, 
PlantGEMs will ask you for every file one by one, which is definitely not convenient. 

_main.ini_ is mandatory, the process will exit as soon as it does not see it in the _main_directory_.

__Example of use :__
```bash
PlantGEMs/python/files/directory$ python main.py path/to/main/directory/
```
__Help displayed with the associated argument :__
```bash
PlantGEMs/python/files/directory$ python main.py -h
usage: main.py [-h] [-v] [-i [0-100]] [-d [0-100]] [-ev [0-1]] [-c [0-100]] [-bs [0-1000]] main_directory

positional arguments:
  main_directory        The path to the main directory where the 'files/' directory is stored

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Toggle the printing of more information
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
```
**NB** : each argument between [brackets] is optional.

## **blasting.py only :**
You will need : (see also _Folders structure_)
* The model's sbml file (.sbml)
* The model's proteome fasta (.faa)
* The subject's proteome fasta (.faa)

You can store the files in a _files/_ directory if you reconstruct several networks, or specify a path for each file needed if you do only one. If you fail to give it on the commandline you'll be prompted to give the different path.
<br />

**NB** : Please note that the "name" of the model is its _**id**_. You will also need a _main.ini_ file stored in the _main_directory_ you provided, with at least the organisms names. Every additional data will not be read (you can use the same file as for _mpwting_ if you manually used this process before the _blasting_ process)

__Example of use for multiple reconstruction with _main.ini_:__
```bash
PlantGEMs/python/files/directory$ python blasting.py path/to/main/directory/
```

__Example of use for single reconstruction without _main.ini_:__
```bash
PlantGEMs/python/files/directory$ python blasting.py unique_pipeline name path/to/main/directory/
```
**NB** : If you put the files correctly in the "files" folder (ie : every file with the appropriate extension and the 
species' name as filename), PlantGEMs will find and use them so that you don't have to specify them. If not, you'll be 
prompted to give an exact path to those files. You can then use a single reconstruction instruction on a previously used
folder with many species by calling only one of them. Or you can use the optional arguments and specify the exact path 
of each needed file (see below).

__Help displayed with the associated argument :__
```bash
PlantGEMs/python/files/directory$ python blasting.py -h
usage: blasting.py [-h] [-v] [-u] [-rr RERUN] [-n NAME] [-m MODEL_FILE_PATH] [-mfaa MODEL_PROTEOMIC_FASTA_PATH] [-sfaa SUBJECT_PROTEOMIC_FASTA_PATH] [-sgff SUBJECT_GFF_PATH]
                   [-i [0-100]] [-d [0-100]] [-ev [0-1]] [-c [0-100]] [-bs [0-1000]]
                   main_directory

positional arguments:
  main_directory        The path to the main directory where the 'files/' directory is stored

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Toggle the printing of more information
  -u, --unique          Specify if the reconstruction is made on a unique species or not
  -rr RERUN, --rerun RERUN
                        Use this option if you want to rerun the blast selection on an existing blasted.pkl object and give its path
  -n NAME, --name NAME  The future draft's name
  -m MODEL_FILE_PATH, --model_file_path MODEL_FILE_PATH
                        Model's file's path, use if 'files/' directory doesn't exist
  -mfaa MODEL_PROTEOMIC_FASTA_PATH, --model_proteomic_fasta_path MODEL_PROTEOMIC_FASTA_PATH
                        Model's proteomic fasta's path, use if 'files/' directory doesn't exist
  -sfaa SUBJECT_PROTEOMIC_FASTA_PATH, --subject_proteomic_fasta_path SUBJECT_PROTEOMIC_FASTA_PATH
                        Subject's proteomic fasta's path, use if 'files/' directory doesn't exist
  -sgff SUBJECT_GFF_PATH, --subject_gff_path SUBJECT_GFF_PATH
                        Subject's gff file's path, use if 'files/' directory doesn't exist
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
```
**NB** : each argument between [brackets] is optional.

## **mpwting.py only :**
This module needs you to create a "files/" directory in a directory of your choice (can be the same as for _blasting.py_
above for single use, **has** to be the same for a complete pipeline use with _main.py_) and put in there the following 
files, for each organism you want to reconstruct :
* all the .fna (e.g.: _grape.fna_, _kiwi.fna_...)
* all the .gff (e.g.: _grape.gff_, _kiwi.gff_...)
* all the .tsv from EggNOG-mapper annotation (you have to do that manually here : http://eggnog-mapper.embl.de/
on the proteomic fasta (.faa) file (e.g.: _grape.tsv_).

**NB** : GFF parsing won't work if the CDS' lines have the protein version number associated to their names (protein.1.**1**, protein.2.**2**...) Issue to be fixed.

The directory in which is stored the _files/_ directory is called the _main_directory_. It is the only argument needed by this module.
\
Then, make a _main.ini_ file corresponding to your needs and save it under the _main_directory_.

__Example of use :__
```bash
PlantGEMs/python/files/directory$ python mpwting.py path/to/main/directory
```

__Help displayed with the associated argument :__
```bash
PlantGEMs/python/files/directory$ python mpwting.py -h
usage: mpwting.py [-h] [-v] main_directory

positional arguments:
  main_directory  The path to the main directory where the 'files/' directory is stored

optional arguments:
  -h, --help      show this help message and exit
  -v, --verbose   Toggle the printing of more information
```
**NB** : each argument between [brackets] is optional.

**NB2** : if you didn't put the files in the _files/_ directory, you will be asked to give the exact path for each file needed. Not recommended if you reconstruct several organisms at once for obvious practicality.

**NB3** : every dependency needed is normally listed in the _requirements.txt_, but you will also need **Pathway-Tools** to be installed. Please see mpwt's GitHub page for more information : https://github.com/AuReMe/mpwt.


## **merging.py only :**
This module merges all the models (sbml/json) and the .dat flat files from Pathway Tools that are present in the merge folder.
This folder is either created automatically if the pipeline is launched with _main.py_, or you need to make it and fill it with the files you want.
See the structure (above) to make it. You can merge as many sbml and json models you want but only one of each flat files (reactions/enzrxns/proteins.dat) from Pathway Tools.

__Example of use :__
```bash
PlantGEMs/python/files/directory$ python merging.py path/to/main/directory
```

__Help displayed with the associated argument :__
```bash
usage: merging.py [-h] [-v] main_directory

positional arguments:
  main_directory  The path to the main directory where the 'files/' directory is stored

optional arguments:
  -h, --help      show this help message and exit
  -v, --verbose   Toggle the printing of more information
```
**NB** : each argument between [brackets] is optional.

## Files description :

### -- Python files --

- ``blasting.py`` -- Creation of plant draft from a template model using Blast.

- ``main.py`` -- Main file to launch all the workflow with a single command line.

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
- A complete pipeline (blast draft (done) & mpwt draft (done) + merging (done) + gap-filling + global analysis of the reconstructed networks).
  - The blast values ~~will be~~ **are** customizable as it suits you when using the rerun function.
  - _Merging.py_ reworking is ~~on its way~~ **done !**
    - A merge of any sbml **and json** model ~~will also be~~ **is** possible in the _merging.py_ module.
  - The compartmentalization of the reactions.
  - The Meneco-based gap-filling (https://github.com/bioasp/meneco) will be done subsequently.
  
Further ideas :
- Multiple models reconstruction in the blast module.
- Graphical interface for a common use and easy metabolic reconstructions.
- A docker version for more compatibility and easy deployment.
