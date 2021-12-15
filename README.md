# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database and an already curated model.

![Alt text](./Flowchart_PlantGEMs.drawio.svg)

# Contents :
 - [Folders structure](#folders-structure-common-to-all-modules-and-main-) 
 - [How to use](#how-to-use)
   - [main.py](#mainpy-)
   - [blasting.py](#blastingpy-only-)
   - [mpwting.py](#mpwtingpy-only-)
 - [Files description](#files-description-)
 - [NEWS](#news-)

Each script can be launched separately. However, the pipeline can be entirely launched with _main.py_.

Every module can be used in CLI or in a Python script.

## Folders structure (common to all modules and main) :

```text
main_directory
 ├── files
 │    └── model.sbml (blasting.py)
 │    └── model.faa (blasting.py)
 │    └── species_1.gff (blasting.py + mpwting.py)
 │    └── species_1.faa (blasting.py + EggNOG Mapper)
 │    └── species_1.fna (mpwting.py)
 │    └── species_1.tsv (=EggNOG annotation on species.faa file) (mpwting.py)
 │    └── species_2.gff
 │    └── species_2.faa
 │    └── species_2.fna
 │    └── species_2.tsv
 │    └── ...
 └── main.ini
```

## HOW TO USE


### **main.py :**

You will need all the files described in "Folders structure" above and the correct folder structure. If not, 
PlantGEMs will ask you for every file one by one, which is definitely not convenient. 

_main.ini_ is mandatory, the process will exit as soon as it does not see it in the _main_directory_.

__Example of use :__
```bash
$ python main.py run "path/to/main/directory/"
```


### **blasting.py only :**
You will need : (see also _Folders structure_)
* The model's sbml file (.sbml)
* The model's proteome fasta (.faa)
* The subject's proteome fasta (.faa)

You can store the files in a _files/_ directory if you reconstruct several networks, or specify a path for each file needed if you do only one. If you fail to give it on the commandline you'll be prompted to give the different path.
<br />

**NB** : Please note that the "name" of the model is its _**id**_. You will also need a _main.ini_ file stored in the _main_directory_ you provided, with at least the organisms names. Every additional data will not be read (you can use the same file as for _mpwting_ if you manually used this process before the _blasting_ process)

__Example of use for multiple reconstruction with _main.ini_:__
```bash
$ python blasting.py pipeline "path/to/main/directory/"
```

__Example of use for single reconstruction without _main.ini_:__
```bash
$ python blasting.py unique_pipeline "name" "path/to/main/directory/" [model_file_path] [model_proteomic_fasta_path] [subject_proteomic_fasta_path] [subject_gff_path]
```
**NB** : the arguments in [brackets] are optional for the launch. If you put the files correctly in the "files" folder, 
PlantGEMs will find and use them, if not, you'll be prompted to give an exact path to those files. You can then use a 
single reconstruction instruction on a previously used folder with many species by calling only one of them.


### **mpwting.py only :**
This module needs you to create a "files/" directory in a directory of your choice (can be the same as for _blasting.py_
above for single use, **has** to be the same for a complete pipeline use with _main.py_) and put in there the following 
files, for each organism you want to reconstruct :
* all the .fna (e.g.: _grape.fna_, _kiwi.fna_...)
* all the .gff (e.g.: _grape.gff_, _kiwi.gff_...)
* all the .tsv from EggNOG-mapper annotation (you have to do that manually here : http://eggnog-mapper.embl.de/
on the proteomic fasta (.faa) file (e.g.: _grape.tsv_).

**NB** : GFF parsing won't work if the CDS' lines have the protein version number associated to their names (protein.1.**1**, protein.2.**2**...).

The directory in which is stored the _files/_ directory is called the _main_directory_. It is the only argument needed by this module.
\
Then, make a _main.ini_ file corresponding to your needs and save it under the _main_directory_.

__Example of use :__
```bash
$ python mpwting.py pipeline "path/to/main/directory"
```

**NB** : if you didn't put the files in the _files/_ directory, you will be asked to give the exact path for each file needed. Not recommended if you reconstruct several organisms at once for obvious practicality.

**NB2** : every dependency needed is normally listed in the _requirements.txt_, but you will also need **Pathway-Tools** to be installed. Please see mpwt's GitHub page for more information : https://github.com/AuReMe/mpwt.

## Files description :

### -- Python files --

- ``blasting.py`` -- Creation of plant draft from a template model using Blast.

- ``main.py`` -- Main file to launch all the workflow with a single command line.

- ``merging.py`` -- Merge two draft of metabolic networks, one from the Metacyc database using Pathway Tools (cf. mpwting.py) and the other one from a homemade reconstruction based on an already curated model (cf. blasting.py).

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
- A complete pipeline (blast draft & mpwt draft + merging + gap-filling + global analysis of the reconstructed networks).
  - The blast values will be customizable as it suits you.
  - _Merging.py_ reworking is on its way (next module to come).
    - A merge of any sbml model will also be possible in the _merging.py_ module.
  - The compartmentalization of the reactions.
  - The Meneco-based gap-filling (https://github.com/bioasp/meneco) will be done subsequently.
  
Further ideas :
- Multiple models reconstruction in the blast module.
- Graphical interface for a common user utility and easy metabolic reconstruction.
