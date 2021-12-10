# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database and an already curated model.

![Alt text](./Flowchart_PlantGEMs.drawio.svg)

## HOW TO USE :

Each script can be launched separately. However, the pipeline is soon to be connected (see _main.py_ WIP) with one simple initialisation file as input like the "example.ini" you can see on this git.

Both **blasting.py** and **mpwting.py** are currently working and can be launched in a bash console or used as python functions.

## Description for CLI usage :

### **For blasting.py** :
You will need :
* The model's sbml file
* The model's proteome fasta (.faa)
* The subject's proteome fasta (.faa)
* The subject's genome fasta (.fna)

Just like for mpwting (see below) you can store the files in a _files/_ directory and name them this way : _organism_name.extension_. 
\
Please note that the "name" of the model is its _**id**_. You will also need a _main.ini_ file stored in the _main_directory_ you provided, with at least the organisms names. Every additional data will not be read (you can use the same file as for _mpwting_ if you manually used this process before the _blasting_ process)

__Example of use :__
```bash
$ python blasting.py pipeline "path/to/main/directory/"
```

[comment]: <> (For a normal use with proper file sorting :)

[comment]: <> (```bash)

[comment]: <> ($ python blasting.py pipeline "name_of_reconstructed_organism" "path/to/main/directory/")

[comment]: <> (```)

[comment]: <> (**NB** : if you didn't put the files in the _files/_ directory, you will be asked to give the exact path for each file needed. You can also choose to specify them following the structure below &#40;every parameter is a string&#41; :)

[comment]: <> (```bash)

[comment]: <> ($ python blasting.py pipeline [name][main_directory][optional=model_file_path][optional=model_fasta_path][optional=subject_fasta_path])

[comment]: <> (```)

<br />

### **For mpwting.py :**
Create a "files/" directory in a directory of your choice (can be the same as for _blasting.py_ above) and put in there, for each organism you want to reconstruct the following files, all named like this : _organism_name.extension_ :
* all the .fasta (ex. : _vitis.fna_, _vitis.faa_, _kiwi.fna_, _kiwi.faa_...)
* all the .gff (ex. : _vitis.gff_)
* all the .tsv from EggNOG-mapper annotation (you have to do that manually here : http://eggnog-mapper.embl.de/) (ex. : _vitis.tsv_) on the proteomic fasta (faa) file.

**NB** : GFF parsing won't work if the CDS have different names than the protein they are linked to in the .faa file.

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

- ``mpwting.py`` -- Preparation of files to run Pathway Tools automatically and create a draft based on Metacyc with the mpwt library (see https://github.com/AuReMe/mpwt).

- ``utils.py`` -- Utility file to avoid code redundancy.

[//]: # (- ``menecoing.py`` -- Performs a gap filling of the model with Meneco &#40;Work In Progress&#41;.)
[//]: # (- ``graph.py`` -- Utility file to create different graphs and statistical analysis on the networks.)

### -- Info files --

- ``example.ini`` -- Example of a file you'd give as an input to launch the pipeline.

<br />

## NEWS :

Some improvements are to come : 
- A complete pipeline (blast draft & mpwt draft + merging + gap-filling + global analysis of the reconstructed networks).
  - The blast values will be customizable as it suits you.
  - _Merging.py_ reworking is on its way (next module to come).
    - A merge of any sbml model will also be possible in the _merging.py_ module.
  - The compartmentalization of the reactions.
  - The Meneco-based gap-filling (https://github.com/bioasp/meneco) will be subsequently.
  
Further ideas :
- Multiple models reconstruction in the blast module.
- Graphical interface for a common user utility and easy metabolic reconstruction.
