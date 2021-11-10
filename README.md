# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database and an already curated model.

![Alt text](./Flowchart_PlantGEMs.drawio.svg)

## Files description :

### -- Python files --

Python_scripts/blasting.py : Creation of plant draft from a template model using Blast.

<!-- Python_scripts/menecoing.py : Performs a gap filling of the model with Meneco. (Work In Progress). -->
<!--  -->
Python_scripts/graph.py : Utilitary file to create different graphs and statistical analysis.

Python_scripts/main.py : Main file to launch all the workflow with a simple command (WIP).

Python_scripts/merging.py : Merge two draft metabolic networks, one from the Metacyc database using Pathway Tools (cf. mpwting.py) and the other one from a homemade reconstruction based on an already curated model (cf. blasting.py).

Python_scripts/module.py : File for the parent class of all the modules, contains useful methods that can be heritated in all module's class.

Python_scripts/mpwting.py : Preparation of files to run Pathway Tools automatically and create a draft based on Metacyc with the mpwt library (see https://github.com/AuReMe/mpwt).

Python_scripts/utils.py : Utilitary file to avoid code redundancy.

### -- Info files --

example.ini : Example of a file you'd give as an input to launch the pipeline.


## HOW TO USE :

Each script can be launched separately. However, the pipeline is soon to be connected (see _main.py_ WIP) with one simple initialisation file as input like the "example.ini" you can see on this git.

Both **blasting.py** and **mpwting.py** are currently working and can be launched in a bash console or used as python functions.

### Description for CLI usage.

**For blasting.py**, you will need :
* The model's sbml file
* The model's fasta (protein)
* The subject's fasta

Just like for mpwting (see below) you can store the files in a _files/_ directory and name them this way : _organism_name.extension_.
Please note that the "name" of the model is its _**id**_.


For a normal use with proper file sorting :
```bash
$ python blasting.py pipeline "name_of_reconstructed_organism" "path/to/main/directory/"
```

**NB** : if you didn't put the files in the _files/_ directory, you will be asked to give the exact path for each file needed. You can also choose to specify them following the structure below (every parameter is a string) :
```bash
$ python blasting.py pipeline [name][main_directory][optional=model_file_path][optional=model_fasta_path][optional=subject_fasta_path]
```

<br />

**For mpwting.py :** create a "files/" directory in a directory of your choice (can be the same as for _blasting.py_ above) and put in there, for each organism you want to reconstruct the following files, all named like this : _organism_name.extension_ :
  * all the .fasta (ex. : _tomato.fasta_, _kiwi.fasta_...)
  * all the .gff (ex. : _tomato.gff_)
  * all the .tsv from EggNOG-mapper annotation (you have to do that manually here : http://eggnog-mapper.embl.de/) (ex. : _tomato.tsv_)  

The directory in which is stored the _files/_ directory is called the _main_directory_. It is the only argument needed by this module.  

Then, make a _main.ini_ file corresponding to your needs and save it under the _main_directory_.
_Further instructions on how to build this file and what the different fields mean are to come. (WIP)_

__Example of use :__
```bash
$ python mpwting.py pipeline "path/to/main/directory"
```

**NB** : if you didn't put the files in the _files/_ directory, you will be asked to give the exact path for each file needed. Not recommended if you reconstruct several organisms at once for obvious practicality.

**NB2** : every dependency needed is normally listed in the _requirements.txt_, but you will also need **Pathway-Tools** to be installed. Please see mpwt's GitHub page for more information : https://github.com/AuReMe/mpwt.

<br />

Further improvements are to come, such as a complete pipeline and an automatic gathering of the mandatory parameters in the main.ini file.
_Merging_ reworking is coming and soon wired to the previous modules. Same thing for a Meneco-based gap-filling (https://github.com/bioasp/meneco).
