# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database and an already curated model.

![Alt text](./Flowchart_PlantGEMs.drawio.svg)

## Files description :

### -- Python files --

Python_scripts/blasting.py : Creation of plant draft from a template model using Blast.

Python_scripts/merging.py : Merge two draft metabolic networks, one from the Metacyc database using Pathway Tools (cf. mpwting.py) and the other one from a homemade reconstruction based on an already curated model (cf. blasting.py).

<!-- Python_scripts/gap_filling.py : Performs a gap filling of the model with Meneco. (Work In Progress). -->
<!--  -->
Python_scripts/graph.py : Utilitary file to create different graphs and statistical analysis.

Python_scripts/mpwting.py : Preparation of files to run Pathway Tools automatically and create a draft based on Metacyc with the mpwt library (see https://github.com/AuReMe/mpwt).

Python_scripts/main.py : Main file to launch all the workflow with a simple command (WIP).

Python_scripts/utils.py : Utilitary file to avoid code repetition.

### -- Info files --

example_Blasting.ini : Example of an ini file for the blasting pipeline (soon obsolete).

example_PT.ini : Example of an ini file for the pathwayToolsPrep pipeline (soon obsolete).

index.txt : Example of index file for the mpwting pipeline (soon obsolete).


## HOW TO USE :

For the moment, each script's output is part of the next script's input (except 1 and 2 which are independent) and you'll need several biologic/bioinformatics files. For the moment, please look at the "pipeline" function at the end of each script to see what is needed, along with the indication in the different .ini files and the worflow's scheme above.

Use in this order (except 1 and 2 that you can reverse) :
1. blasting.py
2. mpwting.py
3. merging.py
4. menecoing.py (next to come, not on GitHub yet, WIP)

Further improvement are to come, such has a complete pipeline (see "main.py") and a gap filling option.
