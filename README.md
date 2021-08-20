# PlantGEMs : Project that aims to automatically reconstruct several plant metabolic networks, using the Metacyc database.


## Files description :

### -- Python files --

Python/blasting.py : Creation of plant draft from a template model.

Python/blastingTK.py : GUI of blasting.py (aborted).

Python/fusion.py : Merge two metabolic networks, one from a Metacyc reconstruction with Pathway Tools (cf. pathwayToolsPrep.py) and the other one from an homemade script which is reconstructed based on an already curated model (cf. blasting.py).

Python/gap_filling.py : Performs a gap filling of the model with Meneco. (Work In Progress).

Python/graph.py : Utilitary file to create different graphs and statistical analysis.

Python/pathwayToolsPrep.py : Preparation of files to run Pathway Tools automatically and create a draft based on Metacyc.

Python/utils.py : Utilitary file to avoid code repetition.

### -- R files --

R/graph.R : Make different graphs for the purpose of a study on the Blast scores to use.

R/regression.R : Regression graph for the same purpose as above.

R/tresholdSearch.R : Research of a treshold value for the score chosen above.

### -- Info files --

example_Blasting.ini : Example of an ini file for the blasting pipeline.

example_PT.ini : Example of an ini file for the pathwayToolsPrep pipeline.

index.txt : Example of index file for the pathwayToolsPrep pipeline.


## HOW TO USE :

For the moment, each script output is part of the next script's input (except 1 and 2 which are independent) and you'll need several biologic/bioinformatics files, which will be soon describe. For the moment, please look at the "pipeline" function in each script to see what is needed, along with the indication in the different .ini files.
Each script has a pipeline-like function wich makes all the things work (normally...).

Use in this order (except 1 and 2 that you can reverse) :
1. blasting.py
2. pathwayToolsPrep.py
3. fusion.py
4. gap_filling.py

Further improvement are to come, such has a complete pipeline and a requirements.txt.
