<<<<<<< HEAD
# Plant-GEMs : Project that aims to the automatic reconstruction of several plants metabolic networks.
=======
# Plant-GEMs : Project that aims to the automatic reconstruction of several plant metabolic networks.
>>>>>>> bbb8783c6c529e5eef0f5e5fc5f75463b93b7672

-- Python files --

blasting.py : Creation of plant draft model from a template model.

blastingTK.py : GUI of blasting.py (work in progress).

fusion.py : Merge two metabolic networks, one from a Metacyc reconstruction with Pathway Tools (cf. pathwayToolsPrep.py) and the other one from an homemade script which is reconstructed based on an already curated model (cf. blasting.py).

gap_filling.py : Performs a gap filling of the model with Meneco. (Work In Progress).

graph.py : Utilitary file to create different graphs and statistical analysis.

pathwayToolsPrep.py : Preparation of files to run Pathway Tools automatically and create a draft based on Metacyc.

utils.py : Utilitary file to avoid code repetition.

-- R files --

R/graph.R : Make different graphs for the purpose of a study on the Balst scores to use.

R/regression.R : Regression graph for the same purpose as above.

R/tresholdSearch.R : Research of a treshold value for the score chosen above.

-- Other files --

example_Blasting.ini : Example of an ini file for the blasting pipeline.

example_PT.ini : Example of an ini file for the pathwayToolsPrep pipeline.

index.txt : Example of index file for the pathwayToolsPrep pipeline.
