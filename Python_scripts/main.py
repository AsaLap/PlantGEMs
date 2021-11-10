# coding: utf8
# python 3.8
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux Aquitaine - DRC
# PhD : Genome scale metabolic reconstruction of several fruits
# 11/2021 - 10/2024

"""This file is the main, calling the others files to create an entirely new metabolic network."""

import blasting
import mpwting
import os
import utils
import sys


def run(main_directory):
    main_parameters = utils.read_config("main.ini" + main_directory)
    species_names = []
    for i in main_parameters.keys():
        if i != "DEFAULT":
            species_names.append(main_parameters[i]["ORGANISM_NAME"])
    
    if os.path.isdir(main_directory):
        for name in species_names:  # Launching blast for every species
            blasting.pipeline(name, main_directory)
        mpwting.pipeline(main_directory)  # Then, launching MPWT


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
