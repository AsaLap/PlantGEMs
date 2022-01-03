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
    print("Proceeding to create all the needed files and checking input files, please stay around...")
    # Launching the first part of Blast (files checking & folder generation)
    list_objects_to_blast = blasting.multirun_first(main_directory)
    # Launching the first part of MPWT (files checking & folder generation)
    list_objects, cpu, input_directory, output_directory, log_directory = mpwting.multirun_first(main_directory)

    # TODO : adding a check-file function here if needed

    # Then, launching the rest of the run in multiprocess without the need of any input from the user
    print("Everything's fine, now launching BLAST and then MPWT processes, it may take some time...")
    blasting.multirun_last(list_objects_to_blast)
    mpwting.multirun_last(list_objects, cpu, input_directory, output_directory, log_directory)

    # TODO
    # Migrate each organism's files to merge directory
    # migrate()  # Declare this function in utils.py to be used in merging.py by the user if wished so

    # launch merge
    # merge()  # uses each subdirectory (species) of the /merging directory


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
