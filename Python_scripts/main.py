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
    print("Proceding to create all the needed files and checking input files, please stay around...")
    # Launching the first part of Blast (files checking & generation)
    list_objects_to_blast = blasting.sub_pipeline_first(main_directory)
    # Launching the first part of MPWT (files checking & generation)
    list_objects, cpu, input_directory, output_directory, log_directory = mpwting.sub_pipeline_first(main_directory)

    # Then, launching the rest of the pipeline in multiprocess without the need of any input from the user
    print("Everything's fine, now launching BLAST and MPWT processes, it may take some time...")
    blasting.sub_pipeline_last(list_objects_to_blast)
    mpwting.sub_pipeline_last(list_objects, cpu, input_directory, output_directory, log_directory)


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
