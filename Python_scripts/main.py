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
    blasting.pipeline(main_directory)  # Launching the Blast
    mpwting.pipeline(main_directory)  # Then, launching MPWT


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
