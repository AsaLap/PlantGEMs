# coding: utf8
# python 3.8
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux Aquitaine - DRC
# PhD : Genome scale metabolic reconstruction of several fruits
# 11/2021 - 10/2024

"""This file is the main, calling the others files to create an entirely new metabolic network."""

import blasting
import merging
import mpwting
import utils
import sys


def migrate(main_directory):
    blast_directory = main_directory + "blast/"
    merge_directory = main_directory + "merge/"
    mpwt_directory = main_directory + "mpwt/"
    utils.make_directory(merge_directory)
    list_species = utils.get_list_directory(blast_directory)
    for species in list_species:
        utils.make_directory(merge_directory + species)
        utils.copy_file(blast_directory + species + "/" + species + "_blast_draft.json",
                        merge_directory + species + "/" + species + "_blast_draft.json")
        utils.copy_file(mpwt_directory + "/output/" + species + "/reactions.dat",
                        merge_directory + species + "/reactions.dat")
        utils.copy_file(mpwt_directory + "/output/" + species + "/proteins.dat",
                        merge_directory + species + "/proteins.dat")
        utils.copy_file(mpwt_directory + "/output/" + species + "/enzrxns.dat",
                        merge_directory + species + "/enzrxns.dat")


def run(main_directory):
    print("Proceeding to create all the needed files and checking input files, please stay around...")
    # Launching the first part of Blast (files checking & folder generation)
    list_objects_to_blast = blasting.blast_multirun_first(main_directory)
    # Launching the first part of MPWT (files checking & folder generation)
    list_objects, cpu, input_directory, output_directory, log_directory = mpwting.mpwt_multirun_first(main_directory)

    # TODO : adding a check-file function here if needed

    # Then, launching the rest of the run in multiprocess without the need of any input from the user
    print("Everything's fine, now launching BLAST and then MPWT processes, it may take some time...")
    blasting.blast_multirun_last(list_objects_to_blast)
    mpwting.mpwt_multirun_last(list_objects, cpu, input_directory, output_directory, log_directory)

    # Merging all the new drafts and pathway tools pgdbs for each organism
    migrate(main_directory)
    merging.run(main_directory)


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
