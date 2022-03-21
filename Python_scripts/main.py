# coding: utf8
# python 3.8
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux Aquitaine - DRC
# PhD : Genome scale metabolic reconstruction of several fruits
# 11/2021 - 10/2024

"""This file is the main, calling the others files to create an entirely new metabolic network."""

import argparse
import blasting
import merging
import mpwting
import utils


def run(args):
    print("Proceeding to create all the needed files and checking input files, please stay around...")
    # Launching the first part of Blast (files checking & folder generation)
    list_objects_to_blast = blasting.blast_multirun_first(args)
    # Launching the first part of MPWT (files checking & folder generation)
    list_objects, cpu, input_directory, output_directory, log_directory = \
        mpwting.mpwt_multirun_first(args.main_directory)

    # TODO : adding a check-file function here if needed (PDP)

    # Then, launching the rest of the run in multiprocess without the need of any input from the user
    print("Everything's fine, now launching BLAST and then MPWT processes, it may take some time...")
    blasting.blast_multirun_last(list_objects_to_blast)
    mpwting.mpwt_multirun_last(list_objects, cpu, input_directory, output_directory, log_directory)

    # Merging all the new drafts and pathway tools pgdbs for each organism
    utils.migrate(args.main_directory)
    merging.run(args.main_directory)


def main_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_directory", help="The path to the main directory where the 'files/' directory is stored",
                        type=str)
    parser.add_argument("-v", "--verbose", help="Toggle the printing of more information", action="store_true")
    parser.add_argument("-i", "--identity", help="The blast's identity percentage tolerated. Default=50",
                        type=int, default=50, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-d", "--difference",
                        help="The tolerated length difference between the two aligned sequences. Default=30",
                        type=int, default=30, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-ev", "--e_val",
                        help="The blast's e-value threshold value. Default=e-100",
                        type=float, default=1e-100, choices=range(0, 1), metavar="[0-1]")
    parser.add_argument("-c", "--coverage", help="The minimum sequence coverage tolerated. Default=20",
                        type=int, default=20, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-bs", "--bit_score", help="The blast's bit-score threshold value. Default=300",
                        type=int, default=300, choices=range(0, 1001), metavar="[0-1000]")
    args = parser.parse_args()
    return args


def main():
    args = main_arguments()
    run(args)


if __name__ == "__main__":
    main()
