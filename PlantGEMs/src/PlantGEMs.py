# coding: utf8
# python 3.8
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux Aquitaine - DRC
# PhD : Genome scale metabolic reconstruction of several fruits
# 11/2021 - 10/2024

"""This file is the main, calling the others files to create an entirely new metabolic network."""

import argparse
import logging
from argparse import RawTextHelpFormatter

import blasting
import merging
import mpwting
import utils


def run(args):
    print("Proceeding to create all the needed files and checking input files, please stay around...")
    # Launching the first part of Blast (files checking & folder generation)
    list_objects_to_blast = blasting.blast_multirun_first(args.main_directory, args.identity, args.difference,
                                                          args.e_val, args.coverage, args.bit_score)
    # Launching the first part of MPWT (files checking & folder generation)
    list_objects, cpu, input_directory, output_directory, log_directory = \
        mpwting.mpwt_multirun_first(args.main_directory)

    # Then, launching the rest of the run in multiprocess without the need of any input from the user
    print("Everything's fine, now launching BLAST and then MPWT processes, it may take some time...")
    blasting.blast_multirun_last(list_objects_to_blast)
    mpwting.mpwt_multirun_last(list_objects, cpu, input_directory, output_directory, log_directory)

    # Merging all the new drafts and pathway tools pgdbs for each organism
    utils.migrate(args.main_directory)
    merging.run(args.main_directory)


def main_arguments():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    # POSITIONAL
    parser.add_argument("module",
                        help="Choice of the module you'd like to use.\nrun = whole pipeline"
                             "\nB or blasting = only blasting module"
                             "\nP or mpwting = only mpwt from AuReMe"
                             "\nM or merging = only merging module",
                        type=str,
                        choices=["run", "blasting", "mpwting", "merging"])
    parser.add_argument("main_directory",
                        help="The path to the main directory where the 'files/' directory is stored",
                        type=str)

    # GENERAL OPTIONS
    parser.add_argument("-le", "--log_erase",
                        help="Erase the existing log file to create a brand new one",
                        action="store_true")
    # TODO : implement verbose option
    # parser.add_argument("-v", "--verbose",
    #                     help="Toggle the printing of more information",
    #                     action="store_true")

    # BLASTING OPTIONS
    parser.add_argument("-rr", "--rerun",
                        help="Use this option if you want to rerun the blast selection on an existing blasted.pkl "
                             "object. The species' name is expected here",
                        type=str)
    parser.add_argument("-i", "--identity",
                        help="The blast's identity percentage tolerated. Default=50",
                        type=int,
                        default=50,
                        choices=range(0, 101),
                        metavar="[0-100]")
    parser.add_argument("-d", "--difference",
                        help="The tolerated length difference between the two aligned sequences. Default=30",
                        type=int,
                        default=30,
                        choices=range(0, 101),
                        metavar="[0-100]")
    parser.add_argument("-ev", "--e_val",
                        help="The blast's e-value threshold value. Default=e-100",
                        type=float,
                        default=1e-100,
                        choices=range(0, 1),
                        metavar="[0-1]")
    parser.add_argument("-c", "--coverage",
                        help="The minimum sequence coverage tolerated. Default=20",
                        type=int,
                        default=20,
                        choices=range(0, 101),
                        metavar="[0-100]")
    parser.add_argument("-bs", "--bit_score",
                        help="The blast's bit-score threshold value. Default=300",
                        type=int,
                        default=300,
                        choices=range(0, 1001),
                        metavar="[0-1000]")

    # MPWTING OPTIONS
    # TODO : add the mpwt options

    # MERGING OPTIONS
    parser.add_argument("-m", "--migrate",
                        help="Take the files previously created by blasting or mpwting modules and store them in a "
                             "new merge folder, ready to be merged.",
                        action="store_true")
    args = parser.parse_args()
    return args


def main():
    args = main_arguments()
    if args.log_erase:
        logging.basicConfig(filename=args.main_directory + '/PlantGEMs.log', filemode='w', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename=args.main_directory + '/PlantGEMs.log', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    if args.module in ["B", "BLASTING", "Blasting", "blasting"]:  # Args for the blasting launch
        logging.info("------ Blasting module started ------")
        if args.rerun:
            blasting.rerun_blast_selection(args.main_directory, args.rerun, args.identity, args.difference, args.e_val,
                                           args.coverage, args.bit_score)
        else:
            blasting.run(args.main_directory, args.identity, args.difference, args.e_val,
                         args.coverage, args.bit_score)
    if args.module in ["P", "MPWTING", "Mpwting", "mpwting"]:  # Args for the mpwting launch
        logging.info("------ Mpwting module started ------")
        mpwting.run(args.main_directory)
    if args.module in ["M", "MERGING", "Merging", "merging"]:  # Args for the merging launch
        logging.info("------ Merging module started ------")
        if args.migrate:
            utils.migrate(utils.slash(args.main_directory))
        merging.run(args.main_directory)
    else:  # Else, run of the whole pipeline
        logging.info("------ PlantGEMs pipeline started ------")
        run(args)


if __name__ == "__main__":
    main()
