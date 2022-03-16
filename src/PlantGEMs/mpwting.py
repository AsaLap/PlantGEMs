# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the preparation of the required files for the 
Pathway Tools software reconstruction and launching of reconstruction 
using mpwt package from AuReMe."""

import argparse
import logging
import mpwt
import multiprocessing
import utils


def make_taxon_file(directory, taxon_name_list):
    """Function to make the taxon_id.tsv file.

    PARAMS:
        taxon_name_list (list) -- the list containing the name of the organism and its taxon id.
    """

    res = [["species", "taxon_id", "element_type"]]
    for i in taxon_name_list:
        res.append([i[0], str(i[1]), str(i[2])])
    utils.write_csv(directory, "taxon_id", res, separator="\t")


def make_mpwt_architecture(main_directory):
    """Prepare the folders and structure needed to run mpwt.

    ARGS:
        main_directory (str) - the folder in which create the structure.
    """

    if utils.check_path(main_directory):
        parameters = utils.read_config(main_directory + "main.ini")
        mpwt_directory = main_directory + "mpwt/"
        input_directory = mpwt_directory + "input/"
        output_directory = mpwt_directory + "output/"
        log_directory = mpwt_directory + "log/"
        utils.make_directory(mpwt_directory)
        utils.make_directory(input_directory)
        utils.make_directory(output_directory)
        utils.make_directory(log_directory)
        taxon_name_list = []
        multirun_list = []
        cpu = len(parameters.keys()) - 1
        count = 0
        for i in parameters.keys():
            if i != "DEFAULT":
                count += 1
                species_name = parameters[i]["ORGANISM_NAME"]
                element_type = parameters[i]["ELEMENT_TYPE"]
                taxon_id = int(parameters[i]["NCBI_TAXON_ID"])
                species_directory = input_directory + species_name + "/"
                utils.make_directory(species_directory)
                taxon_name_list.append([species_name, taxon_id, element_type])
                multirun_list.append([main_directory + "files/", input_directory, 4])
        make_taxon_file(input_directory, taxon_name_list)
        logging.info("Species found : {}".format(count))
        return [cpu, input_directory, output_directory, log_directory, multirun_list]


def make_mpwt_input_files(cpu, multirun_list):
    for i in multirun_list:
        mpwt.to_pathologic.create_pathologic_file(i[0], i[1], i[2])


def run_mpwt(cpu, input_directory, output_directory, log_directory):
    """
    Split of major function 'run', second part = launching the process on each given object with multiprocessing
    and launching the mpwt reconstruction.
    """

    mpwt.multiprocess_pwt(input_folder=input_directory, output_folder=output_directory, patho_inference=True,
                          patho_hole_filler=False, patho_operon_predictor=False, pathway_score=1, flat_creation=True,
                          dat_extraction=True, number_cpu=cpu, size_reduction=False, patho_log=log_directory,
                          taxon_file=input_directory + "taxon_id.tsv", verbose=True)


def run(args):
    cpu, input_directory, output_directory, log_directory, multirun_list = \
        make_mpwt_architecture(utils.slash(args.main_directory))
    make_mpwt_input_files(cpu, multirun_list)
    run_mpwt(cpu, input_directory, output_directory, log_directory)


def mpwt_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_directory", help="The path to the main directory where the 'files/' directory is stored",
                        type=str)
    parser.add_argument("-v", "--verbose", help="Toggle the printing of more information", action="store_true")
    parser.add_argument("-le", "--log_erase", help="Erase the existing log file to create a brand new one",
                        action="store_true")
    args = parser.parse_args()
    return args


def main():
    args = mpwt_arguments()
    if args.log_erase:
        logging.basicConfig(filename=args.main_directory + '/mpwting.log', filemode='w', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename=args.main_directory + '/mpwting.log', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    if args.verbose:
        logging.getLogger().addHandler(logging.StreamHandler())
    logging.info("------ Mpwting module started ------")
    run(args)


if __name__ == "__main__":
    main()
