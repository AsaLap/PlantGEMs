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
import module
import mpwt
import multiprocessing
import numpy as np
import re
import utils


class Mpwting(module.Module):

    def __init__(self, _name, _main_directory, _element_type):
        super().__init__(_name, _main_directory)
        self.element_type = _element_type
        self.genomic_fasta_file_path = self._find_genomic_fasta(self.name)
        self.gff_file_path = self._find_gff(self.name)
        self.eggnog_file_path = self._find_eggnog(self.name)
        self.regions_dict = utils.get_sequence_region(self.gff_file_path)

    @property
    def directory(self):
        return self.main_directory + "mpwt/input/" + self.name + "/"

    def _make_dat_files(self):
        circular = 'N'
        dat_file_str_list = []
        if self.element_type == "NONE":
            for i in self.regions_dict.keys():
                dat_file_str_list.append('ID\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                                         % (i, circular, self.main_directory + i + '.pf',
                                            self.main_directory + i + '.fsa'))
        elif self.element_type == ":CONTIG":
            for i in self.regions_dict.keys():
                dat_file_str_list.append('ID\t%s\nchromosome_type\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                                         % (i, self.element_type, self.main_directory + i + '.pf',
                                            self.main_directory + i + '.fsa'))
        else:
            for i in self.regions_dict.keys():
                dat_file_str_list.append(
                    'ID\t%s\nchromosome_type\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n' % (
                        i, self.element_type, circular, self.main_directory + i + '.pf',
                        self.main_directory + i + '.fsa'))
        utils.write_file(self.directory + "genetic-elements" + ".dat", dat_file_str_list)

    def _make_fsa_files(self):
        with open(self.genomic_fasta_file_path, "r") as file:
            genomic_fasta = file.read()
        genomic_fasta = genomic_fasta.split(">")
        genomic_fasta = list(filter(None, genomic_fasta))
        regions_list = list(self.regions_dict.keys())
        for i in genomic_fasta:
            region = re.search("\w+(\.\w+)*(\-\w+)*", i).group(0)
            if region in regions_list:
                regions_list.remove(region)
                utils.write_file(self.directory + region + ".fsa", i, False)

    def _make_pf_files(self):
        tsv = utils.read_file_listed(self.eggnog_file_path)
        list_index = list(np.arange(0, len(tsv)))
        for region in self.regions_dict.keys():
            sub_pf = []
            for gene in self.regions_dict[region].keys():
                for protein in self.regions_dict[region][gene]["Proteins"].keys():
                    found = False
                    for i in list_index:
                        if found:
                            break
                        if protein in tsv[i]:
                            list_index.remove(i)
                            found = True
                            sub_pf.append(self._eggnog_file_parser(region, gene, protein, tsv[i]))
            if sub_pf:
                f = open(self.directory + region + ".pf", "w")
                for i in sub_pf:
                    for j in i:
                        f.write(j)
                f.close()

    def _eggnog_file_parser(self, region, gene, protein, line):
        """Sub-function of make_pf_files() to write the info in the correct order for each protein.

        PARAMS:
            region (str) -- actual genomic scope.
            gene (str) -- one gene id in the scope.
            protein (str) -- one of the proteins coded by the gene above.
            start (int) -- the start position of the sequence.
            line (str) -- the line corresponding to the protein in the .tsv file.
        RETURNS:
            info (str) -- a string with all the information and with the correct
            file architecture settings for the .pf file.
        """

        start = self.regions_dict[region][gene]["Start"]
        end = self.regions_dict[region][gene]["End"]
        cds_pos = self.regions_dict[region][gene]["Proteins"][protein]

        info = []
        spl = line.split("\t")
        info.append("ID\t" + gene + "\n")
        if spl[8] == "-":
            info.append("NAME\tORF\n")
        else:
            info.append("NAME\t" + spl[8] + "\n")
        info.append("STARTBASE\t" + start + "\n")
        info.append("ENDBASE\t" + end + "\n")
        try:
            spl[7] = spl[7].replace("\n", "")
            if spl[7]:
                info.append("FUNCTION\t" + spl[7] + "\n")
        except IndexError:
            info.append("FUNCTION\tNA\n")
        info.append("PRODUCT-TYPE\tP\n")
        if spl[10] != "-":
            for res in spl[10].split(","):
                info.append("EC\t" + res + "\n")
        for cds in cds_pos:
            info.append("CODING-SEGMENT\t" + str(cds[0]) + "-" + str(cds[1]) + "\n")
        if spl[9] == "-":
            pass
        else:
            go = spl[9].split(",")
            for i in go:
                info.append("DBLINK\t" + i + "\n")
        info.append("//\n")
        return info

    def build(self):
        print(self.name + " : Creating the .dat files...")
        self._make_dat_files()
        print(self.name + " : Creating the .fsa files...")
        self._make_fsa_files()
        print(self.name + " : Creating the .pf files...")
        self._make_pf_files()


def make_taxon_file(directory, taxon_name_list):
    """Function to make the taxon_id.tsv file.

    PARAMS:
        taxon_name_list (list) -- the list containing the name of the organism and its taxon id.
    """

    res = [["species", "taxon_id", "element_type"]]
    for i in taxon_name_list:
        res.append([i[0], str(i[1]), str(i[2])])
    utils.write_csv(directory, "taxon_id", res, separator="\t")


def mpwt_multirun_first(main_directory):
    """
    Split of major function 'run', first part = gathering the parameters, files and candidates' names and
    creating one Mpwting object for each.
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
        list_objects = []
        taxon_name_list = []
        cpu = len(parameters.keys()) - 1
        for i in parameters.keys():
            if i != "DEFAULT":
                logging.info("Species found : {}".format(i))
                species_name = parameters[i]["ORGANISM_NAME"]
                element_type = parameters[i]["ELEMENT_TYPE"]
                taxon_id = int(parameters[i]["NCBI_TAXON_ID"])
                species_directory = input_directory + species_name + "/"
                utils.make_directory(species_directory)
                taxon_name_list.append([species_name, taxon_id, element_type])
                list_objects.append(Mpwting(species_name, main_directory, element_type))
                make_taxon_file(input_directory, taxon_name_list)
        return [list_objects, cpu, input_directory, output_directory, log_directory]


def mpwt_multirun_last(list_objects, cpu, input_directory, output_directory, log_directory):
    """
    Split of major function 'run', second part = launching the process on each given object with multiprocessing
    and launching the mpwt reconstruction.
    """

    p = multiprocessing.Pool(cpu)
    p.map(build_mpwt_objects, list_objects)
    nb_cpu = multiprocessing.cpu_count()
    if nb_cpu <= cpu:
        cpu = nb_cpu - 1
    print("\n------\nNow launching MPWT on %s core(s)\n------" % cpu)
    mpwt.multiprocess_pwt(input_folder=input_directory, output_folder=output_directory, patho_inference=True,
                          patho_hole_filler=False, patho_operon_predictor=False, pathway_score=1, flat_creation=True,
                          dat_extraction=True, number_cpu=cpu, size_reduction=False, patho_log=log_directory,
                          taxon_file=input_directory + "taxon_id.tsv", verbose=True)


def build_mpwt_objects(organism):
    """Small function required for the multiprocessing reconstruction."""

    organism.build()


def run(main_directory):
    """The function to make all the run working."""

    mpwt_multirun_last(*mpwt_multirun_first(main_directory))


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
    run(utils.slash(args.main_directory))


if __name__ == "__main__":
    main()
