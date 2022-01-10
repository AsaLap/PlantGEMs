# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux - DRC
# Reconstruction de réseaux métaboliques
# Novembre 2021

"""This file is the parent class for all the modules, inherited by each of them to centralize repetitive methods."""


import os

import utils


class Module:

    def __init__(self, _name, _main_directory):
        self.name = _name
        self.main_directory = _main_directory.rstrip("/ ") + "/"

    # TODO : take multiple file format in account for the different "find" functions (cf. find_file in utils.py)
    # TODO : log of the file(s) used
    def _find_eggnog(self, target):
        return utils.find_file(self.main_directory + "/files/", target, ".tsv")

    def _find_genomic_fasta(self, target):
        return utils.find_file(self.main_directory + "/files/", target, ".fna")

    def _find_gff(self, target):
        return utils.find_file(self.main_directory + "/files/", target, ".gff")

    def _find_proteomic_fasta(self, target):
        return utils.find_file(self.main_directory + "/files/", target, ".faa")

    def _find_sbml_model(self, files_directory):
        # TODO : log of the file used
        model = utils.find_files(files_directory, "sbml")
        if len(model) == 0:
            print("No SBML file found in the files directory...")
            try:
                model = str(input("SBML path : "))
                if os.path.isfile(model):
                    return model
                else:
                    print("No SBML file found here... Restarting.")
                    self._find_sbml_model(files_directory)
            except ValueError:
                print("Please enter strings only... Restarting.")
                self._find_sbml_model(files_directory)
        elif len(model) >= 2:
            print("More than one SBML file has been found, please select one by entering its corresponding number :")
            for i in range(len(model)):
                print(i + 1, " : ", model[i])
            try:
                res = int(input("Chosen file number : "))
                return files_directory + model[res - 1]
            except IndexError:
                print("Please choose a valid number... Restarting.")
                self._find_sbml_model(files_directory)
            except ValueError:
                print("Please enter a number only... Restarting.")
                self._find_sbml_model(files_directory)
        else:
            return self.main_directory + "files/" + model[0]
