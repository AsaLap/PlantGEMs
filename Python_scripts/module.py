# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux - DRC
# Reconstruction de réseaux métaboliques
# Novembre 2021

"""This file is the parent class for all the modules, inherited by each of them to centralize repetitive methods."""


import os


class Module:

    def __init__(self, _name, _main_directory):
        self.name = _name
        self.main_directory = _main_directory.rstrip("/ ") + "/"

    def _find_file(self, target, extension):
        """Search a file corresponding to the target in the 'files' directory. If no match is found,
        asks the user to input the exact path to the file he wants to use.

        PARAMS:
            target (str) -- the name of the file to search if correctly named.
            extension (str) -- the file's extension.

        RETURNS:
            file_path (str) -- the exact path to the file.
        """
        # TODO : take multiple file format in account
        # TODO : log of the file used

        file_path = self.main_directory + "/files/" + target + extension
        if os.path.isfile(file_path):
            return file_path
        else:
            print("No corresponding file found...")
            try:
                file_path = str(input("Path to " + target + extension + "'s file : "))
                if os.path.isfile(file_path):
                    return file_path
                else:
                    print("No file found with this path, make sure you entered it correctly... Restarting.")
                    self._find_file(target, extension)
            except ValueError:
                print("Please enter a string only... Restarting.")
                self._find_file(target, extension)

    def _find_genomic_fasta(self, target):
        return self._find_file(target, ".fna")

    def _find_proteomic_fasta(self, target):
        return self._find_file(target, ".faa")

    def _find_gff(self, target):
        return self._find_file(target, ".gff")

    def _find_eggnog(self, target):
        return self._find_file(target, ".tsv")

    def _find_model(self):
        # TODO : log of the file used
        files_directory = self.main_directory + "files/"
        model = [i for i in os.listdir(files_directory) if i.endswith("sbml")]
        if len(model) == 0:
            print("No SBML file found in the files directory...")
            try:
                model = str(input("SBML path : "))
                if os.path.isfile(model):
                    return model
                else:
                    print("No SBML file found here... Restarting.")
                    self._find_model()
            except ValueError:
                print("Please enter strings only... Restarting.")
                self._find_model()
        elif len(model) >= 2:
            print("More than one SBML file has been found, please select one by entering its corresponding number :")
            for i in range(len(model)):
                print(i + 1, " : ", model[i])
            try:
                res = int(input("Chosen file number : "))
                return files_directory + model[res - 1]
            except IndexError:
                print("Please choose a valid number... Restarting.")
                self._find_model()
            except ValueError:
                print("Please enter a number only... Restarting.")
                self._find_model()
        else:
            return self.main_directory + "files/" + model[0]