# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the metabolic reconstruction using sequence homology."""

import cobra
import copy
import os
import re
import subprocess
import sys
import time
import utils


class Blasting:

    def __init__(self, _name, _main_directory, _model_file_path=None, _model_fasta_path=None, _subject_fasta_path=None):
        """
        ARGS :
            _name -- name of the subject, must corresponds to the files' names.
            _main_directory -- main directory with the files et subdirectories for the results.
        (optional):
            _model_file_path -- the path to the SBML file containing the model for the reconstruction.
            _model_fasta_path -- the path to the fasta file of the model.
            _subject_fasta_path -- the path to the fasta file of the subject.
        """
        self.name = _name
        self.main_directory = _main_directory.rstrip("/ ") + "/"
        self._check_directories()
        if _model_file_path is not None:
            self.model = cobra.io.read_sbml_model(_model_file_path)
            self.model_fasta_path = _model_fasta_path
        else:
            self.model = cobra.io.read_sbml_model(self._find_model())
        if _model_fasta_path is not None:
            self.model_fasta = open(self.model_fasta_path).read()
        else:
            self.model_fasta = open(self._find_fasta(self.model.id)).read()
        if _subject_fasta_path is not None:
            self.subject_fasta_path = _subject_fasta_path
            self.subject_fasta = open(self.subject_fasta_path).read()
        else:
            self.subject_fasta_path = self._find_fasta(self.name)
            self.subject_fasta = open(self.subject_fasta_path).read()
        self.subject_directory = self.main_directory + "/" + self.model.id + "/"
        self.blast_result = {}
        self.gene_dictionary = {}
        self.identity = 50
        self.difference = 30
        self.e_val = 1e-100
        self.coverage = 20
        self.bit_score = 300

        """
            identity (int) -- the identity's threshold value to select the subject genes.
            difference (int) -- the percentage of length difference tolerated between subject and query.
            e_val (int) -- the minimum E-Value chosen.
            coverage (int) -- the minimum coverage's percentage of the match.
            bit_score (int) -- the minimum Bit-Score chosen.
        """

    @property
    def identity(self):
        print("Getting identity value...")
        return self.__identity

    @identity.setter
    def identity(self, value):
        if 100 >= value >= 0:
            print("Setting identity value to %s" % (str(value)))
            self.__identity = value
        else:
            print("Denied : value must be between 0 and 100 (both included)")

    @property
    def difference(self):
        print("Getting difference value...")
        return self.__difference

    @difference.setter
    def difference(self, value):
        if 100 >= value >= 0:
            print("Setting difference value to %s" % (str(value)))
            self.__difference = value
        else:
            print("Denied : value must be between 0 and 100 (both included)")

    @property
    def e_val(self):
        print("Getting E-value...")
        return self.__e_val

    @e_val.setter
    def e_val(self, value):
        if 10 >= value >= 0:
            print("Setting E-value to %s" % (str(value)))
            self.__e_val = value
        else:
            print("Denied : value must be between 0 and 10 (both included)")

    @property
    def coverage(self):
        print("Getting coverage value...")
        return self.__coverage

    @coverage.setter
    def coverage(self, value):
        if 100 >= value >= 0:
            print("Setting coverage value to %s" % (str(value)))
            self.__coverage = value
        else:
            print("Denied : value must be between 0 and 100 (both included)")

    @property
    def bit_score(self):
        print("Getting Bit-score value...")
        return self.__bit_score

    @bit_score.setter
    def bit_score(self, value):
        if 10000 >= value >= 0:
            print("Setting Bit-score value to %s" % (str(value)))
            self.__bit_score = value
        else:
            print("Denied : value must be between 0 and 10000 (both included)")

    def _check_directories(self):
        if not os.path.isdir(self.main_directory + "blast/"):
            print("Creation of necessary directory blast/")
            subprocess.run(["mkdir", self.main_directory + "blast/"])

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

    def _find_fasta(self, target):
        # TODO : take multiple file format in account
        # TODO : log of the file used
        fasta_path = self.main_directory + "/files/" + target + ".fasta"
        if os.path.isfile(self.main_directory + "/files/" + target + ".fasta"):
            return fasta_path
        else:
            print("No corresponding file found...")
            try:
                fasta_file_path = str(input("Path to " + target + " fasta file : "))
                if os.path.isfile(fasta_file_path):
                    return fasta_file_path
                else:
                    print("No file found with this path, make sure you entered it correctly... Restarting.")
                    self._find_fasta(target)
            except ValueError:
                print("Please enter a string only... Restarting.")
                self._find_fasta(target)

    def set_default_values(self):
        self.identity = 50
        self.difference = 30
        self.e_val = 1e-100
        self.coverage = 20
        self.bit_score = 300

    def _blast_run(self):
        """Runs multiple blasts between the model and the subject."""

        if not self.gene_dictionary:
            print("\nLaunching the blast !")
            i, x = 1, len(self.model.genes)
            total_time = lap_time = time.time()
            tmp_dir = self.main_directory + "/blast/tmp_dir/"
            try:
                subprocess.run(["rm -rf", tmp_dir])
            except FileNotFoundError:
                pass
            subprocess.run(["mkdir", tmp_dir])
            query_file = open(self.model_fasta_path)
            for seq in self.model_fasta.split(">"):
                try:
                    gene_name = re.search('\w+(\.\w+)*(-\w+)*', seq).group(0)
                    f = open(tmp_dir + gene_name + ".fa", "w")
                    f.write(">" + seq)
                    f.close()
                except AttributeError:
                    print("Gene name not found in :", seq)
                    pass
            query_file.close()
            for gene in self.model.genes:
                if i % 10 == 0:
                    print("Protein %i out of %i\nTime : %f s" % (i, x, time.time() - lap_time))
                    lap_time = time.time()
                i += 1
                blast_request = [
                    "blastp",
                    "-subject",
                    self.subject_fasta_path,
                    "-query",
                    tmp_dir + gene.id + ".fa",
                    "-outfmt",
                    "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
                self.blast_result[gene.id] = subprocess.run(blast_request,
                                                            capture_output=True).stdout.decode('ascii').split("\n")[:-1]
            try:
                subprocess.run(["rm", "-rf", tmp_dir])
            except FileNotFoundError:
                print("Temporary folder not found, not erased...")
            print("Blast done !\nTotal time : %f s" % (time.time() - total_time))

    def _select_genes(self):
        """Select the subject organism's genes regarding the different treshold parameters of the Blasting instance."""

        if not self.blast_result:
            print("No blast results found... Please run a blast with blast_run() before launching select_genes()")
        else:
            for key in self.blast_result.keys():
                for res in self.blast_result[key]:
                    spl = res.split(",")
                    for i in range(len(spl)):
                        try:
                            spl[i] = float(spl[i])
                        except ValueError:
                            pass
                    len_subject = spl[3]
                    len_query = [spl[1] * (100 - self.difference) / 100, spl[1] * (100 + self.difference) / 100]
                    min_align = self.coverage / 100 * spl[1]
                    if spl[6] >= self.identity \
                            and len_query[0] <= len_subject <= len_query[1] \
                            and spl[4] >= min_align \
                            and spl[9] >= self.bit_score:
                        try:
                            if spl[8] <= self.e_val:
                                self.gene_dictionary[key].append(spl[2])
                        except KeyError:
                            if spl[8] <= self.e_val:
                                self.gene_dictionary[key] = [spl[2]]

    def _drafting(self):
        """Creates the new COBRA model for the subject organism."""

        self.draft = cobra.Model(self.name)
        for reaction in self.model.reactions:
            to_add = []
            to_search = reaction.gene_reaction_rule.split(" or ")
            for gene in to_search:
                try:
                    to_add += self.gene_dictionary[gene]
                except KeyError:
                    pass
            # TODO : change the proteins' ID in to_add in corresponding genes
            string_reaction_rule = " or ".join(to_add)
            if string_reaction_rule:
                x = copy.deepcopy(reaction)
                x.gene_reaction_rule = string_reaction_rule
                self.draft.add_reactions([x])

    def _history_save(self, step):
        history_directory = self.main_directory + "blast/blast_object_history/"
        if not os.path.isdir(history_directory):
            try:
                subprocess.run(["mkdir", history_directory])
            except PermissionError:
                print("Permission to create this folder :\n" + history_directory + "\nnot granted !")
        utils.save_obj(self, history_directory + self.name + "_" + step)

    # TODO : create a loading function
    # TODO : create a log function to gather all the errors encountered

    def build(self):
        self._blast_run()
        self._history_save("blasted")
        self._select_genes()
        self._history_save("genes_selected")
        self._drafting()
        self._history_save("drafted")
        cobra.io.save_json_model(self.draft, self.subject_directory + self.name + "_blast.json")


def pipeline(*args):
    """Function to use this script in CLI."""
    try:
        cli_blast = Blasting(*args)
        cli_blast.build()
    except TypeError:
        print("Usage : $ python blasting.py pipeline parameter1 parameter2 parameter3 parameter4\nParameters "
              "required : name, main_directory"
              "optional : model_file_path, model_fasta_path, subject_fasta_path")


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
