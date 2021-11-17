# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the metabolic reconstruction using sequence homology."""

import cobra
import copy
import multiprocessing
import os
import re
import subprocess
import sys
import time

import module
import utils


class Blasting(module.Module):

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
        super().__init__(_name, _main_directory)
        utils.make_directory(self.main_directory + "blast/")
        if _model_file_path is not None:
            self.model = cobra.io.read_sbml_model(_model_file_path)
        else:
            self.model = cobra.io.read_sbml_model(self._find_model())
        if _model_fasta_path is not None:
            self.model_fasta_path = _model_fasta_path
            self.model_fasta = utils.read_file_stringed(self.model_fasta_path)
        else:
            self.model_fasta = utils.read_file_stringed(self._find_fasta(self.model.id))
        if _subject_fasta_path is not None:
            self.subject_fasta_path = _subject_fasta_path
            self.subject_fasta = utils.read_file_stringed(self.subject_fasta_path)
        else:
            self.subject_fasta_path = self._find_fasta(self.name)
            self.subject_fasta = utils.read_file_stringed(self.subject_fasta_path)
        self.subject_directory = self.main_directory + "blast/" + self.name + "/"
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

    def set_default_values(self):
        self.identity = 50
        self.difference = 30
        self.e_val = 1e-100
        self.coverage = 20
        self.bit_score = 300

    def _blast_run(self):
        """Runs multiple blasts between the model and the subject."""

        if not self.gene_dictionary:
            print(self.name + " : Launching the blast !")
            i, x = 1, len(self.model.genes)
            total_time = lap_time = time.time()
            tmp_dir = self.subject_directory + "tmp_dir/"
            utils.remove_directory(tmp_dir)
            utils.make_directory(tmp_dir)
            for seq in self.model_fasta.split(">"):
                try:
                    gene_name = re.search('\w+(\.\w+)*(-\w+)*', seq).group(0)
                    utils.write_file(tmp_dir + gene_name + ".fa", [">" + seq])
                except AttributeError:
                    print("Gene name not found in :", seq)
                    pass
            for gene in self.model.genes:
                if i % 10 == 0:
                    print(self.name + " : Protein %i out of %i\nTime : %f s" % (i, x, time.time() - lap_time))
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
            utils.remove_directory(tmp_dir)
            print(self.name + " : Blast done !\nTotal time : %f s" % (time.time() - total_time))

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
        utils.make_directory(history_directory)
        utils.save_obj(self, history_directory + self.name + "_" + step)

    def build(self):
        utils.make_directory(self.subject_directory)
        self._blast_run()
        self._history_save("blasted")
        self._select_genes()
        self._history_save("genes_selected")
        self._drafting()
        self._history_save("drafted")
        cobra.io.save_json_model(self.draft, self.subject_directory + self.name + "_blast.json")


def build_blast_objects(data):
    data.build()


def pipeline(main_directory):
    """The function to make all the pipeline working."""

    main_parameters = utils.read_config(main_directory + "main.ini")
    if os.path.isdir(main_directory):
        list_objects = []
        for i in main_parameters.keys():
            if i != "DEFAULT":
                list_objects.append(Blasting(main_parameters[i]["ORGANISM_NAME"], main_directory))
        cpu = len(list_objects)
        p = multiprocessing.Pool(cpu)
        p.map(build_blast_objects, list_objects)


def pipeline_unique(*args):
    """
    This function allows to launch a blasting process on a unique organism with every argument in command line
    if wished so.
    """

    try:
        unique_blast = Blasting(*args)
        unique_blast.build()
    except TypeError:
        print("Usage : $ python blasting.py pipeline parameter1=value parameter2=value parameter3=value...\n"
              "Parameters (in order) : \n"
              "required : name, main_directory\n"
              "optional : model_file_path, model_fasta_path, subject_fasta_path\n")


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
