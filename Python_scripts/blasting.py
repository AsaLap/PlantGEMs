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

    def __init__(self, _name, _model, _model_fasta_path, _subject_fasta_path):
        """
        ARGS:
            _model -- the path to the SBML file containing the model for the reconstruction.
            _model_fasta_path -- the path to the fasta file of the model.
            _subject_fasta_path -- the path to the fasta file of the subject.
        """
        self.name = _name
        self.model = cobra.io.read_sbml_model(_model)
        self.model_fasta_path = _model_fasta_path
        self.model_directory = os.path.split(self.model_fasta_path)[0]
        self.model_fasta = open(self.model_fasta_path).read()
        self.subject_fasta_path = _subject_fasta_path
        self.subject_fasta = open(self.subject_fasta_path).read()
        self.subject_directory = os.path.split(self.subject_fasta_path)[0]
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
            print("\nLaunching the blast !")
            i, x = 1, len(self.model.genes)
            total_time = lap_time = time.time()
            tmp_dir = self.model_directory + "/tmp_dir/"
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
        history_dir = self.subject_directory + "/blast_history/"
        if not os.path.isdir(history_dir):
            try:
                subprocess.run(["mkdir", history_dir])
            except PermissionError:
                print("Permission to create this folder :\n" + history_dir + "\nnot granted !")
        utils.save_obj(self, history_dir + self.name + "_" + step)

    # TODO : create a loading function

    def build(self):

        self._blast_run()
        self._history_save("blast")
        self._select_genes()
        self._history_save("gene_selected")
        self._drafting()
        self._history_save("drafted")
        cobra.io.save_json_model(self.draft, self.subject_directory + "/" + self.name + "_blast.json")


def cli_blasting(*args):
    """Function to use this script in CLI."""
    try:
        cli_blast = Blasting(*args)
        cli_blast.build()
    except TypeError:
        print("Usage : $ python blasting.py cli_blasting parameter1 parameter2 parameter3 parameter4\nParameters "
              "required : name, model, model_fasta_path, subject_fasta_path")


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
