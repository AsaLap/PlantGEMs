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

    def __init__(self, _name, _main_directory, _model_file_path=None, _model_proteomic_fasta_path=None,
                 _subject_proteomic_fasta_path=None, _subject_gff_path=None):
        """
        ARGS :
            _name -- name of the subject, must corresponds to the files' names.
            _main_directory -- main directory with the files et subdirectories for the results.
        (optional):
            _model_file_path -- the path to the SBML file containing the model for the reconstruction.
            _model_proteomic_fasta_path -- the path to the fasta file of the model.
            _subject_proteomic_fasta_path -- the path to the fasta file of the subject.
        """
        super().__init__(_name, _main_directory)
        utils.make_directory(self.main_directory + "blast/")
        if _model_file_path is not None:
            self.model = cobra.io.read_sbml_model(_model_file_path)
        else:
            self.model = cobra.io.read_sbml_model(self._find_model())
        if _model_proteomic_fasta_path is not None:
            self.model_proteomic_fasta_path = _model_proteomic_fasta_path
            self.model_proteomic_fasta = utils.read_file_stringed(self.model_proteomic_fasta_path)
        else:
            self.model_proteomic_fasta = utils.read_file_stringed(self._find_proteomic_fasta(self.model.id))
        if _subject_proteomic_fasta_path is not None:
            self.subject_proteomic_fasta_path = _subject_proteomic_fasta_path
            self.subject_proteomic_fasta = utils.read_file_stringed(self.subject_proteomic_fasta_path)
        else:
            self.subject_proteomic_fasta_path = self._find_proteomic_fasta(self.name)
            self.subject_proteomic_fasta = utils.read_file_stringed(self.subject_proteomic_fasta_path)
        if _subject_gff_path is not None:
            self.subject_gff_path = _subject_gff_path
        else:
            self.gff_file_path = self._find_gff(self.name)
        self.directory = self.main_directory + "blast/" + self.name + "/"
        self.regions_dict = utils.get_sequence_region(self.gff_file_path)
        self._make_protein_correspondence_file()
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
            tmp_dir = self.directory + "tmp_dir/"
            utils.remove_directory(tmp_dir)
            utils.make_directory(tmp_dir)
            for seq in self.model_proteomic_fasta.split(">"):
                if seq:
                    try:
                        gene_name = re.search('\w+(\.\w+)*(-\w+)*', seq).group(0)
                        utils.write_file(tmp_dir + gene_name + ".fa", [">" + seq])
                    except AttributeError:
                        print("Gene name not found in :", seq)  # TODO : Log this
                        pass
            for gene in self.model.genes:
                if i % 10 == 0:
                    print("\n" + self.name + " : Protein %i out of %i\nTime : %f s\n" % (i, x, time.time() - lap_time))
                    lap_time = time.time()
                i += 1
                blast_request = [
                    "blastp",
                    "-subject",
                    self.subject_proteomic_fasta_path,
                    "-query",
                    tmp_dir + gene.id + ".fa",
                    "-outfmt",
                    "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
                self.blast_result[gene.id] = subprocess.run(blast_request,
                                                            capture_output=True).stdout.decode('ascii').split("\n")[:-1]
            utils.remove_directory(tmp_dir)
            print(self.name + " : Blast done !\nTotal time : %f s" % (time.time() - total_time))

    def _select_genes(self):
        """Select the subject organism's genes regarding the different threshold parameters of the Blasting instance."""

        if not self.blast_result:
            print("No blast results found... Please run a blast with blast_run() before launching select_genes()")
        else:
            for key in self.blast_result.keys():  # key = region name
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
                            and spl[9] >= self.bit_score \
                            and spl[8] <= self.e_val:
                        try:
                            self.gene_dictionary[key].append(spl[2])
                        except KeyError:
                            self.gene_dictionary[key] = [spl[2]]

    def _drafting(self):
        """Creates the new COBRA model for the subject organism."""

        self.draft = cobra.Model(self.name)
        for reaction in self.model.reactions:
            to_add = []
            for gene in reaction.gene_reaction_rule.split(" or "):
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
        utils.make_directory(self.directory)
        self._blast_run()
        self._history_save("blasted")
        self._select_genes()
        self._history_save("genes_selected")
        self._drafting()
        self._history_save("drafted")
        self._protein_to_gene()

    def _make_protein_correspondence_file(self):
        """Function to create a csv file with the correspondence between a protein and the associated gene."""

        correspondence = []
        for region in self.regions_dict.keys():
            for gene in self.regions_dict[region].keys():
                for protein in self.regions_dict[region][gene]["Proteins"].keys():
                    correspondence.append([gene.upper(), protein.upper()])
        utils.write_csv(self.main_directory + "/blast/" + self.name, "protein_gene_correspondence", correspondence,
                        "\t")

    def _protein_to_gene(self):  # TODO : Review this code and use it in the pipeline
        """Function to transform the proteins in gene_reaction_rule into its corresponding genes.
        It creates a new model that will have all the genes' names instead of the proteins' ones."""

        self._make_protein_correspondence_file()
        correspondence_file_path = self.directory + "protein_gene_correspondence.tsv"
        if os.path.isfile(correspondence_file_path):
            gene_protein_dict, gene_protein_dict_reversed = utils.build_correspondence_dict(correspondence_file_path)
            for reaction in self.draft.reactions:
                genes = []
                list_proteins = reaction.gene_reaction_rule.split(" or ")
                for protein in list(filter(None, [trans.upper() for trans in list_proteins])):
                    try:
                        genes.append(gene_protein_dict_reversed[protein])
                    except KeyError:
                        try:
                            if protein in gene_protein_dict.keys():
                                genes.append(protein)
                        except KeyError:
                            print("No match for : ", protein)
                reaction.gene_reaction_rule = " or ".join(set(genes))
            cobra.io.save_json_model(self.draft, self.directory + self.name + "_gene_blast" + ".json")
        else:
            print("No correspondence file found here : " + correspondence_file_path + "\nAborting...")
            sys.exit()


def build_blast_objects(organism_object):
    """Small function required for the multiprocessing reconstruction."""

    organism_object.build()


def sub_pipeline_first(main_directory):
    """
    Split of major function 'pipeline', first part = gathering the files and candidates' names and creating one
    Blasting object for each.
    """

    parameters = utils.read_config(main_directory + "main.ini")
    if os.path.isdir(main_directory):
        list_objects = []
        for i in parameters.keys():
            if i != "DEFAULT":
                list_objects.append(Blasting(parameters[i]["ORGANISM_NAME"], main_directory))
    else:
        sys.exit("Main directory given does not exist : " + main_directory)
    return list_objects


def sub_pipeline_last(list_objects):
    """
    Split of major function 'pipeline', second part = launching the process on each given object with multiprocessing.
    """

    cpu = len(list_objects)
    p = multiprocessing.Pool(cpu)
    p.map(build_blast_objects, list_objects)


def pipeline(main_directory):
    """The function to launch the process when used alone."""

    list_objects = sub_pipeline_first(main_directory)
    sub_pipeline_last(list_objects)


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
              "optional : model_file_path, model_proteomic_fasta_path,"
              "subject_proteomic_fasta_path, subject_gff_path\n")


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
