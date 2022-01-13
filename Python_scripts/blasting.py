# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the metabolic reconstruction using sequence homology."""

import argparse
import cobra
import copy
import getopt
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
                 _subject_proteomic_fasta_path=None, _subject_gff_path=None,
                 identity=50, difference=30, e_val=1e-100, coverage=20, bit_score=300):
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
            self.model = cobra.io.read_sbml_model(self._find_sbml_model(self.main_directory + "/files/"))
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
        self.blast_result = {}
        self.gene_dictionary = {}
        self._identity = identity
        self._difference = difference
        self._e_val = e_val
        self._coverage = coverage
        self._bit_score = bit_score

        """
        identity (int) -- the blast result identity's threshold value to select the subject genes.
        difference (int) -- the percentage of length difference tolerated between subject and query.
        e_val (int) -- the minimum E-Value of each match.
        coverage (int) -- the minimum sequence coverage of the match (percentage).
        bit_score (int) -- the minimum Bit-Score of each match.
        """

    @property
    def identity(self):
        return self._identity

    @identity.setter
    def identity(self, value):
        if 100 >= value >= 0:
            print("Setting identity value to %s" % (str(value)))
            self._identity = value
        else:
            print("Identity value denied : value must be between 0 and 100 (both included), value not changed")

    @property
    def difference(self):
        return self._difference

    @difference.setter
    def difference(self, value):
        if 100 >= value >= 0:
            print("Setting difference value to %s" % (str(value)))
            self._difference = value
        else:
            print("Difference value denied : value must be between 0 and 100 (both included), value not changed")

    @property
    def e_val(self):
        return self._e_val

    @e_val.setter
    def e_val(self, value):
        if 10 >= value >= 0:
            print("Setting E-value to %s" % (str(value)))
            self._e_val = value
        else:
            print("E_val value denied : value must be between 0 and 10 (both included), value not changed")

    @property
    def coverage(self):
        return self._coverage

    @coverage.setter
    def coverage(self, value):
        if 100 >= value >= 0:
            print("Setting coverage value to %s" % (str(value)))
            self._coverage = value
        else:
            print("Coverage value denied : value must be between 0 and 100 (both included), value not changed")

    @property
    def bit_score(self):
        return self._bit_score

    @bit_score.setter
    def bit_score(self, value):
        if 10000 >= value >= 0:
            print("Setting Bit-score value to %s" % (str(value)))
            self._bit_score = value
        else:
            print("Bit_score value denied : value must be between 0 and 10000 (both included), value not changed")

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
                    to_add += self.gene_dictionary[gene]  # TODO : faire le changement de protein/gene ici
                except KeyError:
                    pass
            string_reaction_rule = " or ".join(to_add)
            if string_reaction_rule:
                x = copy.deepcopy(reaction)
                x.gene_reaction_rule = string_reaction_rule
                self.draft.add_reactions([x])

    def _object_history_save(self, step):
        objects_directory = self.directory + "/objects_history/"
        utils.make_directory(objects_directory)
        utils.save_obj(self, objects_directory + step)

    def _make_protein_correspondence_file(self):
        """Function to create a csv file with the correspondence between a protein and the associated gene."""

        correspondence = []
        for region in self.regions_dict.keys():
            for gene in self.regions_dict[region].keys():
                for protein in self.regions_dict[region][gene]["Proteins"].keys():
                    correspondence.append([gene.upper(), protein.upper()])
        utils.write_csv(self.directory, "protein_gene_correspondence", correspondence, "\t")

    def _protein_to_gene(self):  # TODO : Review this code and use it in the run
        """Function to transform the proteins in gene_reaction_rule into their corresponding genes.
        It creates a new model that will have all the genes' names instead of the proteins' ones."""

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
        else:
            print("No correspondence file found here : " + correspondence_file_path + "\nAborting...")
            sys.exit()

    def build(self):
        utils.make_directory(self.directory)
        self._make_protein_correspondence_file()
        self._blast_run()
        self._object_history_save("blasted")
        self._select_genes()
        self._object_history_save("genes_selected")
        self._drafting()
        self._object_history_save("drafted")
        self._protein_to_gene()
        cobra.io.save_json_model(self.draft, self.directory + self.name + "_blast_draft" + ".json")

    def rebuild(self):
        self._select_genes()
        self._object_history_save("genes_selected")
        self._drafting()
        self._object_history_save("drafted")
        self._protein_to_gene()
        cobra.io.save_json_model(self.draft, self.directory + self.name + "_blast_draft_rebuild_" + "_".join(
            (str(self.identity), str(self.difference), str(self.e_val), str(self.coverage), str(self._bit_score))) + ".json")


def blast_multirun_first(main_directory):
    """
    Split of major function 'run', first part = gathering the files and candidates' names and creating one
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


def blast_multirun_last(list_objects):
    """
    Split of major function 'run', second part = launching the process on each given object with multiprocessing.
    """

    cpu = len(list_objects)
    p = multiprocessing.Pool(cpu)
    p.map(build_blast_objects, list_objects)


def build_blast_objects(organism_object):
    """Small function required for the multiprocessing reconstruction."""

    organism_object.build()


def run(main_directory):
    """The function to launch the process when used alone."""

    list_objects = blast_multirun_first(main_directory)
    blast_multirun_last(list_objects)


def run_unique(*args):
    """
    This function allows to launch a blasting process on a unique organism with every argument in command line
    if wished so.
    """

    try:
        unique_blast = Blasting(*args)
        unique_blast.build()
    except TypeError:
        print("Usage : $ python blasting.py run 'parameter1=value parameter2=value parameter3=value...'\n"
              "Parameters (in order) : \n"
              "required : name, main_directory\n"
              "optional : model_file_path, model_proteomic_fasta_path,"
              "subject_proteomic_fasta_path, subject_gff_path\n")


def rerun_blast_selection(blasted_object, identity=50, difference=30, e_val=1e-100, coverage=20, bit_score=300):
    species = utils.load_obj(blasted_object)
    species.identity = int(identity)
    species.difference = int(difference)
    species.e_val = float(e_val)
    species.coverage = int(coverage)
    species.bit_score = int(bit_score)
    species.rebuild()


# def usage():
#     print("blasting.py -n <name> -md <main directory>")


# def main(argv):
#     short_options = "d:hi:"
#     long_options = ["directory=", "help", "identity="]
#     try:
#         opts, args = getopt.getopt(argv, short_options, long_options)
#     except getopt.GetoptError as err:
#         print(err)
#         usage()
#         sys.exit(2)
#     main_directory = ""
#     identity = 50
#     for opt, arg in opts:
#         if opt == "-h":
#             usage()
#         if opt == "-d":
#             main_directory = arg
#         if opt in ("-i", "--identity"):
#             identity = arg

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_directory", type=str, help="The path to the main directory")
    parser.add_argument("-i", "--identity", help="The the blast's identity percentage",
                        action="store_false")
    args = parser.parse_args()
    print(args.main_directory)


if __name__ == "__main__":
    # globals()[sys.argv[1]](*sys.argv[2:])
    # run(sys.argv[1:])
    main()
