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
import logging
import multiprocessing
import os
import re
import subprocess
import sys
import time

import graphs
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
            _subject_gff_path -- the path to the gff file of the subject.
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
        self.regions_dict = utils.get_sequence_region(self.gff_file_path)
        self.blast_result = {}
        self.gene_dictionary = {}
        self.draft = cobra.Model(self.name)
        self._identity = identity
        self._difference = difference
        self._e_val = e_val
        self._coverage = coverage
        self._bit_score = bit_score
        self.version = 1.0

        """
        identity (int) -- the blast result identity's threshold value to select the subject genes.
        difference (int) -- the percentage of length difference tolerated between subject and query.
        e_val (int) -- the minimum E-Value of each match.
        coverage (int) -- the minimum sequence coverage of the match (percentage).
        bit_score (int) -- the minimum Bit-Score of each match.
        """

    @property
    def directory(self):
        return self.main_directory + "blast/" + self.name + "/"

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
                        logging.info("Gene name not found in :", seq)
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
            log_message = self.name + " : Blast done !\nTotal time : %f s" % (time.time() - total_time)
            logging.info(log_message)
            print(log_message)

    def _select_genes(self):
        """Select the subject organism's genes regarding the different threshold parameters of the Blasting instance."""

        if not self.blast_result:
            logging.info(self.name + " : No blast results found... Please run a blast with blast_run() before launching"
                                     " select_genes()")
            print("No blast results found... Please run a blast with blast_run() before launching select_genes()")
        else:
            wrong_identity, wrong_difference, wrong_coverage, wrong_bit_score, wrong_e_val = [], [], [], [], []
            selected_proteins = [["Protein Model\tSize P. Model\tProtein Subject\tSize P. Subject\tAlignment length\t"
                                  "Number of identity\tPercentage of identity\tScore\tEValue\tBitScore"]]
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
                    selected = True
                    if spl[6] < self.identity:
                        wrong_identity.append(res)
                        selected = False
                    if not len_query[0] <= len_subject <= len_query[1]:
                        wrong_difference.append(res)
                        selected = False
                    if spl[4] < min_align:
                        wrong_coverage.append(res)
                        selected = False
                    if spl[9] < self.bit_score:
                        wrong_bit_score.append(res)
                        selected = False
                    if spl[8] > self.e_val:
                        wrong_e_val.append(res)
                        selected = False
                    if selected:
                        selected_proteins.append([str.replace(res, ",", "\t")])
                        try:
                            self.gene_dictionary[key].append(spl[2])
                        except KeyError:
                            self.gene_dictionary[key] = [spl[2]]
            removed_proteins_upsetplot_dict = {"Identity": wrong_identity,
                                               "Difference": wrong_difference,
                                               "Coverage": wrong_coverage,
                                               "Bit_Score": wrong_bit_score,
                                               "E_Value": wrong_e_val}
            graphs.make_upsetplot(self.directory, "removed_proteins_plot", removed_proteins_upsetplot_dict,
                                  "Thresholds responsible for unselected proteins")
            utils.write_csv(self.directory, "selected_proteins", selected_proteins)

    def _drafting(self):
        """Creates the new COBRA model for the subject organism."""

        for reaction in self.model.reactions:
            to_add = []
            for gene in reaction.gene_reaction_rule.split(" or "):
                try:
                    to_add += self.gene_dictionary[gene]
                except KeyError:
                    logging.info("{} gene led to a KeyError, therefore it wasn't selected as a candidate gene ({} model"
                                 " for {} subject)".format(gene, self.model.id, self.name))
                    pass
            string_reaction_rule = " or ".join(to_add)
            if string_reaction_rule:
                x = copy.deepcopy(reaction)
                x.gene_reaction_rule = string_reaction_rule
                self.draft.add_reactions([x])

    def _object_history_save(self, step):
        objects_directory = self.directory + "objects_history/"
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

    def _protein_to_gene(self):
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
                            log_message = self.name + " : No match for : ", protein
                            logging.error(log_message)
                            print(log_message)
                reaction.gene_reaction_rule = " or ".join(set(genes))
        else:
            log_message = self.name + " : No correspondence file found here : " + correspondence_file_path + \
                          "\nAborting..."
            logging.info(log_message)
            print(log_message)
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
        surname = "_".join((str(self.identity), str(self.difference), str(self.e_val), str(self.coverage),
                            str(self._bit_score)))
        self._select_genes()
        self._object_history_save("genes_selected_" + surname)
        self._drafting()
        self._object_history_save("drafted_" + surname)
        self._protein_to_gene()
        cobra.io.save_json_model(self.draft, self.directory + self.name + "_blast_draft_rebuild_" + surname + ".json")


def blast_multirun_first(args):
    """
    Split of major function 'run', first part = gathering the files and candidates' names and creating one
    Blasting object for each.
    """

    parameters = utils.read_config(args.main_directory + "main.ini")
    if os.path.isdir(args.main_directory):
        list_objects = []
        for i in parameters.keys():
            if i != "DEFAULT":
                logging.info("Parameters for : " + parameters[i]["ORGANISM_NAME"] +
                             "\n - Main directory : " + args.main_directory +
                             "\n - Identity : " + str(args.identity) +
                             "\n - Difference : " + str(args.difference) +
                             "\n - E_Value : " + str(args.e_val) +
                             "\n - Coverage : " + str(args.coverage) +
                             "\n - Bit_Score : " + str(args.bit_score))
                list_objects.append(Blasting(parameters[i]["ORGANISM_NAME"], args.main_directory,
                                             identity=args.identity, difference=args.difference, e_val=args.e_val,
                                             coverage=args.coverage, bit_score=args.bit_score))
    else:
        log_message = "Main directory given does not exist : " + args.main_directory
        logging.error(log_message)
        sys.exit(log_message)
    return list_objects


def blast_multirun_last(list_objects):
    """
    Split of major function 'run', second part = launching the process on each given object with multiprocessing.
    """

    cpu = len(list_objects)
    logging.info("Launching %i processes with multiprocess" % cpu)
    p = multiprocessing.Pool(cpu)
    p.map(build_blast_objects, list_objects)


def build_blast_objects(organism_object):
    """Small function required for the multiprocessing reconstruction."""

    organism_object.build()


def run(args):
    """The function to launch the process when used alone."""

    logging.info("\n------ Running multiple species or unique but with a config file ------")
    logging.info("Reading parameters...")
    list_objects = blast_multirun_first(args)
    logging.info("Launching the blast(s) with given parameters...")
    blast_multirun_last(list_objects)


def run_unique(args):
    """
    This function allows to launch a blasting process on a unique organism with every argument in command line
    if wished so.
    """

    logging.info("\n------ Running a unique species ------")
    logging.info("\nParameters for : {}\n - Main directory : {}\n - Model's file's path : {}\n - Model's proteomic "
                 "fasta's path : {}\n - Subject's proteomic fasta's path : {}\n - Subject's gff file's path : {}\n"
                 " - Identity : {}\n - Difference : {}\n - E_Value : {}\n - Coverage : {}\n - Bit_Score : {}"
                 .format(args.name, args.main_directory, args.model_file_path, args.model_proteomic_fasta_path,
                         args.subject_proteomic_fasta_path, args.subject_gff_path,
                         args.identity, args.difference, args.e_val, args.coverage, args.bit_score))
    unique_blast = Blasting(args.name, args.main_directory, args.model_file_path, args.model_proteomic_fasta_path,
                            args.subject_proteomic_fasta_path, args.subject_gff_path,
                            args.identity, args.difference, args.e_val, args.coverage, args.bit_score)
    unique_blast.build()


def rerun_blast_selection(main_directory, name, identity=50, difference=30, e_val=1e-100, coverage=20, bit_score=300):
    logging.info("\n------ Rerunning a species' genes selection ------")
    species = utils.load_obj(utils.slash(main_directory) + "blast/" + name + "/objects_history/blasted.pkl")
    logging.info("Parameters for : {}\n - Main directory : {}\n - Identity : {} -> {}\n - Difference : {} -> {}\n"
                 " - E_Value : {} -> {}\n - Coverage : {} -> {}\n - Bit_Score : {} -> {}"
                 .format(species.name, species.main_directory, species.identity, identity, species.difference,
                         difference, species.e_val, e_val, species.coverage, coverage, species.bit_score, bit_score))
    species.main_directory = main_directory
    species.identity = identity
    species.difference = difference
    species.e_val = e_val
    species.coverage = coverage
    species.bit_score = bit_score
    species.gene_dictionary = {}
    species.draft = cobra.Model(species.name)
    species.rebuild()


def blast_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_directory", help="The path to the main directory where the 'files/' directory is stored",
                        type=str)
    parser.add_argument("-v", "--verbose", help="Toggle the printing of more information", action="store_true")
    parser.add_argument("-le", "--log_erase", help="Erase the existing log file to create a brand new one",
                        action="store_true")
    parser.add_argument("-u", "--unique", help="Specify if the reconstruction is made on a unique species or not",
                        action="store_true")
    parser.add_argument("-rr", "--rerun", help="Use this option if you want to rerun the blast selection on an existing"
                                               " blasted.pkl object. The species' name is expected here", type=str)
    parser.add_argument("-n", "--name", help="The future draft's name", type=str)
    parser.add_argument("-m", "--model_file_path", help="Model's file's path, use if 'files/' directory doesn't exist",
                        type=str)
    parser.add_argument("-mfaa", "--model_proteomic_fasta_path",
                        help="Model's proteomic fasta's path, use if 'files/' directory doesn't exist")
    parser.add_argument("-sfaa", "--subject_proteomic_fasta_path",
                        help="Subject's proteomic fasta's path, use if 'files/' directory doesn't exist")
    parser.add_argument("-sgff", "--subject_gff_path",
                        help="Subject's gff file's path, use if 'files/' directory doesn't exist")
    parser.add_argument("-i", "--identity", help="The blast's identity percentage tolerated. Default=50",
                        type=int, default=50, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-d", "--difference",
                        help="The tolerated length difference between the two aligned sequences. Default=30",
                        type=int, default=30, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-ev", "--e_val",
                        help="The blast's e-value threshold value. Default=e-100",
                        type=float, default=1e-100, choices=range(0, 2), metavar="[0-1]")
    parser.add_argument("-c", "--coverage", help="The minimum sequence coverage tolerated. Default=20",
                        type=int, default=20, choices=range(0, 101), metavar="[0-100]")
    parser.add_argument("-bs", "--bit_score", help="The blast's bit-score threshold value. Default=300",
                        type=int, default=300, choices=range(0, 1001), metavar="[0-1000]")
    args = parser.parse_args()
    return args


def main():
    args = blast_arguments()
    if args.log_erase:
        logging.basicConfig(filename=args.main_directory + '/blasting.log', filemode='w', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename=args.main_directory + '/blasting.log', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    if args.verbose:
        logging.getLogger().addHandler(logging.StreamHandler())
    logging.info("------ Blasting module started ------")
    if args.rerun:
        rerun_blast_selection(args.main_directory, args.rerun, args.identity, args.difference, args.e_val,
                              args.coverage, args.bit_score)
    elif args.unique:
        if args.name:
            run_unique(args)
        else:
            print("Optional argument --name (-n) becomes necessary if you do a unique run")
    else:
        run(args)


if __name__ == "__main__":
    main()
