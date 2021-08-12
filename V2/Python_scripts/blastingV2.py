import cobra
import copy
import os
import re
import subprocess
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
        self.__identity = 50
        self.__difference = 30
        self.__e_val = 1e-100
        self.__coverage = 20
        self.__bit_score = 300
        """
            identity (int) -- the treshold value of identity to select the subject genes.
            difference (int) -- the percentage of length difference tolerated between subject and query.
            e_val (int) -- the minimum E-Value chosen.
            coverage (int) -- the minimum percentage of coverage of the match.
            bit_score (int) -- the minimum Bit-Score chosen.
        """
        # TODO : Bit-score and E-value are mathematically related, investigate the use of both.
        # TODO : Make a rapid explanation of the utility of all those terms.

    @property
    def identity(self):
        print("Getting identity value...")
        return self.__identity

    @identity.setter
    def identity(self, value):
        if value in range(0, 101):
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
        if value in range(0, 101):
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
        if value in range(0, 11):
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
        if value in range(0, 101):
            print("Setting coverage value to %s" % (str(value)))
            self.__coverage = value
        else:
            print("Denied : value must be between 0 and 100 (both included)")

    @property
    def bit_score(self):
        print("Getting Bit-score value...")
        return self.__coverage

    @bit_score.setter
    def bit_score(self, value):
        if value in range(0, 10001):
            print("Setting Bit-score value to %s" % (str(value)))
            self.__bit_score = value
        else:
            print("Denied : value must be between 0 and 10000 (both included)")

        # def model_fasta_cut(self):
        """ Function put aside because the blastp command can only be lauched on files and not strings in a program so 
        every gene is written on a file in the 'blast_run' method.
        Otherwise this function was to store all the genes' strings into an instance attribute named 'model_genes'"""

    #     for seq in self.model_fasta.split(">"):
    #         try:
    #             gene_name = re.search('\w+(\.\w+)*(\-\w+)*', seq).group(0)
    #             self.model_genes[gene_name] = ">" + seq
    #         except AttributeError:
    #             print("Gene name not found in :", seq)
    #             pass

    def _blast_run(self):
        """Runs multiple blasts between the model and the subject."""

        if not self.gene_dictionary:
            print("\nLaunching the blast !")
            i, x = 1, len(self.model.genes)
            total_time = lap_time = time.time()
            tmp_dir = self.model_directory + "/tmp_dir/"
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
            subprocess.run(["rm", "-rf", tmp_dir])
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
        for reac in self.model.reactions:
            to_add = []
            to_search = reac.gene_reaction_rule.split(" or ")
            for gene in to_search:
                try:
                    to_add += self.gene_dictionary[gene]
                except KeyError:
                    pass
            ###TODO : change the proteins's ID in to_add in corresponding genes
            string_reaction_rule = " or ".join(to_add)
            if string_reaction_rule:
                x = copy.deepcopy(reac)
                x.gene_reaction_rule = string_reaction_rule
                self.draft.add_reactions([x])

    def build(self):
        history_dir = os.path.dirname(self.subject_directory) + "/history/" + self.name + "/"
        self._blast_run()
        utils.save_obj(self, history_dir + "blasted")
        self._select_genes()
        utils.save_obj(self, history_dir + "gene_selected")
        self._drafting()
        utils.save_obj(self, history_dir + "drafted")


if __name__ == '__main__':
    # Tests for portable PC
    # test = Blasting("Test", "/home/asa/INRAE/Tests/aracyc.sbml", "/home/asa/INRAE/Tests/query.fasta",
    #                 "/home/asa/INRAE/Tests/subject.fasta")
    # test.blast_run()
    # utils.save_obj(test, "/home/asa/INRAE/Tests/test")
    # test = utils.load_obj("/home/asa/INRAE/Tests/grapeTestGeneSelection")

    # Tests for home computer
    test = Blasting("Test", "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/aracyc.sbml",
                    "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/aracyc.fasta",
                    "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/tomato.fasta")
    test.e_val = 1
    print(test.e_val)
    test.identity = 10
    print(test.identity)
    test.coverage = 80
    print(test.coverage)
    test.bit_score = 700
    print(test.bit_score)
    test.difference = 2
    print(test.difference)

    # test.build()
    # utils.save_obj(test, "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/test")
    # test = utils.load_obj("/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/tomatoBlasting")
    # print(len(test.model.reactions))
    # print(len(test.model.metabolites))
    # print(len(test.model.genes))

    #
    # Common orders
    # test.select_genes()
    # test.drafting()
    # print(len(test.draft.reactions))
    # print(len(test.draft.genes))

    # test = Blasting("Test", "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/aracyc.sbml",
    #                 "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/aracyc.fasta",
    #                 "/home/asa/INRAE/StageMaster_2020/Work/Tests/Tomato_Aracyc/tomato.fasta")
    # utils.save_obj(test, "/home/asa/INRAE/StageMaster_2020/Work/Tests/testpickle")
    # test = utils.load_obj("/home/asa/INRAE/StageMaster_2020/Work/Tests/testpickle")
