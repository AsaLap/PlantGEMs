import cobra
import copy
import os
import re
import subprocess
import time
import utils


class Blasting:
    name = None
    model = None
    draft = None
    model_fasta_path = None
    model_fasta = None
    model_directory = None
    subject_fasta_path = None
    subject_fasta = None
    subject_directory = None
    blast_result = {}
    gene_dictionary = {}

    identity = 50
    difference = 30
    e_val = 1e-100
    coverage = 20
    bit_score = 300
    """
    identity (int) -- the treshold value of identity to select the subject genes.
    difference (int) -- the percentage of length difference tolerated between subject and query.
    e_val (int) -- the minimum E-Value chosen.
    coverage (int) -- the minimum percentage of coverage of the match.
    bit_score (int) -- the minimum Bit-Score chosen.
    """

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

    # def model_fasta_cut(self):
    #     for seq in self.model_fasta.split(">"):
    #         try:
    #             gene_name = re.search('\w+(\.\w+)*(\-\w+)*', seq).group(0)
    #             self.model_genes[gene_name] = ">" + seq
    #         except AttributeError:
    #             print("Gene name not found in :", seq)
    #             pass

    def blast_run(self):
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
                    gene_name = re.search('\w+(\.\w+)*(\-\w+)*', seq).group(0)
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
                self.blast_result[gene.id] = subprocess.run(blast_request, capture_output=True).stdout.decode(
                    'ascii').split(
                    "\n")[:-1]
            subprocess.run(["rm", "-rf", tmp_dir])
            print("Blast done !\nTotal time : %f s" % (time.time() - total_time))

    def select_genes(self):
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

    def drafting(self):
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


if __name__ == '__main__':
    # test = Blasting("Test", "/home/asa/INRAE/Tests/aracyc.sbml", "/home/asa/INRAE/Tests/query.fasta",
    #                 "/home/asa/INRAE/Tests/subject.fasta")
    # test.blast_run()
    # utils.save_obj(test, "/home/asa/INRAE/Tests/test")
    test = utils.load_obj("/home/asa/INRAE/Tests/grapeTestGeneSelection")
    test.select_genes()
    test.drafting()
    print(len(test.draft.reactions))
    print(len(test.draft.genes))
