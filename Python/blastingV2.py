import cobra
import os
import re
import subprocess
import time


class Blasting:
    model = None
    model_fasta_path = None
    model_fasta = None
    model_dir = None
    subject_fasta_path = None
    subject_fasta = None
    subject_dir = None
    model_genes = {}
    blast_res = {}

    def __init__(self, _model, _model_fasta_path, _subject_fasta_path):
        """
        ARGS:
            _model -- the path to the SBML file containing the model for the reconstruction.
            _model_fasta_path -- the path to the fasta file of the model.
            _subject_fasta_path -- the path to the fasta file of the subject.
        """
        self.model = cobra.io.read_sbml_model(_model)
        self.model_fasta_path = _model_fasta_path
        self.model_dir = os.path.split(self.model_fasta_path)[0]
        self.model_fasta = open(self.model_fasta_path).read()
        self.subject_fasta_path = _subject_fasta_path
        self.subject_fasta = open(self.subject_fasta_path).read()
        self.subject_dir = os.path.split(self.subject_fasta_path)[0]

    def model_fasta_cut(self):
        for seq in self.model_fasta.split(">"):
            try:
                gene_name = re.search('\w+(\.\w+)*(\-\w+)*', seq).group(0)
                self.model_genes[gene_name] = ">" + seq
            except AttributeError:
                print("Gene name not found in :", seq)
                pass

    def blast_run(self):
        print("\nLaunching the blast !")
        i, x = 1, len(self.model.genes)
        total_time = lap_time = time.time()
        tmp_dir = self.model_dir + "/tmp_dir/"
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
            self.blast_res[gene.id] = subprocess.run(blast_request, capture_output=True).stdout.decode('ascii').split(
                "\n")[:-1]
        subprocess.run(["rm", "-rf", tmp_dir])
        print("Blast done !\nTotal time : %f s" % (time.time() - total_time))


if __name__ == '__main__':
    test = Blasting("/home/asa/INRAE/Tests/aracyc.sbml", "/home/asa/INRAE/Tests/query.fasta", "/home/asa/INRAE/Tests/subject.fasta")
    test.blast_run()