import cobra
import re
import time


class Blasting:
    model = None
    model_fasta = None
    subject_fasta = None
    model_genes = {}
    blast_res = {}

    def __init__(self, _model, _model_fasta, _subject_fasta):
        """
        ARGS:
            _model -- the path to the SBML file containing the model for the reconstruction.
            _model_fasta -- the path to the fasta file of the model.
            _subject_fasta -- the path to the fasta file of the subject.
        """
        self.model = cobra.io.read_sbml_model(_model)
        self.model_fasta = open(_model_fasta).read()
        self.subject_fasta = open(_subject_fasta).read()

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
        for gene in self.model.genes:
            if i % 10 == 0:
                print("Protein %i out of %i\nTime : %f s" % (i, x, time.time() - lap_time))
                lap_time = time.time()
            i += 1
            requestBlast = [
                "blastp",
                "-subject",
                subjectDir,
                "-query",
                newDir + gene.id + ".fa",
                "-outfmt",
                "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
            blast_res[gene.id] = subprocess.run(requestBlast, capture_output=True).stdout.decode('ascii').split("\n")[
                                 :-1]
        # subprocess.run(["rm", "-rf", newDir])
        # print("Blast done !\nTotal time : %f s" % (time.time() - total_time))