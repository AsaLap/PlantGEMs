# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Aquitaine
# Mars 2020

import cobra
from os.path import join
import subprocess

class BiModel:
    """This class defines a metabolic model with only genes and reactions"""
    def __init__(self, _name):
        self.name = _name
        self.dictionary = {}
        self.genes = []
        self.reactions = []
    
    def buil_model(self, _path):
        """Build a model with only genes and associated reactions.
        
        ARGS : 
        _path -- the GEM where to extract the genes and reactions.
        RETURN : None
        """
        model = cobra.io.read_sbml_model(_path)
        for i in range(len(model.genes)):
            gene = model.genes[i].name
            try:
                reaction = model.reactions[i].name
            except IndexError:
                print("Reaction not found : ",i)
                reaction = "NA"
                pass
            if gene in self.dictionary.keys():
                self.dictionary[gene].append(reaction)
            else:
                self.dictionary[gene] = [reaction]
        self.genes = self.dictionary.keys()
        self.reactions = self.dictionary.values()
        
def blast_p(model, workDir, subjectFile, queryFile):
    """Runs multiple blastp between the subject and each protein of the query file.
    
    ARGS : 
    model -- a BiModel object.
    subjectPATH -- the subject file for the blastp.
    queryPATH -- the query file for the blastp corresponding to the BiModel object given.
    RETURN : 
    dico_res -- a dictionary containing the result of the blastp command for each gene,
    the genes being the keys.
    """
    newDir = workDir + "Proteins/"
    queryDir = workDir + queryFile
    subjectDir = workDir + subjectFile
    subprocess.run(["mkdir", newDir])
    
    ###Creation of one file per protein of the query###
    fileFasta = open(queryDir)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(newDir+geneName+".fa", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    
    ###Blastp###
    dico_res = {}
    for gene in model.genes:
        requestBlastp = [
            "blastp",
            "-subject",
            subjectDir,
            "-query",
            newDir+"in_"+gene+".fa",
            "-outfmt",
            "10 delim=, qseqid sseqid qlen length slen nident pident score evalue bitscore"]
        dico_res[gene] = subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
    subprocess.run(["rm", "-rf", workDir+newDir])
    
    
if __name__=='__main__':
    ###Files and working directory###
    WD = '/home/asa/INRAE/Work/Drafts/Data/Tests/'
    modelGem = 'AraGEM3.xml'
    modelGemFasta = 'genomic.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    
    ###Pipeline###
    modelGemBiModel = BiModel('AraGEM')
    modelGemBiModel.buil_model(WD + modelGem)
    blast_p(modelGemBiModel, WD, tomatoFasta, modelGemFasta)