# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Aquitaine
# Mars 2020

import cobra
from os.path import join
import subprocess
import pickle


class BiModel:
    """This class defines a metabolic model with only genes and reactions"""
    def __init__(self, _name):
        self.name = _name
        self.dictionary = {}
        self.dico_res = {}
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
        workDir -- the working directory where the files are, and where the proteins files
        will be temporarily stored.
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
    print("\nCreating the small fasta files...")
    fileFasta = open(queryDir)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(newDir+geneName+".fa", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    print("...done !")
    
    ###Blastp###
    print("\nLaunching the blastp !")
    i = 1
    x = len(model.genes)
    print(model.genes)
    # for gene in model.genes:
    for gene in model.genes[0:2]:
        if i%50==0:
            print("Protein %i out of %i" %(i,x))
        i+=1
        requestBlastp = [
            "blastp",
            "-subject",
            subjectDir,
            "-query",
            newDir+"in_"+gene+".fa",
            "-outfmt",
            "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
        model.dico_res[gene] = subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
    subprocess.run(["rm", "-rf", newDir])
    print(model.dico_res)

def save_model(workingDir, model):
    with open(workingDir+'model_data.pkl', 'wb') as output:
        pickle.dump(model, output, pickle.HIGHEST_PROTOCOL)
    del model
    
def read_model(modelDir):
    with open(modelDir, 'wb') as input:
        model = pickle.load(input)
    return model


if __name__=='__main__':
    ###Files and working directory###
    WD = '/home/asa/INRAE/Work/Drafts/Data/Tests/'
    modelGem = 'AraGEM3.xml'
    modelGemFasta = 'genomic.in.fasta'
    modelGemFastaTest = 'testBlast.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    tomatoTestFasta = 'TomatoTest.fasta'
    
    ###Pipeline###
    modelGemBiModel = BiModel('AraGEM')
    modelGemBiModel.buil_model(WD + modelGem)
    # blast_p(modelGemBiModel, WD, tomatoFasta, modelGemFasta)
    # save_model(WD, modelGemBiModel)
    