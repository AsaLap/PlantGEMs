# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Mars 2020

import cobra
from os.path import join
import subprocess
import pickle
import time
import matplotlib.pyplot as plt


def blast_run(workDir, model, queryFile, subjectFile):
    """Runs multiple blasts between a subject file and each protein of the query file.
    
    ARGS : 
        model -- the GEM of reference.
        workDir -- the working directory where the files are, and where the proteins files
        will be temporarily stored.
        subjectFile -- the subject file for the blast.
        queryFile -- the query file for the blast corresponding to the reference GEM's CDS fasta.
    RETURN : 
        blast_res -- dictionary containing the result of the blastp command for each gene of the model.
    """
    ###Concatenation of string to get the exact paths###
    newDir = workDir + "Proteins_tmp/"
    queryDir = workDir + queryFile
    subjectDir = workDir + subjectFile
    ###Creation of temporary directory containing one file per protein of the query###
    subprocess.run(["mkdir", newDir])
    print("\nCreating the individual CDS fasta files...")
    fileFasta = open(queryDir)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(newDir+geneName+".fa", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    print("...done !")
    ###Blast###
    print("\nLaunching the blast !")
    i = 1
    x = len(model.genes)
    total_time = lap_time = time.time()
    blast_res = {}
    for gene in model.genes:
        if i%10==0:
            print("Protein %i out of %i\nTime : %f s" %(i,x, time.time() - lap_time))
            lap_time = time.time()
        i+=1
        requestBlast = [
            "blastp",
            "-subject",
            subjectDir,
            "-query",
            newDir+"in_"+gene.name+".fa",
            "-outfmt",
            "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
        blast_res[gene.name] = subprocess.run(requestBlast, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
    subprocess.run(["rm", "-rf", newDir])
    print("Blast done !\nTotal time : %f s" %(time.time() - total_time))
    return blast_res


def select_genes(blast_res, treshold = 80):
    """Select the target organism's genes with a good score.
    
    ARGS:
        blast_res -- the dictionary with the results of the blastp.
        treshold -- the treshold value of identity for selection of the target genes.
    RETURN:
        dico_genes -- a dictionary with model gene as key and corresponding target 
        key and coverage value as value.
    """
    dico_genes = {}
    for key in blast_res.keys():
        for res in blast_res[key]:
            if float(res.split(",")[6]) >= treshold:
                try:
                    if dico_genes[key][1] < res.split(",")[6]:
                        dico_genes[key] = [res.split(",")[2], res.split(",")[6]]
                except KeyError:
                    dico_genes[key] = [res.split(",")[2], res.split(",")[6]]
    return dico_genes


def save_obj(obj, path):
    with open(path + '.pkl', 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(path):
    with open(path + '.pkl', 'rb') as f:
        return pickle.load(f)


def pipeline(WD, ref_gem, queryFile, subjectFile):
    model = cobra.io.read_sbml_model(WD + ref_gem)
    blast_res = blast_run(WD, model, queryFile, subjectFile)
    save_obj(blast_res, WD + "resBlastp")
    # blast_res = load_obj(WD + "resBlastp")
    # draft = select_genes(blast_res)


if __name__=='__main__':
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Drafts/Tomato_Arabidopsis/'
    WDkiw = '/home/asa/INRAE/Work/Drafts/Kiwi_Arabidopsis/'
    WDcuc = '/home/asa/INRAE/Work/Drafts/Cucumber_Arabidopsis/'
    WDche = '/home/asa/INRAE/Work/Drafts/Cherry_Arabidopsis/'
    
    WD = '/home/asa/INRAE/Work/Drafts/Tests/'
    
    aragem = 'AraGEM3.xml'
    aragemFasta = 'genomic.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    kiwiFasta = 'Hongyang_pep_v2.0.fa'
    cucumberFasta = 'Gy14_pep_v2.fa'
    cherryFasta = 'PRUAV_Regina_CDS.fa'
    
    ###Main###
    #For the tomato
    pipeline(WDtom, aragem, aragemFasta, tomatoFasta)
    #For the kiwifruit
    # pipeline(WDkiw, aragem, aragemFasta, kiwiFasta)
    #For the cucumber
    # pipeline(WDcuc, aragem, aragemFasta, cucumberFasta)
    #For the cherry
    # pipeline(WDche, aragem, aragemFasta, cherryFasta)