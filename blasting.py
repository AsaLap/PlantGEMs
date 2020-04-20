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
from statistics import mean


def save_obj(obj, path):
    with open(path + '.pkl', 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(path):
    with open(path + '.pkl', 'rb') as f:
        return pickle.load(f)


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


def select_genes(blast_res, treshold = 75, e_val = 1e-100):
    """Select the target organism's genes with a good score.
    
    ARGS:
        blast_res -- the dictionary with the results of the blastp.
        treshold -- the treshold value of identity to select the target genes.
    RETURN:
        dico_genes -- a dictionary with model gene as key and corresponding target 
        key and coverage value as value.
    """
    dico_genes = {}
    for key in blast_res.keys():
        for res in blast_res[key]:
            if float(res.split(",")[6]) >= treshold:
                try:
                    if float(res.split(",")[8]) <= e_val:
                        dico_genes[key].append(res.split(",")[2])
                except KeyError:
                    if float(res.split(",")[8]) <= e_val:
                        dico_genes[key] = [res.split(",")[2]]
    return dico_genes


def drafting(model, dico_genes, model_name):
    new_model = cobra.Model(model_name)
    #Getting the associated reaction genes
    for key in dico_genes.keys():
        string_reaction_rule = "( "
        for i in dico_genes[key]:
            string_reaction_rule += i + " or "
        string_reaction_rule = string_reaction_rule[:-3] + ")"
        #Getting the reactions for each gene and changing gene_reaction_rule to current genes
        for i in model.genes.get_by_id(key).reactions:
            x = i
            x.gene_reaction_rule = string_reaction_rule
            new_model.add_reactions([x])
    return new_model
        


def pipeline(WD, ref_gem, queryFile, subjectFile, modelName):
    model = cobra.io.read_sbml_model(WD + ref_gem)
    # blast_res = blast_run(WD, model, queryFile, subjectFile)
    # save_obj(blast_res, WD + "resBlastp")
    blast_res = load_obj(WD + "resBlastp")
    dico_genes = select_genes(blast_res, 75, 1e-100)
    
    new_model = drafting(model, dico_genes, modelName)
    new_model.add_metabolites(model.metabolites)
    
    #Printing of verifications
    #Test blast_res, search for the genes with no results
    no_results = []
    for key in blast_res.keys():
        if not blast_res[key]:
            no_results.append(key)
    print("The",len(no_results),"genes that have no results : ", no_results)
    
    #
    test = []
    for val in dico_genes.values():
        for i in val:
            test.append(i)
    print("Nb of genes in ref model : %s\nNb of genes in the new model : %s\nNb of genes in dico_genes : %s\nNb of values in dico_genes : %s (without doublons : %s)" %(len(model.genes), len(new_model.genes), len(dico_genes.keys()), len(test), len(set(test))))
    print("----------------------------------------")
    return test, new_model
    


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
    tomatoDraft = pipeline(WDtom, aragem, aragemFasta, tomatoFasta, "Tomato")
    #For the kiwifruit
    pipeline(WDkiw, aragem, aragemFasta, kiwiFasta, "Kiwi")
    #For the cucumber
    pipeline(WDcuc, aragem, aragemFasta, cucumberFasta, "Cucumber")
    #For the cherry
    pipeline(WDche, aragem, aragemFasta, cherryFasta, "Cucumber")