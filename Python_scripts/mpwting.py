# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the preparation of the required files for the 
Pathway Tools software reconstruction and launching of reconstruction 
using mpwt package from AuReMe."""

import module
import mpwt
import multiprocessing
import numpy as np
import re
import subprocess
import sys
import utils


class Mpwting(module.Module):

    def __init__(self, _name, _main_directory, _m_rna, _chromosome_type):
        super().__init__(_name, _main_directory)
        self.m_rna = _m_rna
        self.chromosome_type = _chromosome_type
        self.regions_dict = self._get_sequence_region()
        self.fasta_file_path = self._find_fasta(self.name)
        self.gff_file_path = self._find_gff(self.name)
        self.eggnog_file_path = self._find_eggnog(self.name)

    def _get_sequence_region(self):
        """Function that browses the .gff file and gets the name of each gene,
        the position in the genome and each corresponding region and protein(s).

        PARAMS:
            m_rna (bool) -- tells the function if it must search the mRNA (True) line or CDS (False).
        RETURNS:
            regions_dict -- a dictionary containing all the gathered
            information (see pipelinePT() for the structure).
        """

        regions_dict = {}
        gff_file = utils.read_file_listed(self.gff_file_path)
        protein_found = False  # Boolean to avoid testing a protein on each line that has already been found.
        for line in gff_file:
            if "\tgene\t" in line:  # Searching the gene's information
                protein_found = False
                spl = line.split("\t")
                region = spl[0]
                try:
                    gene = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
                except AttributeError:
                    try:
                        gene = re.search('(?<=ID=)(gene:)*\w+(\.\w+)*(\-\w+)*', line).group(0)
                        if "gene:" in gene:
                            gene = gene[5:]
                    except AttributeError:
                        print("The gene name hasn't been found...")
                        gene = ""
                        pass
                if region not in regions_dict.keys():
                    regions_dict[region] = {}
                if spl[6] == "+":
                    regions_dict[region][gene] = {"Start": spl[3], "End": spl[4], "Proteins": {}}
                else:
                    regions_dict[region][gene] = {"Start": spl[4], "End": spl[3], "Proteins": {}}
            if self.m_rna:  # Searching the protein's information
                if "RNA\t" in line:
                    try:
                        protein = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
                        regions_dict[region][gene]["Proteins"][protein] = []
                        protein_found = True
                    except AttributeError:
                        print("The mRNA has no attribute 'Name='...")
                        regions_dict[region][gene]["Proteins"]["None"] = []
            else:
                if not protein_found and "CDS\t" in line:  # In case the gff file needs to be looked at on the CDS and
                    # not the mRNA to corresponds to the TSV file
                    try:  # Searching for CDS's ID instead of m_rna.
                        protein = re.search('(?<=ID=)[CcDdSs]*[:-]*\w+(\.\w+)*', line).group(0)[4:]
                        regions_dict[region][gene]["Proteins"][protein] = []
                        protein_found = True
                    except AttributeError:
                        print("The CDS has no attribute 'ID='...")
                        regions_dict[region][gene]["Proteins"]["None"] = []
            if "\tCDS\t" in line:  # Searching the exon's information
                spl = line.split("\t")
                regions_dict[region][gene]["Proteins"][protein].append([int(spl[3]), int(spl[4])])
        return regions_dict

    def _make_protein_correspondence_file(self):
        """Function to create a csv file with the correspondence between a protein and the associated gene.

        PARAMS:
            wd (str) -- the path of the working directory to save the file.
            name (str) -- the name of the model being treated.
            regions_dict -- the dictionary that contains the information needed (see pipelinePT() for the structure).
        """
        correspondence = []
        for region in self.regions_dict.keys():
            for gene in self.regions_dict[region].keys():
                for protein in self.regions_dict[region][gene]["Proteins"].keys():
                    correspondence.append([gene.upper(), protein.upper()])
        utils.write_csv(self.main_directory, correspondence, "protein_corres_" + self.name, "\t")

    def _make_dat_files(self):

        print("\nWARNING ! :\n - If there are circular chromosomes in your data, you have to manually",
              "correct the field 'CIRCULAR?' in the .dat file by changing 'N' (no) with 'Y' (yes).\n")
        circular = 'N'
        dat_file_str_list = []
        if self.chromosome_type == "NONE":
            for i in self.regions_dict.keys():
                dat_file_str_list.append('ID\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                                         % (i, circular, self.main_directory + i + '.pf',
                                            self.main_directory + i + '.fsa'))
        elif self.chromosome_type == ":CONTIG":
            for i in self.regions_dict.keys():
                dat_file_str_list.append('ID\t%s\nchromosome_type\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                                         % (i, self.chromosome_type, self.main_directory + i + '.pf',
                                            self.main_directory + i + '.fsa'))
        else:
            for i in self.regions_dict.keys():
                dat_file_str_list.append(
                    'ID\t%s\nchromosome_type\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n' % (
                        i, self.chromosome_type, circular, self.main_directory + i + '.pf',
                        self.main_directory + i + '.fsa'))
        utils.write_file(self.main_directory, "genetic-elements" + ".dat", dat_file_str_list)

    def _make_fsa_files(self):
        with open(self.fasta_file_path, "r") as file:
            fasta = file.read()
        fasta = fasta.split(">")
        fasta = list(filter(None, fasta))
        regions_list = list(self.regions_dict.keys())
        for i in fasta:
            region = re.search("\w+(\.\w+)*(\-\w+)*", i).group(0)
            if region in regions_list:
                regions_list.remove(region)
                utils.write_file(self.main_directory, region + ".fsa", i)

    def _make_pf_files(self):
        tsv = utils.read_file_listed(self.eggnog_file_path)
        list_index = list(np.arange(0, len(tsv)))
        for region in self.regions_dict.keys():
            sub_pf = []
            for gene in self.regions_dict[region].keys():
                for protein in self.regions_dict[region][gene]["Proteins"].keys():
                    found = False
                    for i in list_index:
                        if found:
                            break
                        if protein in tsv[i]:
                            list_index.remove(i)
                            found = True
                            sub_pf.append(self._eggnog_file_parser(region, gene, protein, tsv[i]))
            if sub_pf:
                f = open(self.main_directory + region + ".pf", "w")
                for i in sub_pf:
                    for j in i:
                        f.write(j)
                f.close()

    def _eggnog_file_parser(self, region, gene, protein, line):
        """Sub-function of make_pf_files() to write the info in the correct order for each protein.

        PARAMS:
            region (str) -- actual genomic scope.
            gene (str) -- one gene id in the scope.
            protein (str) -- one of the proteins coded by the gene above.
            start (int) -- the start position of the sequence.
            line (str) -- the line corresponding to the protein in the .tsv file.
        RETURNS:
            info (str) -- a string with all the information and with the correct
            file architecture settings for the .pf file.
        """

        start = self.regions_dict[region][gene]["Start"]
        end = self.regions_dict[region][gene]["End"]
        exon_pos = self.regions_dict[region][gene]["Proteins"][protein]

        info = []
        spl = line.split("\t")
        info.append("ID\t" + gene + "\n")
        if spl[5]:
            info.append("NAME\t" + spl[5] + "\n")
        else:
            info.append("NAME\tORF\n")
        info.append("STARTBASE\t" + start + "\n")
        info.append("ENDBASE\t" + end + "\n")
        spl[21] = spl[21].replace("\n", "")
        if spl[21]:
            info.append("FUNCTION\t" + spl[21] + "\n")
        else:
            info.append("FUNCTION\tORF\n")
        info.append("PRODUCT-chromosome_type\tP\n")
        if spl[7]:
            for res in spl[7].split(","):
                info.append("EC\t" + res + "\n")
        for exon in exon_pos:
            info.append("CODING-SEGMENT\t" + str(exon[0]) + "-" + str(exon[1]) + "\n")
        if spl[6]:
            go = spl[6].split(",")
            for i in go:
                info.append("DBLINK\t" + i + "\n")
        info.append("//\n")
        return info

    def pipeline(self):
        self._make_protein_correspondence_file()
        self._make_dat_files()
        self._make_fsa_files()
        self._make_pf_files()


def make_taxon_file(wd, taxon_name_list):
    """Function to make the taxon_id.tsv file.

    PARAMS:
        taxon_name_list (list) -- the list containing the name of the organism and its taxon id.
    """

    res = "species\ttaxon_id\n"
    for i in taxon_name_list:
        res += i[0] + "\t" + str(i[1]) + "\n"
    utils.write_file(wd, "taxon_id.tsv", res)

# def make_organism_parameters(wd, species, abbrev, rank, storage="file", private="NIL", tax=2, codon=1, mito_codon=1):
#     # Choose tax = 1(4) for Bacteria, 2(5) for Eukaryota and 3(6) for Archae (2 is default).
#     taxonomy_dict = {1: "TAX-2", 2: "TAX-2759", 3: "TAX-2157",
#                      4: "2", 5: "2759", 6: "2157"}
#     info = []
#     # Making the random ID
#     random_id = random.choice(string.ascii_lowercase)
#     string_choice = string.ascii_lowercase + "0123456789"
#     for loop in range(random.randint(1, 10)):
#         random_id += random.choice(string_choice)
#     info.append("ID\t" + random_id + "\n")
#     info.append("STORAGE\t" + storage + "\n")
#     info.append("NAME\t" + species + "\n")
#     info.append("ABBREV-NAME\t" + abbrev + "\n")
#     info.append("PRIVATE?\t" + private + "\n")
#     info.append("RANK\t" + str(rank) + "\n")
#     info.append("ORG-COUNTER\t\n")  # TODO : Test with or without it
#     info.append("DOMAIN\t" + taxonomy_dict[tax] + "\n")
#     info.append("CODON-TABLE\t" + str(codon) + "\n")
#     info.append("MITO-CODON-TABLE\t" + str(mito_codon) + "\n")
#     info.append("DBNAME\t" + abbrev + "DBcyc\n")
#     info.append("NCBI-TAXON-ID\t" + taxonomy_dict[tax + 3] + "\n")
#     utils.write_file(wd, "organism-params.dat", info)


def pipeline(data):  # TODO : make this function fit with the new architecture
    """The function to make all the pipeline working, from creation of 
    the files to Pathway Tools via mpwt.
    
    PARAMS:
        data (str) -- the path to a txt file containing the path
        to each ini file for each organism to run.
    -- Parameters taken from an ini file:
        gff_filename (str) -- the name of the .gff file for the organism.
        fasta_filename (str) -- the name of the .fasta file for the organism.
        eggnog_filename (str) -- the name of the .tsv file from EggNOG for the organism.
        genetic_type (str) -- indication if the sequences of the organism are assembled
        as chromosomes or contigs (or else, see Pathway Tools user guide).
        taxon_id (int) -- the NCBI taxon ID of the species.
        m_rna (bool) -- decides if the function must search the mRNA line (True)
        or CDS line (False) to get the name of the protein.
        name (str) -- the name of the species database.
    
    -- Structure of regions_dict (created in this function):
    {Region name (str):
        {Gene name (str):
            {"Start": int, "End": int, "Proteins":
                {Protein name (str):[[begin, end] the protein(s)'s CDS positions (list of list of int]}
        }
    }
    """

    index = utils.read_file_listed(data)
    main_directory = index.pop(0).rstrip("/ ") + "/"
    input_directory = main_directory + "input/"
    output_directory = main_directory + "output/"
    log_directory = main_directory + "log/"
    utils.make_directory(input_directory)
    utils.make_directory(output_directory)
    utils.make_directory(log_directory)

    taxon_name_list = []
    cpu = 0
    for ini in index:
        if ini:
            cpu += 1
            parameters = utils.read_config(ini.rstrip())
            # gff_filename = parameters["FILES"]["GFF"]
            # fasta_filename = parameters["FILES"]["FASTA"]
            # eggnog_filename = parameters["FILES"]["EGGNOG"]
            chromosome_type = parameters["INFO"]["chromosome_type"]
            taxon_id = int(parameters["INFO"]["NCBI_TAXON_ID"])
            m_rna = parameters.getboolean("INFO", "mRNA")
            species_name = parameters["INFO"]["DATABASE_NAME"]
            species_directory = input_directory + species_name + "/"
            utils.make_directory(species_directory)

            # Keeping some information for later
            taxon_name_list.append([species_name, taxon_id])

            # Creating the organism's instance and the corresponding pipeline to prepare for MPWT
            organism = Mpwting(species_name, main_directory, m_rna, chromosome_type)
            organism.pipeline()

    # Creating the tsv file for the taxon IDs
    make_taxon_file(input_directory, taxon_name_list)  # Quite unspecific, could be improved (genetic_type, codon table)
    print("------\nCreation of the files finished\n------")

    # Counting the number of cpu to use
    nb_cpu = multiprocessing.cpu_count()
    if nb_cpu <= cpu:
        cpu = nb_cpu - 1
    print("Number of CPU used : ", cpu)

    # Starting the mpwt script
    mpwt.multiprocess_pwt(input_folder=input_directory, output_folder=output_directory, patho_inference=True,
                          patho_hole_filler=False, patho_operon_predictor=False, pathway_score=1, dat_extraction=True,
                          number_cpu=nb_cpu, size_reduction=False, patho_log=log_directory, ignore_error=False,
                          taxon_file=True, verbose=True)


if __name__ == "__main__":
    globals()[sys.argv[1]](*sys.argv[2:])
