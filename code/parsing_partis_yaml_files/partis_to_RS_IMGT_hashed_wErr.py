#!/bin/env python3

# neccesary modules to import 
import yaml
import yamlloader
import pandas as pd
import argparse
from Bio.Seq import Seq
from pathlib import Path
import re
import partis_CDR3_rf as RF
import sys
import logging

# these are the initial columns that are returned from the first part of partis parsing 
out_col_names = ["PRODUCTIVE", "PARTIAL_V", "PARTIAL_J", "VGene", "ProbVAllele", "ProbVGene", "ProbVFam", "DGene", "ProbDAllele", "ProbDGene", "ProbDFam", "JGene", "ProbJAllele", "ProbJFam", \
        "RV", "RD1", "RD2", "RJ", "CDR3", "CDR3Length", "CDR3ReadFrame", \
        "nVDNs", "nDJNs", "AA@Cys2", "AA@JInv", "DGene_RF"]

# this are the columns that are removed after the information comes out of the SMUA 
out_col_names_RS = ["VGene", "ProbVAllele", "ProbVGene", "ProbVFam", "DGene", "ProbDAllele", "ProbDGene", "ProbDFam", "JGene", "ProbJAllele", "ProbJFam", \
        "RV", "RD1", "RD2", "RJ", "CDR3", "CDR3Length", "CDR3ReadFrame", \
        "nVDNs", "nDJNs", "AA@Cys1", "AA@Cys2", "AA@JInv", "#VBases", "#VSubstitutions", "#VInsertions", "#VDeletions", "MuFreq", "MDL"]

# the are the column names that going in the partis evaluation information 
out_col_names_extra = ["VALID", "SMUA", "FUNCTIONAL", "PRODUCTIVE", "PARTIAL_V", "PARTIAL_J"]

# dictionaries to convery between partis and cloanalyst list output
convert_dict_genes = {"v_gene": "VGene", "d_gene": "DGene", "j_gene": "JGene"}
convert_dict_support = {"v_gene": ["Allele", "Gene", "Fam"], "d_gene": ["Allele", "Gene", "Fam"], "j_gene": ["Allele", "Fam"]}
convert_dict_recomb_p1 = {"v_3p_del": "RV", "d_5p_del": "RD1", "d_3p_del": "RD2", "j_5p_del": "RJ"}
convert_dict_recomb_p2 = {"vd_insertion": "nVDNs", "dj_insertion": "nDJNs"}

# dictionary with the IMGT numbers 
imgt_dict = {"heavy": {"cdr1_start": 27, "fr2_start": 39, "cdr2_start": 56, "fr3_start":66, "cdr3_start": 105}, \
        "kappa": {"cdr1_start": 27, "fr2_start": 39, "cdr2_start": 56, "fr3_start": 66, "cdr3_start": 105}, \
        "lambda": {"cdr1_start": 27, "fr2_start": 39, "cdr2_start":  56, "fr3_start": 66, "cdr3_start": 105}}

imgt_dict_cys = {"heavy": 23, "kappa": 23, "lambda": 23}

# the dictionary for the outfile names for the chains 
out_names = {"heavy": "VH", "kappa": "VK", "lambda": "VL"}

imgt_hash_path_shared = "/datacommons/dhvi/partis_germlines/IMGT_V_region_annotation_HASHES/"

imgt_hash_paths = {"human": {"heavy": imgt_hash_path_shared + "partis_human_IGHV.yaml" ,\
        "kappa": imgt_hash_path_shared + "partis_human_IGKV.yaml", "lambda": imgt_hash_path_shared + "partis_human_IGLV.yaml"},\
        "mouse":{"heavy": imgt_hash_path_shared + "partis_mouse_IGHV.yaml",\
        "kappa": imgt_hash_path_shared + "partis_mouse_IGKV.yaml", \
        "lambda": imgt_hash_path_shared + "partis_mouse_IGLV.yaml"}, \
        "rhesus": {"heavy": imgt_hash_path_shared + "KimDB_rhesus_IGHV_noAlleles.yaml", "kappa":"", "lambda":""}}


# function that counts the differences between the observe and naive sequences but ignores missing nucleotides in the observed sequence  
def find_diffs(naive_aligned_V, obs_aligned_V):

    diffs = [1 if (naive_aligned_V[c]!=i and i!="N" and i!="-" and naive_aligned_V[c]!="N" and naive_aligned_V[c]!="-") else 0 for c, i in enumerate(obs_aligned_V)]
    return sum(diffs)

# function for the writing the SMUA file 
def SMUA_write(in_dict, fasta_file_path, chain):

    failed_SMUAs = 0

    # lists that will store the information that is aquired by writing the SMUA 
    cys_out = []
    SMUA_flag_list = []
    V_ins_out = []
    V_dels_out = []
    V_subs_out = []
    V_bases_out = []
    V_MuFreq_out = []

    # open the writing file path 
    with open(fasta_file_path, "w") as S:
        
        # for each key and value pair in the dictionary 
        for k, v in in_dict.items():
            
            # first check if its a valid sequence, if its not a valid sequence then there won't be an SMUA 
            if v["VALID"]:
                
                # if it is a valid sequence then there is going to be an SMUA, unless it gets unset later 
                SMUA_flag = True

                # get the SMUA dictionary from the dictionary associated with each UID 
                SMUA_dict = v["SMUA_dict"]
            
                # make sure that there is something in the SMUA dict for the V part of the sequence. This should happen because we get info for the SMUA dict using naive seqs, 
                #   but if there is only one VGene and that sequence is weird then it won't work 
                if [True if "fw1_start" in v["SMUA_dict"].keys() else False][0]:

                    # create the V part of the markeup string using the SMUA dictionary information 
                    string = SMUA_dict["fw1_len"] * "1" + SMUA_dict["CDR1_len"]* "A" + \
                            SMUA_dict["fw2_len"] * "2" + SMUA_dict["CDR2_len"]* "B" + SMUA_dict["fw3_len"] * "3"

                    # add the CDR3 and FW4 part of the SMUA 
                    string = string + SMUA_dict["end"]

                    # now we check if there are any gapps in the naive alignment so that we can add gaps to the SMUA (i.e. there were insertions in the sequence)
                    if re.search("-", v["fasta_dict"]["naive_aligned"]):

                        # keep track of the insertions 
                        counter_adjust = 0

                        # create a list for the gapped sequence 
                        gapped_list = []
    
                        # do an enumerate loop to go through each nt in the naive aligned sequence and know the nucleotide 
                        for c,s in enumerate(v["fasta_dict"]["naive_aligned"]):
                            
                            # if the nt is a gap then we add the dash to the gapped list and update the counter 
                            if s == "-":
                                counter_adjust += 1
                                gapped_list.append("-")

                            else:
                                # in certain instances of weird CDR3 rearrangements sometimes the lengths won't match
                                try:
                                    gapped_list.append(string[c - counter_adjust])
                                except:
                                    SMUA_flag = False
                                    break 
                        
                        # create string from the gapped list
                        gapped_string = "".join(gapped_list)

                    else: 
                        gapped_string = string


                    # check that the gapped string is the same length as the sequences
                    if (len(gapped_string) != len(v["fasta_dict"]["obs_aligned"])) | (len(gapped_string) != len(v["fasta_dict"]["naive_aligned"])):

                        if SMUA_flag != False:
                            SMUA_flag = False

                    # find the number of start "-"
                    num_start_V = len(re.match('^(N*)', v["fasta_dict"]["obs_aligned"]).group(0))
                    gapped_string = "-"*num_start_V + gapped_string[num_start_V:]

                    # find the # of missing bases and adjust so that the reading frame is correct! 
                    starting_gaps = re.search("^N+", v["fasta_dict"]["obs_aligned"])

                    if starting_gaps:
                        num_gaps = starting_gaps.span()[1]
                        check_mod = num_gaps %3

                        if check_mod:

                            num_gaps = num_gaps + (3 - check_mod)
                            gapped_string = "U"*num_gaps + gapped_string[num_gaps:]
                            
                            # update the obs_aligned so that it is mod3
                            v["fasta_dict"]["obs_aligned"] = "N"*num_gaps + v["fasta_dict"]["obs_aligned"][num_gaps:]
                            
                    # find the V gene

                    if re.search("3+", gapped_string) and SMUA_flag==True:
                            
                        obs_aligned_V = v["fasta_dict"]["obs_aligned"][:re.search("3+", gapped_string).span()[1]]
                        naive_aligned_V = v["fasta_dict"]["naive_aligned"][:re.search("3+", gapped_string).span()[1]]

                        V_subs = find_diffs(naive_aligned_V, obs_aligned_V)

                        V_dels = len(re.findall("-", obs_aligned_V))
                        V_ins = len(re.findall("-", naive_aligned_V))

                        V_short = len(re.findall("N", obs_aligned_V))
                        V_bases = len(obs_aligned_V) - V_dels - V_short

                        V_bases_div = max(V_bases, 1)
                        V_MuFreq = round(V_subs/V_bases_div, 10)
                            
                    else:

                        V_dels = "NaN"
                        V_ins = "NaN"
                        V_bases = "NaN"
                        V_MuFreq = "NaN"
                        V_subs = "NaN" 
                                                                            
                else:
                    print("does this actually happen?!?!")
                    SMUA_flag = False
                    V_dels = "NaN"
                    V_ins = "NaN"
                    V_bases = "NaN"
                    V_MuFreq = "NaN"
                    V_subs = "NaN"
                    
                V_ins_out.append(V_ins)
                V_dels_out.append(V_dels)
                V_bases_out.append(V_bases)
                V_subs_out.append(V_subs)
                V_MuFreq_out.append(V_MuFreq)

                SMUA_flag_list.append(SMUA_flag)
                
                if SMUA_flag:

                    # write the aligned input sequence 
                    S.write(">" + k + "\n")
                    #S.write(v["fasta_dict"]["obs_aligned"].replace("N", "-") + "\n")
                    # Your original string from v["fasta_dict"]["obs_aligned"]
                    original_string = v["fasta_dict"]["obs_aligned"]

                    # Replace starting N's
                    modified_string = re.sub(r'^N+', lambda m: '-' * len(m.group()), original_string)

                    # Replace ending N's
                    modified_string = re.sub(r'N+$', lambda m: '-' * len(m.group()), modified_string)

                    # Write to file
                    S.write(modified_string + "\n")

                    # write the aligned UCA sequence
                    S.write(">" + k + "|UCA" + "\n")
                    S.write(v["fasta_dict"]["naive_aligned"] + "\n")

                   # write the mark up string line
                    S.write(">" + k + "|" + v["VGene"] + "|" + v["DGene"] + "|" + v["JGene"] + "\n")
                    S.write(gapped_string + "\n")

                # now we need to find the Cysteine
                    if re.search("A+", gapped_string):
                        cdr1_start = re.search("A+", gapped_string).span()[0]
                        codon_check = []
                        cys_ind = cdr1_start - 3 * (imgt_dict[chain]["cdr1_start"] - imgt_dict_cys[chain])

                        while len(codon_check)<3:
                            nt = v["fasta_dict"]["obs_aligned"][cys_ind]
                            if nt != "-":
                                codon_check.append(nt)
                            cys_ind += 1
                    
                        codon = "".join(codon_check)
                        
                        cys1 = Seq(codon[:3]).translate()
                        
                        if cys1:
                            cys_out.append(str(cys1))
                        else:
                            cys_out.append("0")
                    else:
                        cys_out.append("0")
        
                else:
                    logging.error("SMUA cannot be created for this ID even though partis found a valid rearrangement: " + k)
                    cys_out.append("0")
            else:
                SMUA_flag_list.append(False)
                cys_out.append("0")
                V_ins_out.append("0")
                V_dels_out.append("0")
                V_bases_out.append("0")
                V_MuFreq_out.append("0")
                V_subs_out.append("0")

    return cys_out, SMUA_flag_list, V_ins_out, V_dels_out, V_bases_out, V_subs_out, V_MuFreq_out


def find_fxnl(dataframe_in, chain):

    if chain == "heavy":
        fxnl_col = [True if (r["CDR3ReadFrame"]==0 and r["AA@Cys1"]=="C" and r["AA@Cys2"]=="C" and r["AA@JInv"] == "W" and r["#VInsertions"]%3==0 and r["#VDeletions"]%3 == 0) else False for i,r in dataframe_in.iterrows()]
    else:
        fxnl_col = [True if (r["CDR3ReadFrame"]==0 and r["AA@Cys1"]=="C" and r["AA@Cys2"]=="C" and r["AA@JInv"] == "F" and r["#VInsertions"]%3==0 and r["#VDeletions"]%3 == 0) else False for i,r in dataframe_in.iterrows()]

    dataframe_in["FUNCTIONAL"] = fxnl_col

    return dataframe_in

def get_anarci_info(aa_load, SMUA_dict, chain):
    
    cdr3_start = imgt_dict[chain]["cdr3_start"]


    # find framework 1 start
    if aa_load.loc[0, "1"] == "-":
        
        fw1_start_aa_ANARCI = re.search("-*", "".join(list(aa_load.loc[:, "1": str(imgt_dict[chain]["cdr1_start"])].iloc[0,:]))).span()[1]
        fw1_start_nt = 3 * fw1_start_aa_ANARCI
        gaps_fw1 = len(re.findall("-", "".join(list(aa_load.loc[:, str(fw1_start_aa_ANARCI+1): str(imgt_dict[chain]["cdr1_start"])].iloc[0,:]))))
        SMUA_dict["fw1_start"] = fw1_start_nt
        SMUA_dict["fw1_len"] = (imgt_dict[chain]["cdr1_start"] - 1 - gaps_fw1 - fw1_start_aa_ANARCI) * 3

    
    elif aa_load.loc[0, "1"] == "X":

        fw1_start_aa_ANARCI = re.search("X*", "".join(list(aa_load.loc[:, "1": str(imgt_dict[chain]["cdr1_start"])].iloc[0,:]))).span()[1]
        fw1_start_nt = 3 * fw1_start_aa_ANARCI
        gaps_fw1 = len(re.findall("-", "".join(list(aa_load.loc[:, str(fw1_start_aa_ANARCI+1): str(imgt_dict[chain]["cdr1_start"])].iloc[0,:]))))
        SMUA_dict["fw1_len"] = (imgt_dict[chain]["cdr1_start"] - 1 - gaps_fw1 - fw1_start_aa_ANARCI) * 3
        SMUA_dict["fw1_start"] = 0

    else:

        fw1_start_nt = 0
        fw1_start_aa_ANARCI = 0
        gaps_fw1 = len(re.findall("-", "".join(list(aa_load.loc[:, str(fw1_start_aa_ANARCI+1): str(imgt_dict[chain]["cdr1_start"])].iloc[0,:]))))
        SMUA_dict["fw1_start"] = fw1_start_nt
        SMUA_dict["fw1_len"] = (imgt_dict[chain]["cdr1_start"] - 1 - gaps_fw1 - fw1_start_aa_ANARCI) * 3
        

        # find framework 2 start
    len_CDR1 =  imgt_dict[chain]["fr2_start"] - imgt_dict[chain]["cdr1_start"] - len(re.findall("-", str(list(aa_load.loc[:, str(imgt_dict[chain]["cdr1_start"]) : str(imgt_dict[chain]["fr2_start"])].iloc[0,:]))))
    SMUA_dict["CDR1_len"] = 3* len_CDR1

        # now find the CDR2 start
    fr2_start = imgt_dict[chain]["fr2_start"]
    cdr2_start = imgt_dict[chain]["cdr2_start"]

    gaps_fw2 = len(re.findall("-", "".join(list(aa_load.loc[:, str(fr2_start): str(cdr2_start)].iloc[0,:]))))
    SMUA_dict["fw2_len"] = 3* (cdr2_start - fr2_start  - gaps_fw2)

        # find framework 3 start
    fr3_start = imgt_dict[chain]["fr3_start"]
    len_CDR2 = fr3_start - cdr2_start - len(re.findall("-", str(list(aa_load.loc[:, str(cdr2_start): str(fr3_start)].iloc[0,:]))))
    SMUA_dict["CDR2_len"] = len_CDR2 * 3

        # now find the CDR3 start 
    #cdr3_start = imgt_dict[chain]["cdr3_start"]
    gaps_fw3 = len(re.findall("-", "".join(list(aa_load.loc[:, str(fr3_start): str(cdr3_start)].iloc[0,:]))))
    SMUA_dict["fw3_len"] = 3* (cdr3_start - fr3_start - gaps_fw3)



    return SMUA_dict

def valid_partis_p1(ids_in, return_dict, dict_part, chain): 
    
    enter_dict = {}

    # check productive
    if dict_part["stops"][0]:
        enter_dict["PRODUCTIVE"] = False
    else:
        enter_dict["PRODUCTIVE"] = True

    # find genes and support
    for k,v in convert_dict_genes.items():
        enter_dict[v] = dict_part[k]
        
        gene_part = k.split("_")[0]

        for i in convert_dict_support[k]:
            enter_dict["Prob" + gene_part.upper() + i] = dict_part[gene_part.lower() + "_per_gene_support"][dict_part[k]]

    # find recombination info part 1
    for k,v in convert_dict_recomb_p1.items():
        enter_dict[v] = dict_part[k]

    # find recombination info part 2
    for k,v in convert_dict_recomb_p2.items():
        enter_dict[v] = len(dict_part[k])

    # find the CDR3 information
    enter_dict["CDR3Length"] = dict_part["cdr3_length"] - 6
    enter_dict["CDR3"] = [i[3:-3] for i in dict_part["cdr3_seqs"]][0]
    
    # find the reading frame
    if chain == "heavy":
        enter_dict["DGene_RF"] = RF.reading_frame(dict_part["regional_bounds"]["d"][0], dict_part["codon_positions"]["v"], dict_part["d_5p_del"])
    else:
        enter_dict["DGene_RF"] = 0

    #check_rf = [i for i in dict_part["in_frames"]][0]
    check_rf = dict_part["cdr3_length"] %3 == 0
    if check_rf == True:
        enter_dict["CDR3ReadFrame"] = 0
    else:
        enter_dict["CDR3ReadFrame"] = 1

    SMUA_dict = {}
    # now create the CDR3 + J part of the SMUA string
    SMUA_CDR3_1 = "V" * (dict_part["regional_bounds"]["v"][1] - (dict_part["codon_positions"]["v"] + 3)) + "n" * enter_dict["nVDNs"]
    SMUA_CDR3_2 = "D" * (dict_part["regional_bounds"]["d"][1] - dict_part["regional_bounds"]["d"][0])
    SMUA_CDR3_3 = "n" * enter_dict["nDJNs"] + "J" * (dict_part["codon_positions"]["j"] - dict_part["regional_bounds"]["j"][0])
    fw4 = "4" * (dict_part["regional_bounds"]["j"][1] - (dict_part["codon_positions"]["j"]))
    SMUA_rest = SMUA_CDR3_1 + SMUA_CDR3_2 + SMUA_CDR3_3 + fw4

    SMUA_dict["end"] = SMUA_rest

    # now do the conserved AA@Cys2 and AA@Jinv
    enter_dict["AA@Cys2"] = str(Seq(dict_part["cdr3_seqs"][0][:3]).translate())
    enter_dict["AA@JInv"] = str(Seq(dict_part["cdr3_seqs"][0][-3:]).translate())

    # set up for the alignment 
    start_V = dict_part["regional_bounds"]["v"][0]
    end_J = dict_part["regional_bounds"]["j"][1]

    # check INDELS
    indels_flag  = dict_part["has_shm_indels"][0]
    re_compile = re.compile("\.")

    if indels_flag:

        obs_seq_gapped = dict_part["qr_gap_seqs"][0]
        naive_seq_gapped = dict_part["gl_gap_seqs"][0]

        all_dels = len(re.findall(re_compile, str(obs_seq_gapped)))
        all_ins = len(re.findall(re_compile, str(naive_seq_gapped)))

        total = all_dels + all_ins
        obs_seq_gapped = obs_seq_gapped[start_V:end_J+total]
        naive_seq_gapped = naive_seq_gapped[start_V:end_J+total]

        # observed replace the "." with gaps
        obs_seq_gapped = str(obs_seq_gapped).replace(".", "-")

        # naive replace the "." with gaps
        naive_seq_gapped = str(naive_seq_gapped).replace(".", "-")

        fasta_dict = {"naive_aligned": naive_seq_gapped, "obs_aligned": obs_seq_gapped}

    else:

        naive_seq  = dict_part["naive_seq"][start_V:end_J]
        obs_seq = dict_part["input_seqs"][0][start_V:end_J]
        fasta_dict = {"naive_aligned": naive_seq, "obs_aligned": obs_seq}
    
    # are there partial V or partial J
    # find the number of Ns and compare to the framework insertions

    start_Ns = re.search("^N+", dict_part["input_seqs"][0])
    end_Ns = re.search("N+$", dict_part["input_seqs"][0])

    partial_V = False
    partial_J = False

    if start_Ns:
        if len(dict_part["fv_insertion"]) - start_Ns.span()[1] <  0:
            partial_V = True
    if end_Ns:
        if len(dict_part["jf_insertion"]) - (end_Ns.span()[1] - end_Ns.span()[0]) < 0:
            partial_J = True

    enter_dict["PARTIAL_V"] = partial_V
    enter_dict["PARTIAL_J"] = partial_J

    # enter the information into the dictionary

    for i in ids_in:
        try:
            return_dict[i] = {}
            return_dict[i]["VALID"]= True
            return_dict[i]["fasta_dict"] = fasta_dict
            return_dict[i]["SMUA_dict"] = SMUA_dict
            
            for k in out_col_names:
                if k in enter_dict.keys():
                    return_dict[i][k] = enter_dict[k]
        except:
            #print("error")
            logging.error("Error with this ID:" + ids_in)
            #print(return_dict)

    return return_dict


def partis_parse2(out_dict_v1, imgt_info_dict, chain, species):

    # for each key in the info dict
    for k, v in out_dict_v1.items():

        if v["VALID"] == True:

            if species == "rearrange" and re.search("Human", v["VGene"]):
                sub_dict = imgt_info_dict["human"]
                v_gene_search = v["VGene"].replace("Human","").replace("*","_").replace("S","-S")
            elif species == "rearrange" and not re.search("Human", v["VGene"]):
                sub_dict = imgt_info_dict["mouse"]
                v_gene_search = v["VGene"].replace("*","_").replace("S","-S")

            elif species == "rhesus":
                sub_dict = imgt_info_dict[species]
                v_gene_search = v["VGene"].replace("*","_")
                if v_gene_search in sub_dict.keys():
                    v_gene_search = v_gene_search
                elif len(v_gene_search.split("_"))>2 and not re.search("NL", v_gene_search):
                    v_gene_search = "_".join(v_gene_search.split("_")[0:2])
                elif len(v_gene_search.split("_"))>3 and re.search("NL", v_gene_search):
                    v_gene_search = "_".join(v_gene_search.split("_")[0:3])

            else:
                sub_dict = imgt_info_dict[species] 
                v_gene_search = v["VGene"].replace("*","_")#.replace("S", "-S")
                
            if v_gene_search in sub_dict.keys():
                
                SMUA_dict_orig = v["SMUA_dict"]
                for kk, vv in sub_dict[v_gene_search].items():
                    out_dict_v1[k]["SMUA_dict"][kk] = vv
            else:
                logging.error("hmmm no info for this V gene:" + v_gene_search + "in the IMGT input file. Original name: " + v["VGene"] + " Exiting!!!")
                sys.exit()

    return out_dict_v1


def partis_parse(partis_dict, species, chain):
    out_dict_p1 = {}
    seq_valid = {}
    for v in partis_dict:
        # check to see if it is valide
        if v["invalid"] == True:
            invalid_ids = v["unique_ids"]

            for inv in invalid_ids:
                out_dict_p1[inv] = {}
                out_dict_p1[inv]["VALID"] = False
                for i in out_col_names:
                    if i == "PRODUCTIVE":
                        out_dict_p1[inv][i] = False
                    elif i in ["PARTIAL_V", "PARTIAL_J"]:
                        out_dict_p1[inv][i] = True
                    else:
                        out_dict_p1[inv][i] = 0
        # if the sequence is valid
        else:

            if len([i for i in v["duplicates"] if i]):
                valid_ids = v["unique_ids"]
                for d in v["duplicates"]:
                    valid_ids = valid_ids + d
            else:
                valid_ids = v["unique_ids"]

            # put all the valid sequences into a dictionary 
            for i in valid_ids:
                seq_valid[i] = v["naive_seq"][v["regional_bounds"]["v"][0]:]

            out_dict_p1 = valid_partis_p1(valid_ids, out_dict_p1, v, chain)
    
    return out_dict_p1, seq_valid
    
def yaml_load(yaml_path):

    stream = open(yaml_path, "r")
    partis_dict_full = yaml.load(stream, Loader=yamlloader.ordereddict.CLoader)
    print("yaml file is loaded")
    return partis_dict_full

def args_parse_function(args_in):

    all_args = args_in
    
    all_args.add_argument("-p","--path_yaml", default = [], required = True, help = "path to the yaml file being parsed a partis RS file")
    all_args.add_argument("-o", "--final_out_path", type = Path, default = Path.cwd(), required = False, help = "absolute path to folder where output goes OR if not set the current working directory")
    all_args.add_argument("-c", "--chain", default = [], required = True, choices = ["heavy", "lambda", "kappa"], help = "chain for the markup strings and functional checks")
    all_args.add_argument("-s", "--species", default = "rearrange", required = True, choices = ["mouse","human","rhesus","rearrange"], help = "species for the markup strings")
    all_args.add_argument("-r", "--rf_flag", default = 0, required = False, choices = [1, 0], type = int, help = "flag for if you want the CDR3 reading frame to be returned.")
    all_args.add_argument("-n", "--ngs_flag", default = 0, required = False, choices = [1, 0], type = int, help = "If set to 1, then the output file names start with the first *dot* separated part of the input filename. Otherwise VH, VK, VL")
    args_return = all_args.parse_args()
    return args_return


def main():

    args_return = args_parse_function(argparse.ArgumentParser("This script takes in a .yaml file after running partis annotation, and returns Cloanalyst style SimpleMarkedUAs and Recombination Summary files. This code requires multiple python modules detailed in the requirements.txt file. It also required having *cython* in your environment to speed up .yaml loading. The RecommbinationSummaries file has a line for all of the sequences that are considered valid by *partis*, and the RS.RF.fxnl.prod files has a line for the sequences that have SMUA and are functional (invariants are invariant, CDR3 length is mod3, there are no non mod3 insertions, and, unless the sequence is short at the beginning, the deletions must also be mod3) and productive. There is also a partis_eval file that has information for every input sequence. An error log file is written to for every sequence that partis finds a valid recombination for but a SMUA string that is the same length as the sequence cannot be generated for. These lines will have NaNs some of the RecombinationSummaries values. If the RF flag is set, then there will be an additional column that is the D gene reading frame. There are two options for output file naming. This script is not currently compatible with rhesus light chains."))
    partis_dict = yaml_load(args_return.path_yaml)
    
    Path(args_return.final_out_path).mkdir(exist_ok = True)

    if args_return.ngs_flag:
        # find the first part of the input name 
        yaml_file = Path(args_return.path_yaml).stem
        out_start = yaml_file.split(".")[0] + "."
    else:
        out_start = out_names[args_return.chain] + "_"


    logging.basicConfig(filename = "error_log.txt", level = logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')

    # get the first part of event information and the sequences
    out_dict_v1, seq_dict = partis_parse(list(partis_dict["events"]), args_return.species, args_return.chain)

    # get the IMGT info
    if args_return.species == "rearrange":
        # need both the human and mouse hashes 
        imgt_v_dict_MOUSE = yaml_load(imgt_hash_paths["mouse"][args_return.chain])
        imgt_v_dict_HUMAN = yaml_load(imgt_hash_paths["human"][args_return.chain])
        imgt_v_dict_DICTS = {"mouse":imgt_v_dict_MOUSE, "human":imgt_v_dict_HUMAN}

    else:
        imgt_v_dict_DICTS = {args_return.species: yaml_load(imgt_hash_paths[args_return.species][args_return.chain])}


    # combine the IMGT info and more event information 
    out_dict = partis_parse2(out_dict_v1, imgt_v_dict_DICTS, args_return.chain, args_return.species)


    out_dict_PD = pd.DataFrame.from_dict(out_dict, orient = "index")
    out_dict_PD.index.name = "UID"


    fasta_file_path = Path(args_return.final_out_path).joinpath(out_start + "SimpleMarkedUAs.fasta")
    
    cys_info, SMUA_info, V_ins_info, V_dels_info, V_bases_info, V_subs_info, V_MuFreq_info = SMUA_write(out_dict, fasta_file_path, args_return.chain)

    out_dict_PD["AA@Cys1"] = cys_info
    out_dict_PD["SMUA"] = SMUA_info
    out_dict_PD["#VDeletions"] = V_dels_info
    out_dict_PD["#VInsertions"] = V_ins_info
    out_dict_PD["#VBases"] = V_bases_info
    out_dict_PD["#VSubstitutions"] = V_subs_info
    out_dict_PD["MuFreq"] = V_MuFreq_info
    out_dict_PD["MDL"] = 0

    out_dict_PD["FUNCTIONAL"] = False

    out_dict_PD_wFXNL = find_fxnl(out_dict_PD, args_return.chain)
    fxnl_prod_seqs = [i for i,r in out_dict_PD_wFXNL.iterrows() if (r["PRODUCTIVE"] == True and r["FUNCTIONAL"] == True and r["SMUA"])]


    # now we want the RS file that only has that are valid

    if args_return.rf_flag == 0:

        out_dict_filt = out_dict_PD_wFXNL.loc[fxnl_prod_seqs, out_col_names_RS]
        out_dict_RS_clean = out_dict_PD.loc[out_dict_PD["VALID"] == 1, out_col_names_extra + out_col_names_RS]

        file_all_path = Path(args_return.final_out_path).joinpath(out_start + "RecombinationSummaries.txt")
        file_filtered_path = Path(args_return.final_out_path).joinpath(out_start + "RS.fxnl.prod.txt")
   
    else:

        out_dict_filt = out_dict_PD_wFXNL.loc[fxnl_prod_seqs, out_col_names_RS + ["DGene_RF"]]
        out_dict_RS_clean = out_dict_PD.loc[out_dict_PD["VALID"] == 1, out_col_names_extra + out_col_names_RS + ["DGene_RF"]]

        file_all_path = Path(args_return.final_out_path).joinpath(out_start + "RecombinationSummaries.RF.txt")
        file_filtered_path = Path(args_return.final_out_path).joinpath(out_start + "RS.RF.fxnl.prod.txt")

    out_dict_extra = out_dict_PD_wFXNL.loc[:, out_col_names_extra]
    extra_info_path = Path(args_return.final_out_path).joinpath(out_start + "partis_eval.txt")

    out_dict_RS_clean.to_csv(file_all_path, index = True, sep = "\t")
    out_dict_filt.to_csv(file_filtered_path, index = True, sep = "\t")    
    out_dict_extra.to_csv(extra_info_path, index = True, sep = "\t")

if __name__ == "__main__":
    main()

