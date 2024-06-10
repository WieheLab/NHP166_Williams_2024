#!/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import substitution_matrices, PairwiseAligner
import gzip
import pandas as pd
from pathlib import Path
import re
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor

aligner = PairwiseAligner()
aligner.match_score = 5
aligner.mismatch_score = -4
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
aligner.mode = "local"

def load_SMUA(SMUA_file):

    SMUA_dict = {}

    # import the SMUA file depending on whether gzip or not
    if (str(SMUA_file).split(".")[-1] == "gz"):
        with gzip.open(SMUA_file, "rt") as handle:
            #SMUA_dict = {rec.description : rec.seq for rec in SeqIO.parse(handle, "fasta")}
            for c, rec in enumerate(SeqIO.parse(handle, "fasta")):
                if c==1:
                    if not re.search("UCA", rec.description):
                        print("bad news! you did not provide an SMUA! try again!")
                        sys.exit(1)
                if c%3 == 0:
                    SMUA_dict[rec.id] = rec.seq
    else:
        with open(SMUA_file, "rt") as handle:
            for c, rec in enumerate(SeqIO.parse(handle, "fasta")):
                if c==1:
                    if not re.search("UCA", rec.description):
                        print("bad news! you did not provide an SMUA! try again!")
                        sys.exit(1)
                if c%3 == 0:
                    SMUA_dict[rec.id] = rec.seq
    return SMUA_dict

def load_fasta(fasta_file):

    # import the fasta  file depending on whether gzip or not
    if (str(fasta_file).split(".")[-1] == "gz"):
        with gzip.open(fasta_file, "rt") as handle:
            fasta_dict = {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}
    else:
        with open(fasta_file, "rt") as handle:
            fasta_dict = {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}

    if len(list(set(fasta_dict.keys()))) != len(list(fasta_dict.keys())):
        
        print("bad news! you fasta file sequences ids are not unique. this could be a problem!")
        sys.exit(1)
    elif sum([1 if re.search("UCA", k) else 0 for k in list(fasta_dict.keys())]) > 0:
        print("it looks like you might have entered an SMUA where it should have been a full length fasta sequence file!")
        sys.exit(1)

    return fasta_dict

def VDJ_constants(SMUA_dict, full_seqs_dict):

    seq_parts_dict = {}
    for k, v in SMUA_dict.items():
        seq_parts_dict[k] = {}

        if k in full_seqs_dict.keys():
            subject = str(v).replace("-","")
            subject_AMBIG = re.escape(subject).replace("\\N", "[ATGCN]")
            full_seq = str(full_seqs_dict[k])
            full_seq_AMBIG = re.escape(full_seq).replace("\\N", "[ATGCN]")

            parts = re.search(subject_AMBIG, full_seq_AMBIG)
            #parts = re.search(str(v).replace("-",""), str(full_seqs_dict[k]).replace("N", ""))

            if parts:
                seq_parts_dict[k]["is_parts"]= True
                seq_parts_dict[k]["3p"] = str(full_seqs_dict[k])[parts.span()[1]:]
            else:
                seq_parts_dict[k]["is_parts"] = False
        else:
            print("An SMUA sequence is not in the full length fasta file. This is a bad sign. EXITING!")
            sys.exit()

    return seq_parts_dict

def find_isotype(seq_parts_dict, constant_fasta_dict, num_return):
    good = 0
    store_all = {}
    for k, v in seq_parts_dict.items():
        
        if v["is_parts"] ==True and len(v["3p"])>0:
            store_score = {}
            good+=1 
            for k_i, v_i in constant_fasta_dict.items():
                out = aligner.align(Seq(v["3p"].lower()), Seq(v_i.lower()))
                store_score[k_i] = {}
                store_score[k_i]["best_score"] = out.score
                store_score[k_i]["best_possible_score"] = len(v["3p"])*int(aligner.match_score)

            store_score_DF = pd.DataFrame.from_dict(store_score, orient = "index").sort_values(by = "best_score", ascending = False)
            store_all[k] = {}
            

            store_all[k]["best_possible_score"] = store_score_DF.loc[store_score_DF.index.tolist()[0], "best_possible_score"]
            for c in range(min(num_return, store_score_DF.shape[0])):
                store_all[k]["match_" + str(c+1)] = store_score_DF.index.tolist()[c]
                store_all[k]["score_" + str(c+1)] = store_score_DF.loc[store_score_DF.index.tolist()[c], "best_score"]
            
            store_all[k]["after_VDJ"] = str(v["3p"])
        
        elif v["is_parts"] == True and len(v["3p"])==0:
            store_all[k] = {}
            store_all[k]["best_possible_score"] = 0
            store_all[k]["match_1"] = "No_Sequence_3p_of_VDJ"
        
        elif v["is_parts"] == False:

            store_all[k] = {}
            store_all[k]["best_possible_score"] = 0
            store_all[k]["match_1"] = "Missing"

        else:
            
            print("Something weird is happening with the regular expression search of the VDJ against the full length sequence!")
            sys.exit(1)

    store_all_DF = pd.DataFrame.from_dict(store_all, orient= "index")
    store_all_DF.index.name = "UID"

    if good == 0:
        print("hrmm. it looks like maybe the fasta sequences that you supplied don't have constant sequence portions?!")
        sys.exit()

    num_actual_return = sum([1 if re.search("match", l) else 0 for l in store_all_DF.columns.tolist()]) 

    return store_all_DF, num_actual_return

def find_isotype_chunk(seq_parts_chunk, constant_fasta_dict, num_return):
    store_all_DF_chunk, c = find_isotype(seq_parts_chunk, constant_fasta_dict, num_return)
    return store_all_DF_chunk, c

def args_parse_function(args_in):

    all_args = args_in
    all_args.add_argument("-s", "--path_SMUA", default = [], required = True, help = "absolute path to the SMUA file")
    all_args.add_argument("-l", "--path_FL_seqs", default = [], required = True, help = "absolute path the fasta file of full length sequences")
    all_args.add_argument("-o", "--output_path", type=Path, default=Path.cwd(), required = False, help = "Defaul is current working directory, but if set is the full path to the output folder where the the isotype result file will be written")
    all_args.add_argument("-f", "--constant_fasta", type = str, required = True, help = "fasta files with the constant regions for the heavy chain (human, rhesus, mouse files can be found here: /datacommons/dhvi/partis_ASAP/speciesCH).")
    all_args.add_argument("-n", "--num_hits", type = int, required = False, default = 3, help = "number of best isotype hits to return (default is 3)")
    all_args.add_argument("-t", "--tasks", type = int, required = False, default = 1, help = "the number of tasks that you want to use. If not set, the value = 1 and thus it is not paralleleize. If >1 then make sure the correct resources are assigned in the job script!")
    all_args.add_argument("-g", "--ngs_flag", default = 0, required = False, choices = [1,0], type = int, help="If set to 1, then the output file names start with the first *dot* separated part of the input full length fasta filename. Otherwise VH")
    args_return = all_args.parse_args()
    return args_return


def main():

    args_return = args_parse_function(argparse.ArgumentParser(description = "This is a script that identifies the isotype of a heavy chain sequence. The script does this by finding the sequence that is 3' of the VDJ sequence, thus, calls are returned for each sequence that has an SMUA line, but not for full length sequences without an SMUA. The returned csv file will have the top 3 hits if a different value is not specified. If the number of hits requested is larger than the number of constant sequences in the fasta file, then the hit information for each sequence in the constant fasta file is returned. To identify the isotype a local alignment is used. The scoring matrix is at the top of the code. This script CAN be parallelized using the TASKS argument. If you set tasks >1 make sure you have assigned at least that many tasks on the same node. This script IS compatible with nucleotide ambiguity codes in that it will allow for mismatches between an N ambiguity and any AGCT character. It is not compatible with non-N ambiguity codes."))

    SMUA_dict = load_SMUA(args_return.path_SMUA)
    full_seqs_dict = load_fasta(args_return.path_FL_seqs)

    # get the name for the output file
    if args_return.ngs_flag:
        # find the first part of the input name
        full_seqs_file = Path(args_return.path_FL_seqs).stem
        out_start = full_seqs_file.split(".")[0] + "."
    else:
        out_start = "VH_"

    # i think the first thing to do is to go through all the sequences in the SMUA and all the sequences in the FASTA and find the end part
    seq_parts_dict = VDJ_constants(SMUA_dict, full_seqs_dict)
    isotypes_dict = load_fasta(args_return.constant_fasta)

# Determine the number of processes to use
    num_processes = args_return.tasks  # or another number that makes sense for your machine

    # Create a list of keys to distribute to workers
    keys = list(seq_parts_dict.keys())
    num_processes_submit = min(len(keys), num_processes)
    # Split the keys into chunks for each process
    chunks = [keys[i::num_processes] for i in range(num_processes_submit)]

    with ProcessPoolExecutor(max_workers=num_processes_submit) as executor:
        # Submit the find_isotype function to the executor, for each chunk of keys
        futures = [executor.submit(find_isotype_chunk, {k: seq_parts_dict[k] for k in chunk}, isotypes_dict, args_return.num_hits) for chunk in chunks]

        # Combine the results as they complete
        isotypes_PD_combined = pd.DataFrame()
        for future in concurrent.futures.as_completed(futures):
            isotypes_PD_chunk, num_real_return = future.result()
            isotypes_PD_combined = pd.concat([isotypes_PD_combined, isotypes_PD_chunk])

    isotypes_PD_combined = isotypes_PD_combined.sort_index()

    if len(keys) == isotypes_PD_combined.shape[0]:
        print("great! information was returned for every VDJ sequence!")
    else:
        print("baddd! informaiton was NOT returend for every VDJ sequence!")
        sys.exit(1)

#    isotypes_PD, num_real_return = find_isotype(seq_parts_dict, isotypes_dict, args_return.num_hits)
    
    isotypes_out = Path(args_return.output_path).joinpath(out_start + "isotypes_top" + str(num_real_return) + ".csv")

    isotypes_PD_combined.to_csv(isotypes_out)

if __name__ == "__main__":
    main()
    
