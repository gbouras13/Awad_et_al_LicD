#!/usr/bin/env python3
import argparse
import gzip
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import os
import shutil
import sys
import glob
import json
from pathlib import Path
import numpy as np
from Bio.Seq import Seq
import pandas as pd


def get_input():
    usage = "python3 rename_preocess_tophits.py ..."
    parser = argparse.ArgumentParser(
        description="script to rename top ranked pdbs and json from colabfold output (e.g. 0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb) to the name as per its multifasta and get some more metadata.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        help="Input file of AAs. Headers with the names, must be gzipped",
        required=True,
    )
    parser.add_argument(
        "-p", "--predictions_dir", action="store", help="Colabfold batch predictions dir", required=True
    )
    parser.add_argument(
        "-o", "--outdir", action="store", help="Output directory", required=True
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    parser.add_argument(
        "--rename", help="Rename MSAs - MSAs generated with Phoenix", action="store_true"
    )
    
    args = parser.parse_args()

    return args


def main():
    args = get_input()

    if args.force == True:
        if os.path.isdir(args.outdir) == True:
            print(
                f"Removing output directory {args.outdir} as -f or --force was specified."
            )
            shutil.rmtree(args.outdir)
        elif os.path.isfile(args.outdir) == True:
            os.remove(args.outdir)
        else:
            print(
                f"--force was specified even though the output directory {args.outdir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(args.outdir) == True or os.path.isfile(args.outdir) == True:
            print(
                f"The output directory {args.outdir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )
            sys.exit()

    # get fasta names
    protein_headers = []

    with open(args.infile, "rt") as handle:  #"rt" mode ensures text reading

        for record in SeqIO.parse(handle, "fasta"):
            header = f"{record.description}"
            protein_headers.append(header)

    # make outdir if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    pdb_out_dir = os.path.join(args.outdir, "structures")
    msa_out_dir = os.path.join(args.outdir, "msas")
    corrupt_msa_out_dir = os.path.join(args.outdir, "corrupt_msas")
    json_out_dir = os.path.join(args.outdir, "jsons")
    pae_out_dir = os.path.join(args.outdir, "paes")

    os.makedirs(pdb_out_dir, exist_ok=True)
    os.makedirs(msa_out_dir, exist_ok=True)
    os.makedirs(json_out_dir, exist_ok=True)
    os.makedirs(pae_out_dir, exist_ok=True)

    # calculate for every MSA file
    msa_files = [f for f in os.listdir(args.predictions_dir) if f.endswith(".a3m")]

    # above_70_plddt_prots = set()

    plddt_ptm_dict = {}

    for f in msa_files:
        if args.rename:
            index = int(f.replace(".a3m", ""))
            if 0 <= index < len(protein_headers):
                prot_name = protein_headers[index]
            prot_name_in = index

        else:
            prot_name = f.replace(".a3m", "")
            prot_name_in = prot_name

        new_path_msa = os.path.join(msa_out_dir, f"{prot_name}.a3m")
        new_path_pdb = os.path.join(pdb_out_dir, f"{prot_name}.pdb")
        new_path_json = os.path.join(json_out_dir, f"{prot_name}_scores.json")
        new_path_pae = os.path.join(pae_out_dir, f"{prot_name}_pae.json")

        for model_num in range(1, 6):  # check models 1 to 5
            pdb_candidate = os.path.join(args.predictions_dir, f"{prot_name_in}_unrelaxed_rank_001_alphafold2_ptm_model_{model_num}_seed_000.pdb")
            json_candidate = os.path.join(args.predictions_dir, f"{prot_name_in}_scores_rank_001_alphafold2_ptm_model_{model_num}_seed_000.json")
            
            if os.path.exists(pdb_candidate) and os.path.exists(json_candidate):
                pdb_file = pdb_candidate
                json_file = json_candidate
                break  # pick the first one found
        
        # pdb_file = f"{prot_name_in}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
        # json_file =  f"{prot_name_in}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"
        pae_file =  f"{prot_name_in}_predicted_aligned_error_v1.json"
        scores = json.loads(Path(json_file).read_text())

        # get the stats
        plddt = np.asarray(scores["plddt"])
        mean_plddt = round(np.mean(plddt), 2)
        ptm = scores["ptm"]

        shutil.copy(os.path.join(args.predictions_dir, f), new_path_msa)
        shutil.copy( pdb_file, new_path_pdb)
        shutil.copy( json_file, new_path_json)
        shutil.copy(os.path.join(args.predictions_dir, pae_file), new_path_pae)
        
            # if mean_plddt >= 70:
            #     above_70_plddt_prots.add(prot_name)
        # except:
        #     print(f"Corrupt data for {prot_name}")
        #     os.makedirs(corrupt_msa_out_dir, exist_ok=True)
        #     prot_name = f.replace(".a3m", "")
        #     new_path_msa = os.path.join(msa_out_dir, f)
        #     # copy the file if it exists
        #     if os.path.exists(os.path.join(args.predictions_dir, f)):
        #         shutil.copy(os.path.join(args.predictions_dir, f), corrupt_msa_out_dir)
        #     mean_plddt = None
        #     ptm = None


        # Add values to the dictionary
        plddt_ptm_dict[prot_name] = {
                    "mean_plddt": mean_plddt,
                    "ptm": ptm,
                    }



    # Initialize an empty set to store protein headers

    protein_len_dict = {}

    # from https://github.com/sokrypton/ColabFold/blob/c21e1768d18e3608e6e6d99c97134317e7e41c75/colabfold/utils.py#L63C1-L64C86
    def safe_filename(file: str) -> str:
        return "".join([c if c.isalnum() or c in ["_", ".", "-"] else "_" for c in file])

    # outfasta
    # outfasta = os.path.join(args.outdir, f"over_70_plddt.fasta")
    # with open(outfasta, "w") as output_file:
    with open(args.infile, "rt", encoding="utf-8") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            header = safe_filename(record.description)
            record.id = header
            protein_len_dict[header] = {"length": len(Seq(record.seq))}
                # if header in above_70_plddt_prots:
                #     SeqIO.write(record, output_file, "fasta")


    len_df = pd.DataFrame.from_dict(protein_len_dict, orient="index")
    len_df = len_df.rename_axis("protein")
    print(len_df)

    plddt_df = pd.DataFrame.from_dict(plddt_ptm_dict, orient="index")
    plddt_df = plddt_df.rename_axis("protein")
    print(plddt_df)

    # Merge the DataFrames on the 'protein' column - check they are all legit
    merged_df = pd.merge(len_df, plddt_df, on="protein", how="inner")

    print(merged_df)

    # Save DataFrame to a CSV file
    outfile_path = os.path.join(args.outdir, "plddt_ptm_len.tsv")
    merged_df.to_csv(outfile_path, sep="\t", index=True)


if __name__ == "__main__":
    main()
