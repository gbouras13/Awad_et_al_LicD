#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: gbouras13 adapting code from https://colab.research.google.com/drive/1A2XJ_oUtlmIfU3XXmev5dzJUNxqR6VV9?usp=sharing#scrollTo=d8bbf189

Runs deposcope on bulk multiFASTA input proteins - both the classifier and the amino-acid level domain assignment

* https://github.com/dimiboeckaerts/DepoScope 
* https://doi.org/10.1371/journal.pcbi.1011831

Models available on zenodo - v2 https://zenodo.org/records/10957073

example command 

python run_deposcope.py -i all_phold_structures_aa.faa -o all_phold_structures_output -e esm2_t12_35M_UR50D__fulltrain__finetuneddepolymerase.2103.4_labels/checkpoint-2255 -d Deposcope.esm2_t12_35M_UR50D.2203.full.model

Classifies depolymerases into 3 classes of depolymerase - the paper suggests that the first two classes are easier to predict (unsurprisingly given the training data)

1. beta-helix
2. beta propeller
3. Triple helix

"""

import os
import sys
import shutil
from argparse import RawTextHelpFormatter
import sys
import torch
import argparse
import pandas as pd
import torch.nn.functional as F
from torch import nn
from Bio import SeqIO
from tqdm import tqdm
from transformers import AutoModelForTokenClassification, AutoTokenizer

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


def get_input():
    parser = argparse.ArgumentParser(
        description="clean_efams.py",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--infile", action="store", help="Input multifasta amino acid format."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Directory to write the output to.",
        default=os.path.join(os.getcwd(), "output/"),
    )
    parser.add_argument(
        "-e",
        "--esm2_model_path",
        action="store",
        help="finetuned esm2_model_path - specify checkpoint-6015",
        default=os.path.join(
            os.getcwd(),
            "esm2_t12_35M_UR50D-finetuned-depolymerase.labels_4/checkpoint-6015",
        ),
    )
    parser.add_argument(
        "-d",
        "--DpoDetection_path",
        action="store",
        help="DepoDetection.T12.4Labels.1908.model path",
        default=os.path.join(os.getcwd(), "DepoDetection.T12.4Labels.1908.model"),
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    parser.add_argument(
         "--cpu", help="CPU only mode", action="store_true"
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

    # make outdir if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)

    # add GPU support
    global device

    if args.cpu:
        print(f"using device cpu as --cpu is {args.cpu}") 
        device = torch.device("cpu")
        dev_name = "cpu"
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda:0")
            dev_name = "cuda:0"
       # check for apple silicon/metal
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
            dev_name = "mps"        
        else:
            device = torch.device("cpu")
            dev_name = "cpu"

    print("Using device: {}".format(dev_name))

    # Load models from their paths
    esm2_model_path = args.esm2_model_path
    DpoDetection_path = args.DpoDetection_path

    tokenizer = AutoTokenizer.from_pretrained(esm2_model_path)
    esm2_finetuned = AutoModelForTokenClassification.from_pretrained(esm2_model_path).to(device)

    # define the class related to the models to run them
    class Dpo_classifier(nn.Module):
        def __init__(self, pretrained_model):
            super(Dpo_classifier, self).__init__()
            self.max_length = 1024
            self.pretrained_model = pretrained_model
            self.conv1 = nn.Conv1d(
                1, 64, kernel_size=5, stride=1
            )  # Convolutional layer
            self.conv2 = nn.Conv1d(
                64, 128, kernel_size=5, stride=1
            )  # Convolutional layer
            self.fc1 = nn.Linear(
                128 * (self.max_length - 2 * (5 - 1)), 32
            )  # calculate the output shape after 2 conv layers
            self.classifier = nn.Linear(32, 1)  # Binary classification
            self.classifier.to(device)

        def make_prediction(self, fasta_txt):
            input_ids = tokenizer.encode(
                fasta_txt, truncation=True, return_tensors="pt"
            ).to(device) 
            with torch.no_grad():
                outputs = self.pretrained_model(input_ids)
                probs = torch.nn.functional.softmax(outputs.logits, dim=-1)
                token_probs, token_ids = torch.max(probs, dim=-1)
                tokens = token_ids.view(1, -1)  # ensure 2D shape
                return tokens

        def pad_or_truncate(self, tokens):
            if tokens.size(1) < self.max_length:
                tokens = F.pad(tokens, (0, self.max_length - tokens.size(1)))
            elif tokens.size(1) > self.max_length:
                tokens = tokens[:, : self.max_length]
            return tokens

        def forward(self, sequences):
            batch_size = len(sequences)
            tokens_batch = []
            for seq in sequences:
                tokens = self.make_prediction(seq)
                tokens = self.pad_or_truncate(tokens)
                tokens_batch.append(tokens)

            outputs = torch.cat(tokens_batch).view(
                batch_size, 1, self.max_length
            )  # ensure 3D shape
            outputs = outputs.float()  # Convert to float

            out = F.relu(self.conv1(outputs))
            out = F.relu(self.conv2(out))
            out = out.view(batch_size, -1)  # Flatten the tensor
            out = F.relu(self.fc1(out))
            out = self.classifier(out)
            return out, outputs

    # define a helper sequence to make predictions with
    def predict_sequence(model, sequence):
        model.eval()
        with torch.no_grad():
            sequence = [
                sequence
            ]  # Wrap the sequence in a list to match the model's input format
            outputs, sequence_outputs = model(sequence)
            probas = torch.sigmoid(
                outputs
            )  # Apply sigmoid activation for binary classification
            predictions = (probas > 0.5).float()  # Convert probabilities to binary predictions
            sequence_outputs_list = sequence_outputs.cpu().numpy().tolist()[0][0]
            prob_predicted = probas[0].item()
            return (predictions.item(), prob_predicted), sequence_outputs_list

    # instantiate the model
    model_classifier = Dpo_classifier(
        esm2_finetuned
    ).to(device)  # Create an instance of Dpo_classifier and GPU support
    model_classifier.load_state_dict(
        torch.load(DpoDetection_path), strict=False
    )  # Load the saved weights



    model_classifier.eval()  # Set the model to evaluation mode for inference

    ####### run predictions #######

    # get the ids and records
    def read_multi_fasta_to_df(filename):
        data = []
        for record in SeqIO.parse(filename, "fasta"):
            data.append([record.id, str(record.seq)])
        df = pd.DataFrame(data, columns=["protein_id", "protein_sequence"])
        return df

    filename = args.infile
    phage_genes_df = read_multi_fasta_to_df(filename)
    phage_data = list(zip(phage_genes_df["protein_id"], phage_genes_df["protein_sequence"]))


    print("Running Deposcope")
 
    # # ********************************************
    # token-level classes correspond from here https://github.com/dimiboeckaerts/DepoScope/blob/main/Training/2.PT1_esm2_FineTuning_Token_class_review_server_fulltrain.ipynb
    # beta-helix = 1
    # beta propeller = 2
    # triple helix = 3

    def find_ranges(lst):
        ranges_1 = []
        ranges_2 = []
        ranges_3 = []

        start_1 = start_2 = start_3 = None
        start = None

        for i, num in enumerate(lst):
            if num == 1:
                if start_1 is None:
                    start_1 = i + 1  # start position (1-indexed)
            else:
                if start_1 is not None:
                    ranges_1.append((start_1, i))
                    start_1 = None

            if num == 2:
                if start_2 is None:
                    start_2 = i + 1  # start position (1-indexed)
            else:
                if start_2 is not None:
                    ranges_2.append((start_2, i))
                    start_2 = None

            if num == 3:
                if start_3 is None:
                    start_3 = i + 1  # start position (1-indexed)
            else:
                if start_3 is not None:
                    ranges_3.append((start_3, i))
                    start_3 = None

        if start_1 is not None:  # handle the case if the last element is part of a range
            ranges_1.append((start_1, len(lst)))

        if start_2 is not None:  # handle the case if the last element is part of a range
            ranges_2.append((start_2, len(lst)))

        if start_3 is not None:  # handle the case if the last element is part of a range
            ranges_3.append((start_3, len(lst)))

        return ranges_1, ranges_2, ranges_3

    results = []

    progress_bar = tqdm(phage_data, total=len(phage_genes_df["protein_sequence"]))
    for protein_id, protein_seq in progress_bar:
        prediction, sequence_outputs = predict_sequence(model_classifier, protein_seq)
        if prediction[0] == 0:
            depolymerase = False
            depolymerase_types = ""
            depolymerase_ranges = ""
        else:
            depolymerase = True
            depolymerase_types_list = []
            depolymerase_ranges_list = []
            print(f"{protein_id} predicted depolymerase {depolymerase} with probability {round(prediction[1], 6)}")
        if depolymerase:
            #print(sequence_outputs)
            ranges_1, ranges_2, ranges_3 = find_ranges(sequence_outputs)
            if ranges_1:
                depolymerase_types_list.append("beta-helix")
                depolymerase_ranges_list.append(ranges_1)
            if ranges_2:
                depolymerase_types_list.append("beta propeller")
                depolymerase_ranges_list.append(ranges_2)
            if ranges_3:
                depolymerase_types_list.append("triple helix")
                depolymerase_ranges_list.append(ranges_3)
                # Append the result to the list

            # Join the lists into comma-separated strings
            depolymerase_types = ", ".join(depolymerase_types_list)
            depolymerase_ranges = ", ".join(str(r) for r in depolymerase_ranges_list)
    

        results.append({
            "protein_id": protein_id,
            "depolymerase": depolymerase,
            "depolymerase_probability": round(prediction[1], 6),
            "depolymerase_type": depolymerase_types,
            "depolymerase_range": depolymerase_ranges
        })        


    df = pd.DataFrame(results)

    df.to_csv(os.path.join(args.outdir, "deposcope_output.csv"), index=False)


if __name__ == "__main__":
    main()