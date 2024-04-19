# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2024-04-06 12:53:55
# @Last Modified by:   chem
# @Last Modified time: 2024-04-10 16:11:17


import os, sys
import re
import json
import pandas as pd
import subprocess
import fastaparser

# from cleavage_pred.protai.annotation.data import DatasetGenerator as ADG
from pathlib import Path

folder = os.path.abspath(__file__)
sys.path.append(os.path.join(os.path.dirname(folder), "cleavage_pred"))
from protai.annotation.data import DatasetGenerator as ADG


# def split_sequence(sequence):
#     # Search for the first occurrence of "GG", "GA", or "GS" within coordinates 10-40
#     match = re.search(r"(GG|GA|GS)", sequence[10:41])

#     # If any of "GG", "GA", or "GS" is found within the specified range
#     if match:
#         # Find the index of the match within the original sequence
#         min_index = match.start() + 10  # Add 10 to account for the offset
#         # Split the sequence into two fragments
#         fragment1 = sequence[:min_index + 2]  # Include the matched segment in the first fragment
#         fragment2 = sequence[min_index + 2:]  # Skip the two characters of the pattern
#         return fragment1, fragment2
#     else:
#         # None of "GG", "GA", or "GS" found within the specified range
#         return None, None


def split_leader(sequence, cleavageSite):
    # For instance: CS pos: 30-31. Pr: 0.7732
    leader, core = None, None
    if not pd.isna(cleavageSite):
        posit = re.search("CS pos: (.+)-", cleavageSite).group(1)
        leader = sequence[: int(posit)]
        core = sequence[int(posit) :]
    return leader, core


def leader(infa, injsonSec, injsonGG, outf, threshold):
    """
    infa: predicted precursors with Sec leader
    outf: output folder
    threshold: processor core
    """
    with open(injsonSec) as f1, open(injsonGG) as f2:
        dict_pre_sec = json.load(f1)
        dict_pre_gg = json.load(f2)
    list_pre_sec = [(key, value) for key, value in dict_pre_sec.items()]
    list_pre_gg = [(key, value) for key, value in dict_pre_gg.items()]
    # Convert list to DataFrame
    df_pre_sec = pd.DataFrame(list_pre_sec, columns=["ID", "Sequence"])
    df_pre_gg = pd.DataFrame(list_pre_gg, columns=["ID", "Sequence"])

    ###################################
    # run signalp6
    ###################################
    if os.path.getsize(infa) > 0:
        command = f"signalp6 --fastafile {infa} --organism other --output_dir {outf} --write_procs {threshold} --format txt --mode fast"
        subprocess.run(command, shell=True)
        # parase prediction outcome
        out_prediction = os.path.join(outf, "prediction_results.txt")
        df_signalp6 = pd.read_table(out_prediction, skiprows=1, header=0)
        df_signalp6_m = pd.merge(
            df_pre_sec,
            df_signalp6[["# ID", "CS Position"]],
            how="left",
            left_on="ID",
            right_on="# ID",
        )
        df_signalp6_m[["leader_sec", "core_sec"]] = df_signalp6_m.apply(
            lambda row: split_leader(row["Sequence"], row["CS Position"]),
            axis=1,
            result_type="expand",
        )
        df_signalp6_f = df_signalp6_m[["ID", "Sequence", "leader_sec", "core_sec"]]
        df_signalp6_f = df_signalp6_f.rename(columns={"ID": "Uniq_ID"})
    else:
        df_signalp6_f = pd.DataFrame(
            columns=["Uniq_ID", "Sequence", "leader_sec", "core_sec"]
        )

    ###################################
    ## double-glycine leader prediction using nlpprecursor
    ####################################
    current_path = os.path.split(os.path.realpath(__file__))[0]
    out_gg_json = os.path.join(outf, "prediction_double_GG_leader.json")
    train_code = os.path.join(current_path, "cleavage_pred/train.py")
    models_dir = Path(current_path + "/cleavage_pred/training_data/annotation/")
    annot_model_path = models_dir / "model.p"
    annot_vocab_path = models_dir / "vocab.pkl"

    if list_pre_gg:  # judge whether the input is empty
        # judge whether the model is available
        if not os.path.exists(annot_model_path):
            train_command = f"python {train_code}"
            subprocess.run(train_command, shell=True)

        pred_in = [{"name": k, "sequence": v} for k, v in dict_pre_gg.items()]

        cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, pred_in)
        with open(out_gg_json, "w") as fo:
            json.dump(cleavage_predictions, fo, indent=4)

        list_gg = []

        for d in cleavage_predictions:
            seqname = d["name"]
            sequence = dict_pre_gg[seqname]
            core_gg, start_p, leader_gg = None, None, None
            if d["cleavage_prediction"]["status"] == "success":
                core_gg = d["cleavage_prediction"]["sequence"]
                start_p = int(d["cleavage_prediction"]["start"])
                leader_gg = sequence[:start_p]
            list_gg.append(
                {
                    "Uniq_ID": seqname,
                    "Sequence": sequence,
                    "leader_gg": leader_gg,
                    "core_gg": core_gg,
                }
            )

        df_gg = pd.DataFrame(list_gg)
    else:
        df_gg = pd.DataFrame(columns=["Uniq_ID", "Sequence", "leader_gg", "core_gg"])

    ###################################
    ## combine two predictions
    ####################################
    df_c = pd.merge(df_signalp6_f, df_gg, on=["Uniq_ID", "Sequence"], how="outer")

    return df_c
