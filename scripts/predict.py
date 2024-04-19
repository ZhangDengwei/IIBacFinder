# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2023-04-21 20:12:48
# @Last Modified by:   chem
# @Last Modified time: 2024-04-13 10:50:52


import re
import os
import sys
import json
import time
import logging
import argparse
import functools
import subprocess
import warnings
import pandas as pd

# import numpy as np
from calpro import calproperty
from leader_predict import leader
from diamond import alnDiamond
from tqdm import tqdm
from pathlib import Path
from region_annotation import *
from tqdm.contrib.concurrent import process_map
from gbkAnnotation import read_region

warnings.filterwarnings("ignore")


def purge(d, pattern):
    """removes files matching a pattern
    parameters
    ----------
    d
            string, directory path
    pattern
            string, regex
    returns
    ----------
    """
    for f in os.listdir(d):
        if re.search(pattern, f):
            os.remove(os.path.join(d, f))


def scanDir(inDir):
    file_list = []
    suffix = []
    fasta_suffix = [".fas", ".fa", ".fasta", ".faa", ".fna"]
    # Detect the suffix of file
    logging.info("Examining the suffix of input fasta file")
    i = 0
    for file in os.listdir(inDir):
        if Path(file).suffix in fasta_suffix:
            suffix.append(Path(file).suffix)
        i += 1
        if i >= 1000:
            break
    final_suffix = max(suffix, key=suffix.count)
    logging.warning(f"Fasta files with suffiex of '{final_suffix}' would be loaded!!!")

    for root, dirs, files in os.walk(inDir):
        for file in files:
            if file.endswith(final_suffix):
                # check whether the file is existing and not empty
                try:
                    file_path = os.path.abspath(os.path.join(root, file))
                    if os.path.getsize(file_path):
                        file_list.append(file_path)
                except:
                    pass
    return file_list, final_suffix


def run_prodigal(prodigal_exe, model, outDir, fasta):
    """
    model: single or meta
    """
    fasta_prefix = Path(fasta).stem
    out_f1 = os.path.join(outDir, fasta_prefix + ".faa")
    out_f2 = os.path.join(outDir, fasta_prefix + ".fna")
    out_f3 = os.path.join(outDir, fasta_prefix + ".gbk")
    cmd_d = f"{prodigal_exe} -i {fasta} -a {out_f1} -d {out_f2} -o {out_f3} -p {model} > /dev/null 2>&1"
    p = os.system(cmd_d)
    if p == 1:
        sys.exit("Failed gene calling when running prodigal-short!!!")


def region(outRegion):
    df_list_region = []
    for root, dirs, files in os.walk(outRegion):
        for file in files:
            if file.endswith(".tsv"):
                file_path = os.path.abspath(os.path.join(root, file))
                df_region = pd.read_table(file_path)
                # all CDs type in the region
                q = list(set(df_region.CDs_Class))
                # remove "others"
                q[:] = [x for y in q if y != "others" for x in y.split("/")]
                q = list(set(q))
                q.sort()
                df_region_s = df_region[
                    [
                        "Genome",
                        "Region",
                        "Contig",
                        "CDs",
                        "Start",
                        "End",
                        "Strand",
                        "Partial_index",
                        "Start_type",
                        "RBS_motif",
                    ]
                ]

                # df_region_s.loc[:,"Including_elements"] = ",".join(q)
                df_region_s.insert(
                    df_region_s.shape[1], "Including_elements", ",".join(q)
                )
                df_list_region.append(df_region_s)

    result_region = pd.concat(df_list_region)
    return result_region


def leader_combine(leaderType, fullseq, CoreSec, CoreGG):
    core = None
    if leaderType == "Leaderless":
        core = fullseq
    else:
        if isinstance(CoreSec, str):
            core = CoreSec
        else:
            if isinstance(CoreGG, str):
                core = CoreGG
            else:
                pass

    return core


def confidenceP(rbs, partial, rule, leaderType, matureP, leaderSec):
    confidence = None
    if isinstance(rbs, str):
        if partial == "complete":
            if rule == "Both":
                confidence = "High"
            elif rule == "Domain":
                if leaderType in [
                    "Sec-dependent",
                    "Double-glycine",
                    "Double-glycine or Sec-dependent",
                ]:
                    if matureP:
                        confidence = "High"
                    else:
                        confidence = "Low"
                else:
                    confidence = "High"
            else:
                if leaderType in ["Sec-dependent", "Double-glycine"]:
                    if matureP:
                        confidence = "High"
                    else:
                        confidence = "Low"
                elif leaderType == "Double-glycine or Sec-dependent":
                    if leaderSec:
                        confidence = "High"
                    else:
                        confidence = "Low"
                else:
                    confidence = "Low"
        else:
            confidence = "Low"
    else:
        confidence = "Low"
    return confidence


def predict(
    hmmexcute, ampep, Rcmd, hmmscan, pfamDomain, ncbiDomain, outDir, threshold, fasta
):
    identifier = Path(fasta).stem
    # ---------------------------- prediction based on domains
    current_path = os.path.split(os.path.realpath(__file__))[0]
    predict_d = os.path.join(current_path, "predict_based_on_domains.py")
    out_d = os.path.join(outDir, "prediction_domain")
    # Create folder
    out_hmmsearch = os.path.join(out_d, "hmmsearch")
    out_hmmsearch_parase = os.path.join(out_d, "hmmsearch_parase")

    cmd_d = f"python {predict_d} --hmmexcute {hmmexcute} --hmmscanExc {hmmscan} --pfamDomain {pfamDomain} --ncbiDomain {ncbiDomain} --inFasta {fasta} --outDir {out_d}"
    p = os.system(cmd_d)
    if p == 1:
        sys.exit("Failed to prediction based on domain rule!!!")
    result_d = os.path.join(out_d, "hmmsearch_parase", identifier + ".parase")

    # ---------------------------- prediction based on gene context
    predict_g = os.path.join(current_path, "predict_based_on_gene_context.py")
    out_g = os.path.join(outDir, "prediction_geneContext")
    out_summary = os.path.join(out_g, "summary")
    out_hmmsearch_gene = os.path.join(out_g, "hmm_geneContext")
    out_hmmsearch_parase_gene = os.path.join(out_g, "hmm_geneContext_parase")

    cmd_g = f"python {predict_g} --infasta {fasta} --ampep {ampep} --Rcmd {Rcmd} --hmmscan {hmmscan} --pfamDomain {pfamDomain} --ncbiDomain {ncbiDomain} --outDir {out_g}"
    q = os.system(cmd_g)
    if q == 1:
        sys.exit("Failed to prediction based on gene context rule!!!")
    result_g = os.path.join(out_g, "summary", identifier + ".result")

    # ---------------------------- combine result
    out_s = os.path.join(outDir, "results")
    if not os.path.exists(out_s):
        os.makedirs(out_s)
    combine_s = os.path.join(current_path, "combine.py")
    out_file = os.path.join(out_s, identifier + ".prediction")
    cmd_c = f"python {combine_s} --domainRes {result_d} --geneContextRes {result_g} --out {out_file}"
    os.system(cmd_c)

    # ---------------------------- annotating regions
    geneType_d = os.path.join(out_hmmsearch_parase, identifier + ".json")
    geneType_g = os.path.join(out_hmmsearch_parase_gene, identifier + ".json")
    f_d = open(geneType_d, "r")
    f_g = open(geneType_g, "r")
    dict_geneTypes_d = json.load(f_d)
    dict_geneTypes_g = json.load(f_g)
    f_d.close()
    f_g.close()
    dict_geneTypes = {**dict_geneTypes_d, **dict_geneTypes_g}
    annotate_code = os.path.join(current_path, "region_annotation.py")

    df_final_out = pd.read_table(out_file)
    final_gene_list = list(df_final_out.CDs)
    """
	cmd_a = f"python {annotate_code} --geneList {final_gene_list} --fasta {fasta} --outdir {outDir} \
				--out_prefix {identifier} --geneTypes {dict_geneTypes} --hmmscan {hmmscan} --pfamdomain {pfamDomain} --core {threshold}"
	os.system(cmd_a)
	"""
    annotate_regions(
        final_gene_list,
        fasta,
        outDir,
        identifier,
        dict_geneTypes,
        hmmscan,
        pfamDomain,
        threshold,
    )


def main():
    parse = argparse.ArgumentParser(
        prog="IIBacFinder", description="Detecting class II bacteriocins from genomes"
    )

    parse.add_argument(
        "-i",
        "--inDir",
        help="Input path to folder containg FASTA file, whose suffix could in list of ['.fas', '.fa', '.fasta', '.faa', '.fna']",
        required=True,
    )
    parse.add_argument(
        "-e",
        "--hmmexcute",
        help="The excutive hmmsearch, defualt: hmmsearch",
        default="hmmsearch",
        required=False,
    )
    parse.add_argument(
        "-a",
        "--ampep",
        help="The excutive ampep, defualt: ampep",
        default="ampep",
        required=False,
    )
    parse.add_argument(
        "-r",
        "--Rcmd",
        help="The excutive Rscript, defualt: Rscript",
        default="Rscript",
        required=False,
    )
    parse.add_argument(
        "-s",
        "--hmmscan",
        help="The excutive hmmscan, defualt: hmmscan",
        default="hmmscan",
        required=False,
    )
    # parse.add_argument("--pfamdomain", help="PFAM domain", required=True)
    parse.add_argument(
        "-o",
        "--outDir",
        help="The path to the folder storing output files",
        required=True,
    )
    parse.add_argument(
        "-t",
        "--threshold",
        help="Number of threshols used, default: 20",
        type=int,
        default=20,
        required=False,
    )
    # parse.add_argument("--signalp6", default=False, action="store_true", help="Whether perform signal peptide prediction using signalp6, default: FALSE")
    parse.add_argument(
        "-p",
        "--prodigal_short",
        default=True,
        action="store_false",
        help="Whether perform gene prediction using prodigal short (default: TRUE), toggle to close. This is indispensable for FASTA files.",
    )
    parse.add_argument(
        "-m",
        "--prodigal_p",
        default="single",
        choices=["single", "meta"],
        help="Select procedure model (single or meta) in 'prodigal_short', identical to the parameter '-p'.  Default is single",
    )
    parse.add_argument(
        "-v",
        "--version",
        action="version",
        version="IIBacFinder: v1.2",
        help="Print out the version and exit.",
    )

    args = parse.parse_args()

    ##############################################################
    # Creat folders
    ##############################################################
    current_path = os.path.split(os.path.realpath(__file__))[0]
    cp = Path(current_path)
    pfamdomain = os.path.join(cp.parent, "domains", "Pfam-A.hmm")
    ncbidomain = os.path.join(cp.parent, "domains", "hmm_PGAP.LIB")
    out_d = os.path.join(args.outDir, "prediction_domain")
    out_annotation = os.path.join(args.outDir, "region_annotation")
    out_plot = os.path.join(args.outDir, "region_plot")
    out_diamond = os.path.join(args.outDir, "diamond_alingment")
    # Create folder
    if not os.path.exists(out_d):
        os.makedirs(out_d)
    if not os.path.exists(out_annotation):
        os.makedirs(out_annotation)
    if not os.path.exists(out_plot):
        os.makedirs(out_plot)
    out_hmmsearch = os.path.join(out_d, "hmmsearch")
    if not os.path.exists(out_hmmsearch):
        os.makedirs(out_hmmsearch)
    out_hmmsearch_parase = os.path.join(out_d, "hmmsearch_parase")
    if not os.path.exists(out_hmmsearch_parase):
        os.makedirs(out_hmmsearch_parase)
    if not os.path.exists(out_diamond):
        os.makedirs(out_diamond)

    out_g = os.path.join(args.outDir, "prediction_geneContext")
    # Create folder
    if not os.path.exists(out_g):
        os.makedirs(out_g)
    out_summary = os.path.join(out_g, "summary")
    if not os.path.exists(out_summary):
        os.makedirs(out_summary)
    out_hmmsearch_gene = os.path.join(out_g, "hmm_geneContext")
    if not os.path.exists(out_hmmsearch_gene):
        os.makedirs(out_hmmsearch_gene)
    out_hmmsearch_parase_gene = os.path.join(out_g, "hmm_geneContext_parase")
    if not os.path.exists(out_hmmsearch_parase_gene):
        os.makedirs(out_hmmsearch_parase_gene)

    ##############################################################
    # Log file
    ##############################################################
    ### define the format of log
    log_file_path = os.path.join(args.outDir, "my.log")
    LOG_FORMAT = "[%(asctime)s][%(levelname)s][%(module)s] %(message)s"
    DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
    logging.basicConfig(
        filename=log_file_path,
        level=logging.DEBUG,
        format=LOG_FORMAT,
        datefmt=DATE_FORMAT,
        filemode="w",
    )

    ### start time
    start_time = time.time()
    start_time_out = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    logging.info("Starting at: " + start_time_out)

    ### Scan input files
    fasta_lists, genome_suffix = scanDir(args.inDir)
    logging.warning(f"Total fasta files for prediction: {len(fasta_lists)}")

    ##############################################################
    # Prediction
    ##############################################################
    # define partial function
    out_abs = os.path.abspath(args.outDir)
    partial_predict = functools.partial(
        predict,
        args.hmmexcute,
        args.ampep,
        args.Rcmd,
        args.hmmscan,
        pfamdomain,
        ncbidomain,
        out_abs,
        args.threshold,
    )

    if args.prodigal_short:
        prodigal_out = os.path.join(args.outDir, "prodigal_out")
        if not os.path.exists(prodigal_out):
            os.makedirs(prodigal_out)
        prodigal_exe = os.path.join(current_path, "prodigal-short")
        logging.info(f"Run prodigal-short for gene calling")
        partial_prodigal = functools.partial(
            run_prodigal, prodigal_exe, args.prodigal_p, prodigal_out
        )
        print("##############################")
        print("Run prodigal-short:")
        print("##############################")
        process_map(
            partial_prodigal, fasta_lists, max_workers=int(args.threshold), chunksize=1
        )
        # prediction list
        prediction_lists, fasta_suffix = scanDir(prodigal_out)
        print("##############################")
        print("Run class II bacteriocin prediction:")
        print("##############################")
        process_map(
            partial_predict,
            prediction_lists,
            max_workers=int(args.threshold),
            chunksize=1,
        )
    else:
        print("##############################")
        print("Run class II bacteriocin prediction:")
        print("##############################")
        process_map(
            partial_predict, fasta_lists, max_workers=int(args.threshold), chunksize=1
        )

    # Generate 'GenBank' annotation file
    if args.prodigal_short:
        for file in os.listdir(out_plot):
            if file.endswith("tsv"):
                genome_id = re.sub("_region.+", genome_suffix, file)
                genome_path = os.path.join(args.inDir, genome_id)
                outfile = re.sub("tsv", "gbk", file)
                read_region(
                    os.path.join(out_plot, file),
                    genome_path,
                    os.path.join(out_plot, outfile),
                )

    ##############################################################
    # Overall summary
    ##############################################################
    out_overall = os.path.join(args.outDir, "overall_result.tsv")
    out_s = os.path.join(args.outDir, "results")
    df_list = []
    for root, dirs, files in os.walk(out_s):
        for file in files:
            file_path = os.path.abspath(os.path.join(root, file))
            genome = Path(file).stem
            df = pd.read_table(file_path)
            df["Genome"] = genome
            if not df.empty:
                df_list.append(df)

    # check whether the predictions are empty
    if not df_list:
        logging.info("No putative class II bacteriocins have been identified!!!")
        sys.exit(0)

    result = pd.concat(df_list)

    ##############################################################
    # add region information
    ##############################################################
    out_region = os.path.join(args.outDir, "region_plot")
    result_region = region(out_region)
    df_combine = pd.merge(
        result,
        result_region,
        how="left",
        left_on=["Genome", "CDs"],
        right_on=["Genome", "CDs"],
    )
    df_combine["Uniq_ID"] = (
        df_combine["Genome"].astype(str) + "_+_" + df_combine["CDs"].astype(str)
    )
    logging.info(
        f"Total predicted precusors of class II bacteriocins: {df_combine.shape[0]}"
    )

    ##############################################################
    # Leader prediction
    ##############################################################
    # write out all precursors into FASTA file
    out_fa = os.path.join(args.outDir, "all.precusors.fa")
    out_sec_fa = os.path.join(args.outDir, "all.precusors.sec.fa")
    out_sec_json = os.path.join(args.outDir, "all.precusors.sec.json")
    out_gg_fa = os.path.join(args.outDir, "all.precusors.gg.fa")
    out_gg_json = os.path.join(args.outDir, "all.precusors.gg.json")

    json_dict_sec, json_dict_gg = {}, {}
    with open(out_fa, "w") as fo, open(out_sec_json, "w") as jse, open(
        out_gg_json, "w"
    ) as jgg, open(out_sec_fa, "w") as fos, open(out_gg_fa, "w") as fog:
        # fin.readline()
        # for line in fin:
        # 	genome = line.split("\t")[12]
        # 	cds = line.split("\t")[0]
        # 	seq = line.split("\t")[8]
        # 	print(">"+genome+"_+_"+cds, file=fo)
        # 	print(seq, file=fo)
        # 	json_dict[genome+"_+_"+cds] = seq
        for index, row in df_combine.iterrows():
            uniqid = row["Uniq_ID"]
            seq = row["Sequence"]
            print(">" + uniqid, file=fo)
            print(seq, file=fo)
            if row["Potential_Leader_Type"] == "Sec-dependent":
                json_dict_sec[uniqid] = seq
                print(">" + uniqid, file=fos)
                print(seq, file=fos)
            elif row["Potential_Leader_Type"] == "Double-glycine":
                json_dict_gg[uniqid] = seq
                print(">" + uniqid, file=fog)
                print(seq, file=fog)
            elif row["Potential_Leader_Type"] in [
                "Double-glycine or Sec-dependent",
                "Unclear",
            ]:
                json_dict_sec[uniqid] = seq
                print(">" + uniqid, file=fos)
                print(seq, file=fos)
                json_dict_gg[uniqid] = seq
                print(">" + uniqid, file=fog)
                print(seq, file=fog)
            else:
                pass
        json.dump(json_dict_sec, jse, indent=4)
        json.dump(json_dict_gg, jgg, indent=4)

    print("##############################")
    print("Run leader prediction:")
    print("##############################")
    logging.info("Run leader prediction")

    # perform prediction using deeep learning methods
    sigoutdir = os.path.join(args.outDir, "leader_prediction")
    if not os.path.exists(sigoutdir):
        os.makedirs(sigoutdir)

    # run prediction based on the premise that at least one file is not empty
    if os.path.getsize(out_sec_fa) > 0 or os.path.getsize(out_gg_fa) > 0:
        df_leader = leader(
            out_sec_fa, out_sec_json, out_gg_json, sigoutdir, args.threshold
        )
        df_combine_leader = pd.merge(
            df_combine, df_leader, on=["Uniq_ID", "Sequence"], how="left"
        )

        # Combine the leader prediction
        df_combine_leader["Predicted_mature_peptide"] = df_combine_leader.apply(
            lambda row: leader_combine(
                row["Potential_Leader_Type"],
                row["Sequence"],
                row["core_sec"],
                row["core_gg"],
            ),
            axis=1,
            result_type="expand",
        )
        # add confience level
        df_combine_leader["Confidence"] = df_combine_leader.apply(
            lambda row: confidenceP(
                row["RBS_motif"],
                row["Partial_index"],
                row["Rules"],
                row["Potential_Leader_Type"],
                row["Predicted_mature_peptide"],
                row["leader_sec"],
            ),
            axis=1,
            result_type="expand",
        )

        ##############################################################
        # Calculate physicochemical properties for core peptides
        ##############################################################
        df_p = df_combine_leader[["Uniq_ID", "Predicted_mature_peptide"]]
        df_p2 = df_p.set_index("Uniq_ID")
        df_p2 = df_p2[pd.notnull(df_p2["Predicted_mature_peptide"])]
        seq_lists = list(df_p2.to_dict()["Predicted_mature_peptide"].items())
        pro_r = process_map(calproperty, seq_lists, max_workers=args.threshold)
        pro_df = pd.DataFrame(pro_r)
        df_add_pro = pd.merge(df_combine_leader, pro_df, how="left", on="Uniq_ID")
        # df_add_pro.to_csv(out_overall, sep="\t", index=False)
    else:
        df_add_pro = df_combine

    ##############################################################
    # Compare candidates against known AMP (antimicrobial peptides) sequences
    ##############################################################
    out_alin = os.path.join(out_diamond, "alignment.tsv")
    amp_db = os.path.join(cp.parent, "AMP_database", "amp.index.dmnd")
    amp_annotation = os.path.join(cp.parent, "AMP_database", "annotation.tsv")
    df_amp = alnDiamond(out_fa, amp_db, out_alin, args.threshold, amp_annotation)
    if not df_amp.empty:
        df_add_amp = pd.merge(df_add_pro, df_amp, how="left", on=["Genome", "CDs"])
    else:
        df_add_amp = df_add_pro
    df_add_amp.to_csv(out_overall, sep="\t", index=False)

    # clean output folder
    logging.info("Purge temporary files")
    purge(".", "scored_features")

    ### end time
    end_time = time.time()
    end_time_out = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
    logging.info("Ending at: " + end_time_out)

    interval_time = end_time - start_time
    m, s = divmod(interval_time, 60)
    h, m = divmod(m, 60)
    logging.info("Total running time: %02d:%02d:%02d" % (h, m, s))


if __name__ == "__main__":
    main()
