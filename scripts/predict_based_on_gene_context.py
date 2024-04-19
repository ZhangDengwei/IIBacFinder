# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2023-04-21 10:51:49
# @Last Modified by:   chem
# @Last Modified time: 2024-04-09 20:40:57


# Import statements
import os
import re
import sys
import json
import argparse
import subprocess
import fastaparser
import pandas as pd
import networkx as nx
from pathlib import Path


# Excuting HMMsearch
def hmmExcute(HMMexcute, HMMdomain, outDir, query):
    """
    HMMexcute: executive hmmsearch
    HMMdomain: class II bacteriocin-related domains
    query: query FASTA file
    outDir: output dir
    """
    identifier = Path(query).stem
    out_domain = os.path.join(outDir, identifier + ".tblout")
    command = f"{HMMexcute} --cpu 5 --tblout {out_domain} -Z 61295632 -E 0.001 {HMMdomain} {query} > /dev/null"
    subprocess.run(command, shell=True)


# Access the gene id of ±8 CDs
def surrounding(cdsID):
    # cdsID must be assigned by Prodigal
    contig = "_".join(cdsID.split("_")[:-1])
    gene_order = cdsID.split("_")[-1]
    # Including ±8 CDs
    down_order = int(gene_order) - 8 if int(gene_order) - 8 > 0 else 1
    up_order = int(gene_order) + 8
    gene_list = ["_".join([contig, str(i)]) for i in range(down_order, up_order + 1)]
    return gene_list


# Merge the regions found
def merge(inList):
    G = nx.Graph()
    G.add_nodes_from(sum(inList, []))
    q = [[(s[i], s[i + 1]) for i in range(len(s) - 1)] for s in inList]
    for i in q:
        G.add_edges_from(i)
    merge_list = [list(i) for i in nx.connected_components(G)]
    merge_list_ordered = [
        sorted(x, key=lambda y: int(y.split("_")[-1])) for x in merge_list
    ]
    return merge_list_ordered


# determine ambiguous amino acids
def has_unambiguous_amino_acids(sequenceIN):
    normal_aa = [
        "Z",
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    seq_list = list(set(list(sequenceIN)))
    result = set(seq_list).issubset(normal_aa)
    return result


# Summarizing the result of hmmsearch
def summarize(hmmResult, outDir):
    dict_domain = dict()
    dict_geneTypes = {}
    identifier = Path(hmmResult).stem
    out_file = os.path.join(outDir, identifier + "_hmmsearch.out")

    inclusion_rules = [
        "Peptidase_C39",
        "EntA_Immun",
        "DUF5841",
        "bPH_2",
        "DUF1430",
        "Rce1-like",
        "DUF5976",
        "DUF5836",
        "imm_Aci221A",
        "imm_Aci221B",
        "imm_BM1I",
        "imm_ORF16",
        "imm_ORF2",
        "imm_OrfB",
        "imm_abpIM",
        "imm_bacB_32",
        "imm_bacB_43",
        "imm_bimlp",
        "imm_bimlp2",
        "imm_brcI",
        "imm_cbiA",
        "imm_cbiB2",
        "imm_dviA",
        "imm_entCI",
        "imm_entiP",
        "imm_entqC",
        "imm_gakI",
        "imm_lafI",
        "imm_lagC",
        "imm_lciM",
        "imm_lebI",
        "imm_papB",
        "imm_plnc8C",
        "imm_saiA",
        "imm_sakTIM",
        "imm_spiA",
        "imm_spiQ",
    ]
    protease_related = ["Peptidase_C39"]
    immunity_related = [
        "EntA_Immun",
        "DUF1430",
        "Rce1-like",
        "imm_Aci221A",
        "imm_Aci221B",
        "imm_BM1I",
        "imm_ORF16",
        "imm_ORF2",
        "imm_OrfB",
        "imm_abpIM",
        "imm_bacB_32",
        "imm_bacB_43",
        "imm_bimlp",
        "imm_bimlp2",
        "imm_brcI",
        "imm_cbiA",
        "imm_cbiB2",
        "imm_dviA",
        "imm_entCI",
        "imm_entiP",
        "imm_entqC",
        "imm_gakI",
        "imm_lafI",
        "imm_lagC",
        "imm_lciM",
        "imm_lebI",
        "imm_papB",
        "imm_plnc8C",
        "imm_saiA",
        "imm_sakTIM",
        "imm_spiA",
        "imm_spiQ",
    ]
    other_related = ["DUF5841", "bPH_2", "DUF5976", "DUF5836"]
    transporter = ["ABC_tran", "MFS_1"]
    typeI_lantibiotic = "LANC_like"

    transporter_cds = []
    lantibiotic_cds = []
    """
	#----rule 1:
	1. Peptidase_C39 present;
	2. transporter present;
	3. Lanthipeptides-related domains (LANC_like.hmm) absent;
	#----rule 2:
	1. Other related domains present;
	2. transporter present;
	"""
    rule1, rule2 = [], []
    hits_genes = []
    hit_rules = dict()

    with open(hmmResult, "r") as fin:
        for line in fin:
            if not line.startswith("#"):
                items = re.split(r"\s+", line)
                target_id = items[0]
                domain_id = items[2]
                evalue = items[7]
                hit_score = items[8]

                if domain_id in protease_related:
                    rule1.append(target_id)
                    hit_rules[target_id] = "Peptidase_C39"
                    hits_genes.append(target_id)
                    if target_id not in dict_geneTypes:
                        dict_geneTypes[target_id] = ["peptidase"]
                    else:
                        dict_geneTypes[target_id].append("peptidase")
                elif domain_id in immunity_related:
                    rule2.append(target_id)
                    # hit_rules[target_id] = "Immunity_related"
                    hit_rules[target_id] = domain_id
                    hits_genes.append(target_id)
                    if target_id not in dict_geneTypes:
                        dict_geneTypes[target_id] = ["immunity"]
                    else:
                        dict_geneTypes[target_id].append("immunity")
                elif domain_id in other_related:
                    rule2.append(target_id)
                    hit_rules[target_id] = "Others_related"
                    hits_genes.append(target_id)
                    if target_id not in dict_geneTypes:
                        dict_geneTypes[target_id] = ["others"]
                    else:
                        dict_geneTypes[target_id].append("others")
                elif domain_id in transporter:
                    transporter_cds.append(target_id)
                    if target_id not in dict_geneTypes:
                        dict_geneTypes[target_id] = ["transporter"]
                    else:
                        dict_geneTypes[target_id].append("transporter")
                elif domain_id == "LANC_like":
                    lantibiotic_cds.append(target_id)

                if target_id not in dict_domain:
                    dict_domain[target_id] = [[domain_id, evalue, hit_score]]
                else:
                    dict_domain[target_id].append([domain_id, evalue, hit_score])

    # re-formate dict_geneTypes
    dict_geneTypes_re = dict()
    for gene, types in dict_geneTypes.items():
        if len(types) == 1:
            dict_geneTypes_re[gene] = types[0]
        else:
            dict_geneTypes_re[gene] = "/".join(sorted(types))

    # Determine whether the gene context encoding potential novel class II bacteriocins
    gene_candidates = []
    gene_rules = dict()

    # Merge regions found
    regions_list = [surrounding(i) for i in hits_genes]
    merged_region = merge(regions_list)

    # Filter lists
    for region in merged_region:
        including_rule = [hit_rules.get(gene) for gene in region]
        including_rules = list(filter(lambda item: item is not None, including_rule))
        # 1. Case 1: "Peptidase_C39" domain present
        if "Peptidase_C39" in including_rules:
            rec_1_1 = list(set(region).intersection(set(lantibiotic_cds)))
            if not rec_1_1:
                # Check whether "ABC_tran" is present
                rec_1_2 = list(set(region).intersection(set(transporter_cds)))
                if rec_1_2:
                    for geneid in region:
                        gene_rules[geneid] = list(set(including_rules))
                    gene_candidates.extend(region)
        else:
            # 2. Case 2: Other domain present
            # other_four_domain = ["EntA_Immun", "DUF5841", "bPH_2", "DUF1430"]
            # intect = list(set(including_rules).intersection(set(other_four_domain)))
            intect = list(set(including_rules))
            if intect:
                # examine the presence of transporter
                rec_2 = list(set(region).intersection(set(transporter_cds)))
                if rec_2:
                    for geneid in region:
                        gene_rules[geneid] = list(set(including_rules))
                    gene_candidates.extend(region)

    with open(out_file, "w") as fo:
        for target, values in dict_domain.items():
            out = [target] + values
            print(*out, sep="\t", file=fo)

    return (list(set(gene_candidates)), gene_rules, dict_geneTypes_re)


# retrive FASTA sequences of candidate genes
def fastaSeq(fasta, geneList, out):
    dict_fasta = dict()
    gene_for_prediction = []
    with open(fasta) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            dict_fasta[seq.id] = seq.sequence_as_string()
    # CDs in geneList may not be present when it exceed the contig boundary
    gene_final = list(set(geneList).intersection(set(list(dict_fasta.keys()))))

    with open(out, "w") as fo:
        for gene in gene_final:
            # sequence length should be < 150
            gene_seq = dict_fasta[gene].rstrip("*")
            if (len(gene_seq) < 150) and (has_unambiguous_amino_acids(gene_seq)):
                gene_for_prediction.append(gene)
                print(">" + gene, file=fo)
                print(gene_seq, file=fo)
    return (dict_fasta, gene_for_prediction)


# Excuting amPEP
def ampepExcute(amPEP, amPEPmodel, outDir, inFasta):
    """
    amPEP: executive amPEP, default "ampep"
    HMMdomain: class II bacteriocin-related domains
    inFasta: input FASTA file
    outDir: output dir
    """
    identifier = Path(inFasta).stem
    out_prediction = os.path.join(outDir, identifier + ".ampep")
    command = f"{amPEP} predict -t 5 -m {amPEPmodel} -i {inFasta} -o {out_prediction} --seed 2012 > /dev/null"
    # subprocess.run(command, shell=True)
    proc = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = proc.communicate()
    # Retrive the prediction result
    df = pd.read_table(out_prediction)
    # Sequences with probability_AMP> 0.8 are considered AMP
    ampep_AMP = list(df[(df.predicted == "AMP") & (df.probability_AMP > 0.8)].seq_id)
    return ampep_AMP


# Excuting ampir
def ampirExcute(Rcmd, outDir, inFasta):
    """
    Rcmd: executive Rscript, default "Rscript"
    inFasta: input FASTA file
    outDir: output dir
    """
    identifier = Path(inFasta).stem
    current_path = os.path.split(os.path.realpath(__file__))[0]
    script_path = os.path.join(current_path, "ampir.R")
    out_prediction = os.path.join(outDir, identifier + ".ampir")
    command = f"{Rcmd} {script_path} {inFasta} {out_prediction} > /dev/null"
    subprocess.run(command, shell=True)
    # Retrive the prediction result
    df = pd.read_table(out_prediction)
    # Sequences with probability_AMP> 0.8 are considered AMP
    ampir_AMP = list(df[df.prob_AMP > 0.8].seq_name)
    return ampir_AMP


# Excuting AmpGram
# Ref: https://github.com/michbur/AmpGram
def ampgramExcute(Rcmd, outDir, inFasta):
    # Rcmd: executive Rscript, default "Rscript"
    # inFasta: input FASTA file
    # outDir: output dir

    identifier = Path(inFasta).stem
    current_path = os.path.split(os.path.realpath(__file__))[0]
    script_path = os.path.join(current_path, "AmpGram.R")
    out_prediction = os.path.join(outDir, identifier + ".ampgram")
    command = f"{Rcmd} {script_path} {inFasta} {out_prediction} > /dev/null"
    # subprocess.run(command, shell=True)
    proc = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = proc.communicate()
    # Retrive the prediction result
    df = pd.read_table(out_prediction)
    # Sequences with probability_AMP> 0.8 are considered AMP
    ampgram_AMP = list(df[df.prob_AMP > 0.8].seq_name)
    return ampgram_AMP


# Excuting APIN
# Ref: https://github.com/zhanglabNKU/APIN
def apinExcute(currentPath, outDir, inFasta):
    """
    amPEP: executive amPEP, default "ampep"
    HMMdomain: class II bacteriocin-related domains
    inFasta: input FASTA file
    outDir: output dir
    """
    identifier = Path(inFasta).stem
    out_prediction = os.path.join(outDir, identifier + ".apin")
    # command = f"{amPEP} predict -t 5 -m {amPEPmodel} -i {inFasta} -o {out_prediction} --seed 2012 > /dev/null"
    script_path = os.path.join(currentPath, "APIN/proposed_fusion_model.py")
    command = f"python {script_path} -epochs 5 -test_file {inFasta} -prediction_file {out_prediction} > /dev/null"
    # subprocess.run(command, shell=True)
    proc = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = proc.communicate()
    # Retrive the prediction result
    # df = pd.read_table(out_prediction)
    seq_id = []
    apin_AMP = []
    with open(inFasta, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                seq_id.append(line.rstrip("\n").lstrip(">"))
    index = 0
    with open(out_prediction, "r") as fin:
        for line in fin:
            if line.rstrip("\n") == "1":
                apin_AMP.append(seq_id[index])
            else:
                pass
            index += 1
    return apin_AMP


# Excuting hmmscan
def hmmscan(hmmcmd, hmmdomain, outDir, inFasta, database):
    """
    hmmcmd: executive hmmscan, default "hmmscan"
    hmmdomain: PFAM domain
    inFasta: input FASTA file
    outDir: output dir
    """
    identifier = Path(inFasta).stem
    out_put = identifier + f".{database}.hmmscan"
    out_prediction = os.path.join(outDir, out_put)
    command = (
        f"{hmmcmd} --tblout {out_prediction} --cut_ga {hmmdomain} {inFasta} > /dev/null"
    )
    subprocess.run(command, shell=True)

    # Parase the result of hmmscan
    dict_hmmscan = dict()
    with open(out_prediction, "r") as fin:
        for line in fin:
            if not line.startswith("#"):
                items = re.split(r"\s+", line)
                domain_name = items[0]
                query = items[2]
                if query not in dict_hmmscan:
                    dict_hmmscan[query] = [domain_name]
                else:
                    dict_hmmscan[query].append(domain_name)

    return dict_hmmscan


# Performing gene context rule-based prediction
def geneContextPredict(
    inFasta, HMMexcute, amPEP, Rcmd, hmmcmd, pfamDomain, ncbiDomain, outDir
):
    """
    inFasta: input FASTA
    outDir: output dir
    """
    current_path = os.path.split(os.path.realpath(__file__))[0]
    cp = Path(current_path)
    HMMdomain = os.path.join(cp.parent, "domains", "rule.gene.context.hmm")
    amPEPmodel = os.path.join(cp.parent, "models", "amPEP.model")
    pfam_filter_pfam = os.path.join(cp.parent, "filter", "filter_pfam_35.tsv")
    pfam_filter_ncbi = os.path.join(cp.parent, "filter", "filter_ncbi.tsv")

    # exclude PFAM domains
    df_pfam = pd.read_table(pfam_filter_pfam)
    df_ncbi = pd.read_table(pfam_filter_ncbi)
    exclusion_domains_pfam = df_pfam.loc[df_pfam.Related == "No", "#Name"].to_list()
    exclusion_domains_ncbi = df_ncbi.loc[df_ncbi.Related == "No", "#Name"].to_list()
    exclusion_domains = exclusion_domains_pfam + exclusion_domains_ncbi

    # -------------1. Excute hmmsearch against query FASTA
    identifier = Path(inFasta).stem
    out_hmmsearch = os.path.join(outDir, "hmm_geneContext")
    hmmExcute(HMMexcute, HMMdomain, out_hmmsearch, inFasta)
    output_hmmsearch = os.path.join(outDir, "hmm_geneContext", identifier + ".tblout")

    # -------------2. Parase hmmseach output
    out_hmmsearch_parase = os.path.join(outDir, "hmm_geneContext_parase")
    gene_candidates, gene_rules, dict_geneTypes = summarize(
        output_hmmsearch, out_hmmsearch_parase
    )

    out_summary = os.path.join(outDir, "summary")
    summary_final = os.path.join(out_summary, identifier + ".result")
    with open(summary_final, "w") as fo:
        header = [
            "CDs",
            "Prediction_rules",
            "PFAM_domain",
            "NCBI_domain",
            "Sequence",
            "Length",
        ]
        print(*header, sep="\t", file=fo)
        # Candidates may not be found based on gene context rules
        if gene_candidates:
            # -------------3. Retrive FASTA sequences of candidate genes
            out_candidates = os.path.join(
                outDir, "hmm_geneContext_parase", identifier + ".fa"
            )
            fasta_seq, gene_prediction = fastaSeq(
                inFasta, gene_candidates, out_candidates
            )
            # Length of flanking genes may be >150aa
            if gene_prediction:
                # -------------4. Excuting amPEP
                ampep_AMP_list = ampepExcute(
                    amPEP, amPEPmodel, out_hmmsearch_parase, out_candidates
                )

                # -------------5. Excuting ampir
                ampir_AMP_list = ampirExcute(Rcmd, out_hmmsearch_parase, out_candidates)

                # -------------6. Excuting APIN
                # apin_AMP_list = apinExcute(current_path, out_hmmsearch_parase, out_candidates)
                # ampgram_AMP_list = ampgramExcute(Rcmd, out_hmmsearch_parase, out_candidates)

                # overall_AMP = list(set(ampep_AMP_list + ampir_AMP_list + ampgram_AMP_list))
                overall_AMP = list(set(ampep_AMP_list + ampir_AMP_list))

                # -------------7. Run hmmscan to annotate the PFAM domians of candidates
                dict_hmmscan_pfam = hmmscan(
                    hmmcmd, pfamDomain, out_hmmsearch, out_candidates, "pfam"
                )
                dict_hmmscan_ncbi = hmmscan(
                    hmmcmd, ncbiDomain, out_hmmsearch, out_candidates, "ncbi"
                )

                for classII_candidate in overall_AMP:
                    predict_rule = sorted(gene_rules[classII_candidate])
                    pfam_domain = dict_hmmscan_pfam.get(classII_candidate, "*")
                    ncbi_domain = dict_hmmscan_ncbi.get(classII_candidate, "*")
                    both_domain = list(pfam_domain) + list(ncbi_domain)
                    # Excluding the sequence containing domains not associated with class II bacteriocins
                    if not list(set(exclusion_domains).intersection(set(both_domain))):
                        seq = fasta_seq[classII_candidate].rstrip("*")
                        print(
                            *[
                                classII_candidate,
                                predict_rule,
                                pfam_domain,
                                ncbi_domain,
                                seq,
                                len(seq),
                            ],
                            sep="\t",
                            file=fo,
                        )

                        dict_geneTypes[classII_candidate] = "precursor"

    # -------------3. Output gene types
    out_json = os.path.join(out_hmmsearch_parase, identifier + ".json")
    out_file = open(out_json, "w")
    json.dump(dict_geneTypes, out_file, indent=6)
    out_file.close()


def main():
    parse = argparse.ArgumentParser(
        description="Annotating class II bacteriocins based on gene context rules"
    )

    parse.add_argument(
        "--infasta", help="Input FASTA file for annotation", required=True
    )
    parse.add_argument(
        "--hmmexcute",
        help="The excutive hmmsearch, defualt: hmmsearch",
        default="hmmsearch",
        required=False,
    )
    parse.add_argument(
        "--ampep",
        help="The excutive ampep, defualt: ampep",
        default="ampep",
        required=False,
    )
    parse.add_argument(
        "--Rcmd",
        help="The excutive Rscript, defualt: Rscript",
        default="Rscript",
        required=False,
    )
    parse.add_argument(
        "--hmmscan",
        help="The excutive hmmscan, defualt: hmmscan",
        default="hmmscan",
        required=False,
    )
    parse.add_argument("--pfamDomain", help="PFAM domain", required=True)
    parse.add_argument("--ncbiDomain", help="NCBI domain", required=True)
    parse.add_argument("--outDir", help="Output dir", required=True)

    args = parse.parse_args()

    geneContextPredict(
        args.infasta,
        args.hmmexcute,
        args.ampep,
        args.Rcmd,
        args.hmmscan,
        args.pfamDomain,
        args.ncbiDomain,
        args.outDir,
    )


if __name__ == "__main__":
    main()
