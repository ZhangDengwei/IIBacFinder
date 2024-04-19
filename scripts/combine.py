# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2023-04-21 17:12:03
# @Last Modified by:   chem
# @Last Modified time: 2024-04-06 14:05:17


import os
import argparse
import pandas as pd
from pathlib import Path


def combine(domainRes, geneContextRes, out):
    # Parease annotation file
    current_path = os.path.split(os.path.realpath(__file__))[0]
    cp = Path(current_path)
    description = os.path.join(cp.parent, "domains", "curated_domain_description.tsv")
    dict_desp = dict()
    with open(description, "r") as fin:
        fin.readline()
        for line in fin:
            items = line.rstrip("\n").split("\t")
            domain_name = items[0]
            domain_description = items[4]
            leader_type = items[5]
            dict_desp[domain_name] = [domain_description, leader_type]

    # Combine results based on two rules
    df_d = pd.read_table(domainRes)
    df_g = pd.read_table(geneContextRes)
    gene_list_d = list(df_d.CDs)
    gene_list_g = list(df_g.CDs)

    # genes found by both rule
    list_both = list(set(gene_list_d).intersection(set(gene_list_g)))
    # genes found exculsively by domian rule
    list_d = list(set(gene_list_d).difference(set(gene_list_g)))
    # genes found exculsively by gene context rule
    list_g = list(set(gene_list_g).difference(set(gene_list_d)))

    with open(out, "w") as fo:
        header = [
            "CDs",
            "Rules",
            "Domain_rule",
            "Domain_Evalue",
            "Domain_Bitscore",
            "Context_rule",
            "PFAM_domain",
            "NCBI_domain",
            "Sequence",
            "Length",
            "Description",
            "Potential_Leader_Type",
        ]
        print(*header, sep="\t", file=fo)
        if list_both:
            for gene in list_both:
                domain_rule = df_d.loc[df_d.CDs == gene, "Annotation_Domain"].values[0]
                domain_evalue = df_d.loc[df_d.CDs == gene, "E_value"].values[0]
                domain_bitscore = df_d.loc[df_d.CDs == gene, "Bit_Score"].values[0]
                context_rule = df_g.loc[df_g.CDs == gene, "Prediction_rules"].values[0]
                pfam_domain = df_g.loc[df_g.CDs == gene, "PFAM_domain"].values[0]
                ncbi_domain = df_g.loc[df_g.CDs == gene, "NCBI_domain"].values[0]
                seq = df_d.loc[df_d.CDs == gene, "Sequence"].values[0]
                seq_len = df_d.loc[df_d.CDs == gene, "Length"].values[0]
                des_domain, leader_pre = dict_desp[domain_rule]
                out_c = [
                    gene,
                    "Both",
                    domain_rule,
                    domain_evalue,
                    domain_bitscore,
                    context_rule,
                    pfam_domain,
                    ncbi_domain,
                    seq,
                    seq_len,
                    des_domain,
                    leader_pre,
                ]
                print(*out_c, sep="\t", file=fo)
        if list_d:
            for gene in list_d:
                domain_rule = df_d.loc[df_d.CDs == gene, "Annotation_Domain"].values[0]
                pfam_domain = df_d.loc[df_d.CDs == gene, "PFAM_domain"].values[0]
                ncbi_domain = df_d.loc[df_d.CDs == gene, "NCBI_domain"].values[0]
                domain_evalue = df_d.loc[df_d.CDs == gene, "E_value"].values[0]
                domain_bitscore = df_d.loc[df_d.CDs == gene, "Bit_Score"].values[0]
                seq = df_d.loc[df_d.CDs == gene, "Sequence"].values[0]
                seq_len = df_d.loc[df_d.CDs == gene, "Length"].values[0]
                des_domain, leader_pre = dict_desp[domain_rule]
                out_c = [
                    gene,
                    "Domain",
                    domain_rule,
                    domain_evalue,
                    domain_bitscore,
                    "*",
                    pfam_domain,
                    ncbi_domain,
                    seq,
                    seq_len,
                    des_domain,
                    leader_pre,
                ]
                print(*out_c, sep="\t", file=fo)
        if list_g:
            for gene in list_g:
                context_rule = df_g.loc[df_g.CDs == gene, "Prediction_rules"].values[0]
                pfam_domain = df_g.loc[df_g.CDs == gene, "PFAM_domain"].values[0]
                ncbi_domain = df_g.loc[df_g.CDs == gene, "NCBI_domain"].values[0]
                seq = df_g.loc[df_g.CDs == gene, "Sequence"].values[0]
                seq_len = df_g.loc[df_g.CDs == gene, "Length"].values[0]
                if "Peptidase_C39" in context_rule:
                    des_rule = "Peptidase_C39 related"
                    leader_pre = "Double-glycine"
                else:
                    if "EntA_Immun" in context_rule:
                        des_rule = "Enterocin A immunity related"
                        leader_pre = "Double-glycine or Sec-dependent"
                    elif "DUF5841" in context_rule:
                        des_rule = "Enterocin A biosynthesis related"
                        leader_pre = "Double-glycine or Sec-dependent"
                    elif "bPH_2" in context_rule:
                        des_rule = "Potentially involve in the biosynthesis of leaderless peptides"
                        leader_pre = "Leaderless"
                    elif "DUF1430" in context_rule:
                        des_rule = (
                            "Association with lactococcin 972 family of bacteriocins"
                        )
                        leader_pre = "Sec-dependent"
                    elif "DUF5976" in context_rule:
                        des_rule = "Association with salivaricin CRL1328"
                        leader_pre = "Double-glycine"
                    elif "DUF5836" in context_rule:
                        des_rule = "Association with ABP-118"
                        leader_pre = "Double-glycine"
                    else:
                        des_rule = "Other immunity related"
                        leader_pre = "Unclear"

                out_c = [
                    gene,
                    "GeneContext",
                    "*",
                    "*",
                    "*",
                    context_rule,
                    pfam_domain,
                    ncbi_domain,
                    seq,
                    seq_len,
                    des_rule,
                    leader_pre,
                ]
                print(*out_c, sep="\t", file=fo)

    putative_classII_genes = list_both + list_d + list_g
    return putative_classII_genes


def main():
    parse = argparse.ArgumentParser(
        description="Annotating class II bacteriocins with curated HMM domains"
    )

    parse.add_argument(
        "--domainRes", help="The prdiction results based on domain rule", required=True
    )
    parse.add_argument(
        "--geneContextRes",
        help="The prdiction results based on gene context rule",
        required=True,
    )
    parse.add_argument("--out", help="The combined prdiction results", required=True)

    args = parse.parse_args()

    combine(args.domainRes, args.geneContextRes, args.out)


if __name__ == "__main__":
    main()
