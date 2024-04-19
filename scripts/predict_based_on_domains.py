# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2023-04-20 20:30:09
# @Last Modified by:   chem
# @Last Modified time: 2024-04-12 21:43:40


# Import statements
import os
import re
import sys
import json
import argparse
import fastaparser
import subprocess
from pathlib import Path
import pandas as pd


# Excuting HMMsearch
def hmmExcute(HMMexcute, HMMdomain, output, query):
    """
    HMMexcute: executive hmmsearch
    HMMdomain: class II bacteriocin-related domains
    query: query FASTA file
    outDir: output dir
    """
    command = f"{HMMexcute} --cpu 5 --domtblout {output} -Z 61295632 -E 0.01 {HMMdomain} {query} > /dev/null"
    subprocess.run(command, shell=True)


# Summarizing the result of hmmsearch
def summarize(hmmResult):
    dict_domain = dict()

    with open(hmmResult, "r") as fin:
        for line in fin:
            if not line.startswith("#"):
                items = re.split("\s+", line)
                target_id = items[0]
                target_len = items[2]
                domain_id = items[3]
                evalue = items[6]
                hit_score = items[7]
                if int(target_len) < 150:
                    if target_id not in dict_domain:
                        dict_domain[target_id] = [domain_id, evalue, hit_score]
                    else:
                        if hit_score > dict_domain[target_id][2]:
                            dict_domain[target_id] = [domain_id, evalue, hit_score]
                        else:
                            pass
    return dict_domain


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


# retrive FASTA sequences of candidate genes
def fastaSeq(fasta, geneList, out):
    dict_fasta = dict()
    with open(fasta) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            dict_fasta[seq.id] = seq.sequence_as_string()
    # CDs in geneList may not be present when it exceed the contig boundary
    gene_final = list(set(geneList).intersection(set(list(dict_fasta.keys()))))

    with open(out, "w") as fo:
        for gene in gene_final:
            # sequence length should be < 150
            if len(dict_fasta[gene]) < 150:
                print(">" + gene, file=fo)
                print(dict_fasta[gene], file=fo)
    return dict_fasta


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
                items = re.split("\s+", line)
                domain_name = items[0]
                query = items[2]
                if query not in dict_hmmscan:
                    dict_hmmscan[query] = [domain_name]
                else:
                    dict_hmmscan[query].append(domain_name)

    return dict_hmmscan


# Performing domain rule-based prediction
def domainPredict(inFasta, HMMexcute, hmmscanExc, pfamDomain, ncbiDomain, outDir):
    dict_geneTypes = {}
    current_path = os.path.split(os.path.realpath(__file__))[0]
    cp = Path(current_path)
    HMMdomain = os.path.join(cp.parent, "domains", "classII-related.hmm")
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
    out_hmmsearch = os.path.join(outDir, "hmmsearch")
    out_search = os.path.join(out_hmmsearch, identifier + ".domtblout")
    # if not os.path.exists(out_hmmsearch):
    # os.makedirs(out_hmmsearch)
    hmmExcute(HMMexcute, HMMdomain, out_search, inFasta)

    # -------------2. Summarizing hmmsearch
    out_hmmsearch_parase = os.path.join(outDir, "hmmsearch_parase")
    # if not os.path.exists(out_hmmsearch_parase):
    # os.makedirs(out_hmmsearch_parase)
    out_parase_file = os.path.join(out_hmmsearch_parase, identifier + ".parase")

    dict_domain = summarize(out_search)

    if dict_domain:
        # -------------3. write out FASTA sequences
        fasta_out = os.path.join(out_hmmsearch_parase, identifier + ".fa")
        dict_fasta = fastaSeq(inFasta, list(dict_domain.keys()), fasta_out)

        # -------------4. HMMscan for domain annotation
        dict_hmmscan_pfam = hmmscan(
            hmmscanExc, pfamDomain, out_hmmsearch_parase, fasta_out, "pfam"
        )
        dict_hmmscan_ncbi = hmmscan(
            hmmscanExc, ncbiDomain, out_hmmsearch_parase, fasta_out, "ncbi"
        )

        # -------------5. Output results
        """
		dict_fasta = dict()
		with open(inFasta) as fasta_file:
			parser = fastaparser.Reader(fasta_file)
			for seq in parser:
				dict_fasta[seq.id] = seq.sequence_as_string()
	"""
    with open(out_parase_file, "w") as fo:
        header = [
            "CDs",
            "Annotation_Domain",
            "PFAM_domain",
            "NCBI_domain",
            "E_value",
            "Bit_Score",
            "Sequence",
            "Length",
        ]
        print(*header, sep="\t", file=fo)
        if dict_domain:
            for target, values in dict_domain.items():
                domain_id, evalue, hit_score = values
                pfam_domain = dict_hmmscan_pfam.get(target, "*")
                ncbi_domain = dict_hmmscan_ncbi.get(target, "*")
                both_domain = list(pfam_domain) + list(ncbi_domain)
                # Excluding the sequence containing domains not associated with class II bacteriocins
                if not list(
                    set(exclusion_domains).intersection(set(list(both_domain)))
                ):
                    out = [
                        target,
                        domain_id,
                        pfam_domain,
                        ncbi_domain,
                        evalue,
                        hit_score,
                    ]
                    seq = dict_fasta[target].rstrip("*")
                    # check whether the sequence contains any ambigious amino acids
                    if has_unambiguous_amino_acids(seq):
                        out += [seq, len(seq)]
                        print(*out, sep="\t", file=fo)
                        dict_geneTypes[target] = "precursor"
                    else:
                        pass

    # write out the gene types
    out_json = os.path.join(out_hmmsearch_parase, identifier + ".json")
    out_file = open(out_json, "w")
    json.dump(dict_geneTypes, out_file, indent=6)
    out_file.close()


def main():
    parse = argparse.ArgumentParser(
        description="Annotating class II bacteriocins with curated HMM domains"
    )

    parse.add_argument("--inFasta", help="The input FASTA files", required=True)
    parse.add_argument(
        "--hmmexcute",
        help="The excutive hmmsearch, defualt: hmmsearch",
        default="hmmsearch",
        required=False,
    )
    parse.add_argument(
        "--hmmscanExc",
        help="The excutive hmmscan, defualt: hmmscan",
        default="hmmscan",
        required=False,
    )
    parse.add_argument("--pfamDomain", help="The Pfam domain", required=True)
    parse.add_argument("--ncbiDomain", help="The NCBI domain", required=True)
    parse.add_argument(
        "--outDir", help="The path to the folder storing output files", required=True
    )

    args = parse.parse_args()

    domainPredict(
        args.inFasta,
        args.hmmexcute,
        args.hmmscanExc,
        args.pfamDomain,
        args.ncbiDomain,
        args.outDir,
    )


if __name__ == "__main__":
    main()
