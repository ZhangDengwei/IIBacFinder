# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2024-04-05 11:15:03
# @Last Modified by:   chem
# @Last Modified time: 2024-04-05 12:47:07


import os
import subprocess
import pandas as pd


def alnDiamond(infa, db, out, core, annotation):
    # excute diamond against AMP database
    command = f"diamond blastp --threads {core} --db {db} --ultra-sensitive --out {out} --outfmt 6 --query {infa} \
				--max-target-seqs 1 --max-target-seqs 1 --id 20 --query-cover 20 --masking 0 \
				--outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue full_qseq full_sseq  > /dev/null 2>&1"
    subprocess.run(command, shell=True)

    # parase result
    if os.path.getsize(out):
        df = pd.read_table(
            out,
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "slen",
                "qcovhsp",
                "scovhsp",
                "evalue",
                "full_qseq",
                "full_sseq",
            ],
        )
        df[["Genome", "CDs"]] = df["qseqid"].str.split("_\+_", expand=True)
        df["query_coverage"] = (
            df["qcovhsp"].astype(str)
            + " ("
            + df["length"].astype(str)
            + "/"
            + df["qlen"].astype(str)
            + ")"
        )
        df["target_coverage"] = (
            df["scovhsp"].astype(str)
            + " ("
            + df["length"].astype(str)
            + "/"
            + df["slen"].astype(str)
            + ")"
        )
        # merge annotation
        df_anno = pd.read_table(annotation)

        df_m = pd.merge(df, df_anno, how="left", left_on="sseqid", right_on="Accession")
        df_m_s = df_m[
            [
                "Genome",
                "CDs",
                "length",
                "pident",
                "query_coverage",
                "target_coverage",
                "full_sseq",
                "Accession",
                "Name",
                "Database",
            ]
        ]
        df_m_s = df_m_s.rename(
            {
                "pident": "Hsp_identity",
                "length": "Hsp_len",
                "query_coverage": "Coverage_q",
                "target_coverage": "Coverage_s",
                "full_sseq": "AMP_seq",
                "Accession": "AMP_accession",
                "Name": "AMP_name",
                "Database": "AMP_database",
            },
            axis=1,
        )
    else:
        df_m_s = pd.DataFrame()
    return df_m_s
