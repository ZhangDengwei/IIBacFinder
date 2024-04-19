# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2024-04-02 20:18:16
# @Last Modified by:   chem
# @Last Modified time: 2024-04-03 13:01:03


import argparse
import pyfastx
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition


def read_region(region, fasta, outFile):
    df = pd.read_table(region)
    genome = df["Genome"][0]
    region_id = df["Region"][0]
    contig = df["Contig"][0]
    start_p = df["Start"].min()
    end_p = df["End"].max()

    # retrive the DNA sequence from fasta file
    fa = pyfastx.Fasta(fasta)
    interval = (int(start_p), int(end_p))
    sequence_string = fa.fetch(contig, interval)

    # create a sequence
    sequence_object = Seq(sequence_string)

    # Create a record
    record = SeqRecord(
        sequence_object,
        id=region_id,  # random accession number
        name=region_id,
        description="Putative biosynthetic gene cluster for class II bacteriocin",
    )
    record.annotations["molecule_type"] = "DNA"

    # Add annotation
    with open(region, "r") as fin:
        fin.readline()
        cds_record = []
        for line in fin:
            items = line.rstrip("\n").split("\t")
            cds_id = items[3]
            # whether the CDS is complete or truncated
            partial_index = items[9]
            # position
            strand = int(items[8])
            strand_re = 1 if strand == 1 else -1
            cds_start = (
                int(items[5]) - start_p + 1
                if strand_re == -1
                else int(items[5]) - start_p
            )
            cds_end = int(items[6]) - start_p + 1
            cds_seq = items[-1]
            domain_description = items[15]

            cds_class = items[-2]
            # gene_kind = "biosynthetic" if cds_class != "others" else "others"
            gene_kind = cds_class

            # add CDs annotation
            cds_feature = SeqFeature(
                FeatureLocation(
                    ExactPosition(cds_start), ExactPosition(cds_end), strand=strand_re
                ),
                type="CDS",
                qualifiers={
                    "locus_tag": cds_id,
                    "gene_kind": gene_kind,
                    "product": cds_class,
                    "partial_index": partial_index,
                },
            )
            cds_feature.qualifiers["transl_table"] = [11]
            cds_feature.qualifiers["translation"] = [cds_seq]
            if cds_id not in cds_record:
                record.features.append(cds_feature)
                cds_record.append(cds_id)

            if domain_description != "NA":
                domain_start = int(items[13]) + cds_start
                domain_end = int(items[14]) + cds_start
                pfam_id = items[12]
                pfam_des = items[15]
                domain_feature = SeqFeature(
                    FeatureLocation(
                        ExactPosition(domain_start), ExactPosition(domain_end)
                    ),
                    type="domain",
                    id=pfam_id,
                    qualifiers={"domain_name": pfam_id, "description": pfam_des},
                )
                record.features.append(domain_feature)

    # Save as GenBank file
    SeqIO.write(record, outFile, "genbank")
