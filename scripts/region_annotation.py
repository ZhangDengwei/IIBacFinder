# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2023-05-10 11:58:26
# @Last Modified by:   chem
# @Last Modified time: 2024-04-03 20:55:31


import re
import os
import argparse
import subprocess
import fastaparser
from pathlib import Path
from predict_based_on_gene_context import surrounding, merge
from predict import purge


# retrive FASTA sequences of candidate genes
def fastaSeq(fasta, geneList, out):
    dict_fasta = dict()
    gene_for_annotation = []
    with open(fasta) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            sequence = seq.sequence_as_string()
            description = seq.description
            contig = "_".join(seq.id.split("_")[:-1])
            start_p = description.split("#")[1].strip()
            end_p = description.split("#")[2].strip()
            strand = description.split("#")[3].strip()
            partial_index = re.search("partial=(\d+);", description).group(1)
            if partial_index == "00":
                partial = "complete"
            elif partial_index == "01":
                partial = "right-incomplete"
            elif partial_index == "10":
                partial = "left-incomplete"
            elif partial_index == "11":
                partial = "both-incomplete"
            else:
                partial = "Code: " + partial_index

            start_type = re.search("start_type=(.+?);", description).group(1)
            rbs_motif = re.search("rbs_motif=(.+?);", description).group(1)
            dict_fasta[seq.id] = [
                sequence,
                contig,
                start_p,
                end_p,
                strand,
                partial,
                start_type,
                rbs_motif,
            ]
    # CDs in geneList may not be present when it exceed the contig boundary
    gene_final = list(set(geneList).intersection(set(list(dict_fasta.keys()))))
    gene_final_order = sorted(gene_final, key=lambda y: int(y.split("_")[-1]))

    with open(out, "w") as fo:
        for gene in gene_final_order:
            gene_for_annotation.append(gene)
            print(">" + gene, file=fo)
            print(dict_fasta[gene][0], file=fo)
    return (dict_fasta, gene_for_annotation)


"""
def parase_mmseqs(inPut):
	dict_mmseqs = {}
	with open(inPut, "r") as fin:
		for line in fin:
			query,target,evalue,pident,fident,nident,qlen,tlen,bits,theader = line.rstrip("\n").split("\t")
			annotation = " ".join(re.search("(.+?) n=\d+", theader).group(1).split()[1:])
			dict_mmseqs[query] = [target, evalue, annotation]
	return(dict_mmseqs)
"""


def overlap(region1, region2):
    # determine whether two regions are overlapped
    region1_start, region1_end = region1
    region2_start, region2_end = region2

    overlap_start = max(region1_start, region2_start)
    overlap_end = min(region1_end, region2_end)

    if overlap_end > overlap_start:
        overlap_length = overlap_end - overlap_start
    else:
        overlap_length = 0
    return overlap_length


# Excuting hmmscan
def hmmscan(hmmcmd, hmmdomain, outDir, inFasta, core):
    """
    hmmcmd: executive hmmscan, default "hmmscan"
    hmmdomain: PFAM domain
    inFasta: input FASTA file
    outDir: output dir
    """
    identifier = Path(inFasta).stem
    out_prediction = os.path.join(outDir, identifier + ".hmmscan")
    command = f"{hmmcmd} --domtblout {out_prediction} --cpu {core} --cut_ga {hmmdomain} {inFasta} > /dev/null"
    subprocess.run(command, shell=True)

    # Parase the result of hmmscan
    dict_hmmscan = dict()
    with open(out_prediction, "r") as fin:
        for line in fin:
            if not line.startswith("#"):
                items = re.split("\s+", line)
                target_name = items[0]
                target_acc = items[1]
                query = items[3]
                score = float(items[13])
                coord_from = int(items[17])
                coord_to = int(items[18])
                target_description = " ".join(items[22:])
                full_domain_annotation = (
                    f"{target_acc}: {target_description} ({target_name})"
                )
                if query not in dict_hmmscan:
                    dict_hmmscan[query] = []
                    dict_hmmscan[query].append(
                        [score, coord_from, coord_to, full_domain_annotation]
                    )
                else:
                    append = True
                    for x in dict_hmmscan[query]:
                        # determine whether the domain is overlaped
                        overlap_length = overlap([x[1], x[2]], [coord_from, coord_to])
                        if overlap_length > 10:
                            append = False
                            if score > x[0]:
                                dict_hmmscan[query].remove(x)
                                dict_hmmscan[query].append(
                                    [
                                        score,
                                        coord_from,
                                        coord_to,
                                        full_domain_annotation,
                                    ]
                                )
                                break
                            else:
                                pass
                    if append:
                        dict_hmmscan[query].append(
                            [score, coord_from, coord_to, full_domain_annotation]
                        )

    return dict_hmmscan


def annotate_regions(
    geneList, fasta, outDir, out_prefix, gene_types, hmmcmd, hmmdomain, core=10
):
    """
    geneList: Candidate genes of class II bacteriocins
    fasta: Fasta file of CDs
    """
    out_annotation = os.path.join(outDir, "region_annotation")
    out_plot = os.path.join(outDir, "region_plot")

    # Merge regions found
    regions_list = [surrounding(i) for i in geneList]
    merged_regions = merge(regions_list)

    # --- write out CDs within a region for annotation
    index = 1
    for region in merged_regions:
        out_fa = os.path.join(
            out_annotation, out_prefix + "_region_" + "{:0>3d}".format(index) + ".fa"
        )
        dict_fasta, gene_for_annotation = fastaSeq(fasta, region, out_fa)

        """
		#---run mmseq2 for protein annotation
		search_out = os.path.join(out_annotation, out_prefix+"_region_"+"{:0>3d}".format(index)+".mmseqs")
		tmp_dir = os.path.join(out_annotation, out_prefix+"temp")
		command = f"mmseqs easy-search {out_fa} {unirefdb} {search_out} {tmp_dir} --threads {core} --sort-results 1 --max-seqs 1 \
					--format-output query,target,evalue,pident,fident,nident,qlen,tlen,bits,theader > /dev/null"
		subprocess.run(command, shell=True)

		#--- parase mmseqs result
		dict_mmseqs = parase_mmseqs(search_out)
		"""
        # --- run hmmscan for protein annotation
        dict_hmmscan = hmmscan(hmmcmd, hmmdomain, out_annotation, out_fa, core)
        # print(dict_hmmscan)

        # --- write out gene annotaion
        out_annotation_file = os.path.join(
            out_plot, out_prefix + "_region_" + "{:0>3d}".format(index) + ".tsv"
        )
        with open(out_annotation_file, "w") as fo:
            headers = [
                "Genome",
                "Region",
                "Contig",
                "CDs",
                "Domain_index",
                "Start",
                "End",
                "Strand",
                "Orientation",
                "Partial_index",
                "Start_type",
                "RBS_motif",
                "Pfam_domain",
                "Domain_start",
                "Domain_end",
                "Domain_description",
                "CDs_Class",
                "Sequence",
            ]
            print(*headers, sep="\t", file=fo)

            # Some ±8 flanking CDs might be not present in FASTA as they exceed the contig region
            region_filter = list(set(region).intersection(set(list(dict_fasta.keys()))))
            region_order = sorted(region_filter, key=lambda y: int(y.split("_")[-1]))

            dict_domain_index = dict()
            domain_index = 1
            for gene in region_order:
                # print(gene)
                file_name = out_prefix + "_region_" + "{:0>3d}".format(index)
                locus_tag = gene
                (
                    sequence,
                    contig,
                    start_p,
                    end_p,
                    raw_orientation,
                    partial_index,
                    start_type,
                    rbs_motif,
                ) = dict_fasta[locus_tag]
                orientation = "1" if raw_orientation == "1" else 0
                strand = "forward" if orientation == "1" else "reverse"
                if gene not in dict_hmmscan:
                    pfam_domain, domain_start, domain_end, domain_description = (
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                    )
                    if gene in gene_types:
                        gene_class = gene_types[gene]
                    else:
                        gene_class = "others"
                    domain_index_out = "NA"
                    out_c = [
                        out_prefix,
                        file_name,
                        contig,
                        locus_tag,
                        domain_index_out,
                        start_p,
                        end_p,
                        strand,
                        orientation,
                        partial_index,
                        start_type,
                        rbs_motif,
                        pfam_domain,
                        domain_start,
                        domain_end,
                        domain_description,
                        gene_class,
                        sequence,
                    ]
                    print(*out_c, sep="\t", file=fo)
                else:
                    for domain in dict_hmmscan[gene]:
                        domain_start, domain_end, domain_full_description = (
                            domain[1],
                            domain[2],
                            domain[3],
                        )
                        pfam_domain = domain_full_description.split(":")[0]
                        if gene in gene_types:
                            gene_class = gene_types[gene]
                        else:
                            if re.search(r"regulator", domain_full_description, re.I):
                                gene_class = "regulator"
                            elif re.search(r"immunity", domain_full_description, re.I):
                                gene_class = "immunity"
                            elif re.search(
                                r"transporter", domain_full_description, re.I
                            ):
                                gene_class = "transporter"
                            elif re.search(
                                r"peptidase|protease", domain_full_description, re.I
                            ):
                                gene_class = "peptidase"
                            else:
                                gene_class = "others"
                        if pfam_domain not in dict_domain_index:
                            dict_domain_index[pfam_domain] = domain_index
                            domain_index_out = domain_index
                            domain_index += 1
                        else:
                            domain_index_out = dict_domain_index[pfam_domain]
                        out_c = [
                            out_prefix,
                            file_name,
                            contig,
                            locus_tag,
                            domain_index_out,
                            start_p,
                            end_p,
                            strand,
                            orientation,
                            partial_index,
                            start_type,
                            rbs_motif,
                            pfam_domain,
                            domain_start,
                            domain_end,
                            domain_full_description,
                            gene_class,
                            sequence,
                        ]
                        print(*out_c, sep="\t", file=fo)

        # --- plot region
        out_plot_file = os.path.join(
            out_plot, out_prefix + "_region_" + "{:0>3d}".format(index) + ".svg"
        )
        current_path = os.path.split(os.path.realpath(__file__))[0]
        script_path = os.path.join(current_path, "region_plot.R")
        command = (
            f"Rscript {script_path} {out_annotation_file} {out_plot_file} > /dev/null"
        )
        # subprocess.run(command, shell=True)
        proc = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        output, error = proc.communicate()

        index += 1

    # clean temp
    purge(out_annotation, "temp")


def main():
    parse = argparse.ArgumentParser(
        description="Summarizing, annotating, ploting regions of class II bacteriocins"
    )

    parse.add_argument(
        "--geneList", help="List of putative class II bacteriocin genes", required=True
    )
    parse.add_argument("--fasta", help="Input FASTA file for annotation", required=True)
    parse.add_argument("--outdir", help="Path to output", required=True)
    parse.add_argument("--out_prefix", help="Prefix of output files", required=True)
    parse.add_argument("--geneTypes", help="Gene types", required=True)
    parse.add_argument(
        "--hmmscan",
        help="The excutive hmmscan, defualt: hmmscan",
        default="hmmscan",
        required=False,
    )
    parse.add_argument("--pfamdomain", help="PFAM domain", required=True)
    parse.add_argument(
        "--core",
        help="Threshold used, default: 10",
        default=10,
        type=int,
        required=False,
    )

    args = parse.parse_args()

    annotate_regions(
        args.geneList,
        args.fasta,
        args.outDir,
        args.out_prefix,
        args.gene_types,
        args.hmmscan,
        args.pfamdomain,
        args.core,
    )


if __name__ == "__main__":
    main()
