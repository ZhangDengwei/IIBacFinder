import re
import argparse
from pysam import FastaFile

def tidy(input_file, out_prefix):
    out_combine = out_prefix + "_overall.tsv"
    out_sep = out_prefix + "_one_vs_one.tsv"
    out_pairs = out_prefix + "_pairs.tsv"  # New output file for pairs
    with open(input_file, 'r') as fin, open(out_combine, 'w') as fo1, open(out_sep, "w") as fo2, open(out_pairs, "w") as fo3:
        dic_cluster = dict()
        for line in fin:
            if re.search("^>", line):
                cluster_id = line.rstrip("\n")[1:]
                cluster_last = cluster_id
                dic_cluster[cluster_last] = [0, "*", "*", []] # gene number in cluster, representative gene, gene length, all genes in cluster
            else:
                dic_cluster[cluster_last][0] += 1
                gene_id = re.search(">(.+?)\.\.\.", line).group(1)
                if line.rstrip("\n").endswith("%"):
                    dic_cluster[cluster_last][3].append(gene_id)
                else:
                    dic_cluster[cluster_last][3].append(gene_id)
                    length = line.rstrip("\n").split()[1][:-1]
                    dic_cluster[cluster_last][1] = gene_id
                    dic_cluster[cluster_last][2] = length
        print("Cluster_ID\tGene_number\tRepresentative_Gene\tGene_Length\tGene_List", file=fo1)
        
        # Write header for pairs file
        print("Source\tTarget", file=fo3)
        
        for k, v in dic_cluster.items():
            cont = k + "\t" + str(v[0]) + "\t" + v[1] + "\t" + v[2] + "\t" + ",".join(v[3])
            print(cont, file=fo1, flush=True)
            for x in v[3]:
                cout_c = [x, k]
                print(*cout_c, sep="\t", file=fo2, flush=True)
            
            # Generate pairs: Representative_Gene (Source) vs all other genes (Target)
            rep_gene = v[1]  # Representative gene
            other_genes = [gene for gene in v[3] if gene != rep_gene]  # All genes except the representative
            
            for target in other_genes:
                print(f"{rep_gene}\t{target}", file=fo3, flush=True)
                

def main():
    parse = argparse.ArgumentParser(description="remodle the output of CD-hit")
    parse.add_argument("--input", help="the input file", required=True)
    parse.add_argument("--out_prefix", "-p", help="the prefix of output file", required=True)

    args = parse.parse_args()

    tidy(args.input, args.out_prefix)

if __name__ == "__main__":
    main()