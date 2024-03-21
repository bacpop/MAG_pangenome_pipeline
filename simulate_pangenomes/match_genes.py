import argparse


def get_options():
    description = 'Matches genes pre and post simulation using remove_sequences.py'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python match_genes.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--input',
                    default=None,
                    required=True,
                    help='Output .tsv from mmeseqs2. Should be run on concatenated mmseqs2 all_seqs.fasta files from pre and post simulation using remove_sequence.py.')
    IO.add_argument('--output',
                    default="matched_genes",
                    help='Output prefix. Default = "matched_genes"')

    return parser.parse_args()

def main():
    options = get_options()
    input = options.input
    output = options.output

    #testing
    #input = "data/mmseqs_comp_cluster.tsv"
    #output = "test_gene_matching"

    # first entry in tuple is pre-sim gene, second is post-sim gene
    gene_mappings = []

    pre_sim_genes = {}
    post_sim_genes = set()

    with open(input, "r") as f:
        cluster = None
        for line in f.readlines():
            split_line = line.rstrip().split("\t")
            cluster_name = split_line[0]
            gene_name = split_line[1]

            if cluster != cluster_name:
                if cluster != None:
                    # new cluster found
                    for entry in post_sim_genes:
                        genome_name = entry.split("_comp_")[0]
                        if genome_name in pre_sim_genes:
                            if len(pre_sim_genes[genome_name]) > 0:
                                match = pre_sim_genes[genome_name].pop()
                                gene_mappings.append((match, entry))
                            else:
                                gene_mappings.append((None, entry))
                        else:
                            gene_mappings.append((None, entry))
                    # get genes that are no longer matched pre-sim
                    for genome, gene_set in pre_sim_genes.items():
                        for entry in gene_set:
                            gene_mappings.append((entry, None))
                    # reset for new cluster
                    pre_sim_genes = {}
                    post_sim_genes = set()
                cluster = cluster_name
            
            # add gene names 
            if "_comp_" in gene_name:
                post_sim_genes.add(gene_name)
            else:
                genome_name = gene_name.rsplit('_', 1)[0]
                if genome_name not in pre_sim_genes:
                    pre_sim_genes[genome_name] = set()
                pre_sim_genes[genome_name].add(gene_name)
    
    # output completeness file
    with open(output + ".tsv", "w+") as f:
        f.write("Pre-sim_gene\tPost-sim_gene\n")
        for entry in gene_mappings:
            f.write(str(entry[0]) + "\t" + str(entry[1]) + "\n")
            

if __name__ == "__main__":
    main()