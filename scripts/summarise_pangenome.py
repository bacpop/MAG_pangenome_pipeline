"""
script to summarise pangenomes
"""

import sys
import pandas as pd

breaks = sys.argv[1]
rare_perc = float(breaks.split(",")[0])
core_perc = float(breaks.split(",")[1])
gene_desc = sys.argv[2]
matrix = sys.argv[3]
outputfile = sys.argv[4]

# import matrix as pandas df
gene_matrixDf = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
n_MAG = len(gene_matrixDf.columns)
#gene_matrixDf = gene_matrixDf.T
#genes = (list(gene_matrixDf)[0:])

# put gene descriptions in dictionary:
gene_description_dictionary = {}
with open(gene_desc, "r") as desc_in:
    for line in desc_in:
        gene_name = line.split()[0]
        gene_name = gene_name.replace(">", "")
        description = line.split()[1:]
        description = "_".join(description)
        gene_description_dictionary[gene_name] = description

for gene, row in gene_matrixDf.iterrows():
    # sum number of MAGs the gene appears in:
    count = row.sum()
    prevalence = float(count) / float(n_MAG)
    pan_type = "rare"
    if prevalence >= core_perc:
        pan_type = "core"
    else:
        if prevalence >= rare_perc:
            pan_type = "middle"
    with open(outputfile, "a") as fout:
        fout.write("%s\t%s\t%s\t%s\n" % (gene, gene_description_dictionary[gene], prevalence, pan_type))


