import pandas as pd


def make_presence_absence_matrix(adjacency_file, outputfilename):
    dictionary_of_gene_presence_absence = {}
    with open(adjacency_file, "r") as fin:
        for line in fin:
            edge1 = line.split()[0]
            edge2 = line.split()[1]
            sample_1_name = edge1.rsplit('_', 1)[0]
            sample_2_name = edge2.rsplit('_', 1)[0]
            if edge1 not in dictionary_of_gene_presence_absence:
                dictionary_of_gene_presence_absence[edge1] = {}
            dictionary_of_gene_presence_absence[edge1][sample_1_name] = 1
            dictionary_of_gene_presence_absence[edge1][sample_2_name] = 1
    matrixDf = pd.DataFrame.from_dict(dictionary_of_gene_presence_absence)
    matrixDf = matrixDf.reindex(sorted(matrixDf.columns), axis=1)
    matrixDf = matrixDf.sort_index()
    matrixDf.fillna(0, inplace=True)
    matrixDf = matrixDf.T.astype(int)
    matrixDf.to_csv(outputfilename, sep='\t', index_label = "Gene")


make_presence_absence_matrix(snakemake.input.clusters, snakemake.output.matrix)
