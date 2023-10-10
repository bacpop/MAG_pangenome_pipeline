import pandas as pd


def make_presence_absence_matrix(adjacency_file, representative_sequences, outputfilename):
    reference_list = []
    with open(representative_sequences, "r") as reps_in:
        for name in reps_in:
            name = name.rstrip()
            reference_list.append(name)
    dictionary_of_gene_presence_absence = {}
    for rep in reference_list:
        dictionary_of_gene_presence_absence[rep] = {}
    with open(adjacency_file, "r") as fin:
        for line in fin:
            edge1 = line.split()[0]
            edge2 = line.split()[1]
            sample_1_name = edge1.split("_")[0]
            sample_2_name = edge2.split("_")[0]
            dictionary_of_gene_presence_absence[edge1][sample_1_name] = 1
            dictionary_of_gene_presence_absence[edge1][sample_2_name] = 1
    matrixDf = pd.DataFrame.from_dict(dictionary_of_gene_presence_absence)
    matrixDf = matrixDf.reindex(sorted(matrixDf.columns), axis=1)
    matrixDf = matrixDf.sort_index()
    matrixDf.fillna(0, inplace=True)
    matrixDf = matrixDf.T.astype(int)
    matrixDf.to_csv(outputfilename, sep='\t', index_label = "Gene")


make_presence_absence_matrix(snakemake.input.clusters, snakemake.input.rep_list, snakemake.output.matrix)
