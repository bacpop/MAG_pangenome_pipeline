import glob
from Bio import SeqIO
import os
# inputdir and outputdir need to be without trailing "/"


def fix_faa_files(list_of_annotation_dirs, outputdir):
    try:
        os.mkdir(outputdir)
    except OSError:
        pass
    for ann_dirs in list_of_annotation_dirs:
        file = list(set(glob.glob(ann_dirs + "/*.faa")) - set(glob.glob(ann_dirs + "/*.hypotheticals.faa")))[0] 
        filename = file.split("/")[-1]
        filename = filename.replace(".faa", "")
        records = SeqIO.parse(file, "fasta")
        outfilename = outputdir + "/" + filename + "_fixed.faa"
        all_records = []
        for record in records:
            desc = record.description.split("_", 1)[1]
            record.id = filename + "_" + desc
            all_records.append(record)
        SeqIO.write(all_records, outfilename, "fasta")


fix_faa_files(snakemake.input.annotations, snakemake.output.fixed_annotations)
