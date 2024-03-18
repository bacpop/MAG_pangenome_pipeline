import glob
from Bio import SeqIO
import os
# inputdir and outputdir need to be without trailing "/"


def symlink(list_of_dirs, file_ext, outputdir):
    try:
        os.mkdir(outputdir)
    except OSError:
        pass
    for ann_dirs in list_of_dirs:
        file = glob.glob(ann_dirs + "/*." + file_ext)[0]
        os.symlink(file, outputdir)

symlink(snakemake.input.indir, snakemake.input.file_ext, snakemake.output.outputdir)
