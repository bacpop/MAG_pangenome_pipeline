import glob

configfile: "config.yaml"

rule all:
    input:
        matrix = f"{config['output_dir']}/presence_absence_matrix.txt",
        summary_file = f"{config['output_dir']}/pangenome_summary.tsv",
        checkm_file = f"{config['output_dir']}/checkm_out.tsv",
        cgt_output = f"{config['output_dir']}/cgt_output.txt"

rule mkdir:
    output:
        ann_dir = directory(f"{config['output_dir']}/annotated")
    shell:
        """
        mkdir {output.ann_dir}
        """

rule bakta:
    input:
        genome = f"{config['genome_fasta']}/{{sample}}.fasta"
    output:
        ann_dir = directory(f"{config['output_dir']}/annotated/{{sample}}_ann")
    conda:
        "celebrimbor"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000
    params:
        DB = directory(f"{config['bakta_db']}")
    log:
        f"{config['output_dir']}/logs/bakta/{{sample}}.log"
    shell:
        """
        bakta {input.genome} --db {params.DB} --prefix {wildcards.sample} --translation-table 11 --skip-plot --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-crispr --skip-ori --threads {threads} --output {output.ann_dir} >{log} 2>&1
        """

def get_samples(genome_dir):
    list_of_samples = glob.glob(genome_dir + "/*.fasta")
    new_list_of_samples = []
    for sam in list_of_samples:
        sam = sam.split("/")[-1]
        sam = sam.replace(".fasta", "")
        new_list_of_samples.append(sam)
    return new_list_of_samples

if config['clustering_method'] in ["mmseqs2"]:
    rule fix_ffn_file:
        input:
            annotations = expand(f"{config['output_dir']}/annotated/{{sample}}_ann", sample=get_samples(config['genome_fasta']))
        output:
            fixed_annotations = directory(f"{config['output_dir']}/all_ffn")
        conda: #just needs biopython
            "celebrimbor"
        script: "scripts/fix_ffn_files.py"

    rule concat:
        input:
            ffn_files = f"{config['output_dir']}/all_ffn"
        output:
            f"{config['output_dir']}/all_samples.concat.ffn"
        shell:
            "cat {input.ffn_files}/* > {output}"

    rule mmseqs2:
        input:
            f"{config['output_dir']}/all_samples.concat.ffn"
        output:
            all_seqs = f"{config['output_dir']}/mmseqs/mmseqs_all_seqs.fasta",
            clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.tsv",
            rep_seq = f"{config['output_dir']}/mmseqs/mmseqs_rep_seq.fasta"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 16000 * attempt
        log:
            f"{config['output_dir']}/logs/mmseqs2.log"
        params:
            seq_id = 0.9,
            cov_mode = 0,
            c = 0.8,
            output_prefix = f"{config['output_dir']}/mmseqs/mmseqs",
            tmp_dir = f"{config['output_dir']}/mmseqs/tmp"
        conda:
            "celebrimbor"
        shell:
            "mmseqs easy-cluster {input} {params.output_prefix} {params.tmp_dir} --min-seq-id {params.seq_id} \
            --cov-mode {params.cov_mode} -c {params.c} --threads {threads} >{log} 2>&1"

    rule rep_seq_list:
        input:
            f"{config['output_dir']}/mmseqs/mmseqs_rep_seq.fasta"
        output:
            f"{config['output_dir']}/mmseqs/rep_sequences.list"
        shell:
            "grep '>' {input} | cut -f2 -d'>' | cut -f1 -d' ' > {output}"

    rule sort_mmseqs2_clusters:
        input:
            clusters = rules.mmseqs2.output.clusters
        output:
            sorted_clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.sorted.tsv"
        threads: 1
        resources:
            mem_mb=2000
        shell: "sort {input.clusters} > {output.sorted_clusters}"

    rule build_matrix_mmseqs:
        input:
            clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.sorted.tsv"
        output:
            matrix = f"{config['output_dir']}/presence_absence_matrix.txt"
        threads: 1
        resources:
            mem_mb=5000
        conda:
            # env has biopython and pandas
            "celebrimbor"
        script: "scripts/make_presence_absence_matrix.py"

    rule summarise_pangenome_mmseqs:
        input:
            rep_seq = f"{config['output_dir']}/mmseqs/mmseqs_rep_seq.fasta",
            matrix= f"{config['output_dir']}/presence_absence_matrix.txt"
        output:
            gene_descriptions = f"{config['output_dir']}/mmseqs/rep_seq_descriptions.txt",
            summary_file = f"{config['output_dir']}/pangenome_summary.tsv"
        threads: 1
        resources:
            mem_mb=5000
        params:
            breaks = config['cgt_breaks']
        conda:
            # env has biopython and pandas
            "celebrimbor"
        shell:
            """
            grep ">" {input.rep_seq} > {output.gene_descriptions}
            python scripts/summarise_pangenome_mmseqs.py {params.breaks} {output.gene_descriptions} {input.matrix} {output.summary_file}
            """

elif config['clustering_method'] in ["panaroo"]:
    rule symlink:
            input:
                indir = expand(f"{config['output_dir']}/annotated/{{sample}}_ann", sample=get_samples(config['genome_fasta']))
            output:
                outputdir = directory(f"{config['output_dir']}/all_gff")
            params:
                file_ext = "gff3"
            conda: #just needs biopython
                "celebrimbor"
            script: "scripts/create_symlink.py"

    rule panaroo:
        input:
            directory(f"{config['output_dir']}/all_gff")
        output:
            summary_stats = f"{config['output_dir']}/panaroo/summary_statistics.txt",
            rtab = f"{config['output_dir']}/panaroo/gene_presence_absence.Rtab",
            matrix = f"{config['output_dir']}/presence_absence_matrix.txt"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 16000 * attempt
        params:
            stringency = config['panaroo_stringency'],
            outdir = f"{config['output_dir']}/panaroo"
        conda:
            "celebrimbor"
        shell:
            """
            panaroo -i {input}/*.gff3 -o {params.outdir} --clean-mode {params.stringency} --remove-invalid-genes -t {threads}
            cp {output.rtab} {output.matrix}
            """

    rule summarise_pangenome_panaroo:
        input:
            matrix= f"{config['output_dir']}/presence_absence_matrix.txt"
        output:
            summary_file = f"{config['output_dir']}/pangenome_summary.tsv"
        threads: 1
        resources:
            mem_mb=5000
        params:
            breaks = config['cgt_breaks']
        conda:
            # env has biopython and pandas
            "celebrimbor"
        shell:
            """
            python scripts/summarise_pangenome_panaroo.py {params.breaks} {input.matrix} {output.summary_file}
            """


# checkm analysis
rule fix_faa_file:
        input:
            annotations = expand(f"{config['output_dir']}/annotated/{{sample}}_ann", sample=get_samples(config['genome_fasta']))
        output:
            fixed_annotations = directory(f"{config['output_dir']}/all_faa")
        conda: #just needs biopython
            "celebrimbor"
        script: "scripts/fix_faa_files.py"

if config['checkm_method'] in ["checkm1"]:
    rule run_checkm1:
        input:
            fixed_annotations = f"{config['output_dir']}/all_faa"
        output:
            workdir = directory(f"{config['output_dir']}/checkm1_out"),
            checkm_file = f"{config['output_dir']}/checkm_out.tsv"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 16000 * attempt
        shell:
            """
            checkm lineage_wf -q --genes -t {threads} -x faa --tab_table -f {output.checkm_file} {input.fixed_annotations} {output.workdir}
            sed 's/# //g' -i {output.checkm_file}
            sed 's/ /_/g' -i {output.checkm_file}
        """

elif config['checkm_method'] in ["checkm2"]:
    rule run_checkm2:
        input:
            fixed_annotations = directory(f"{config['output_dir']}/all_faa")
        output:
            workdir = directory(f"{config['output_dir']}/checkm2_out")
        params:
            DB = directory(f"{config['checkm2_db']}")
        threads: 1
        shell:
            """
            checkm2 predict --threads {threads} --genes -x .faa --database_path {params.DB} --input {input.fixed_annotations} --output-directory {output.workdir}
        """

# cgt analysis
rule run_cgt:
    input:
        matrix= f"{config['output_dir']}/presence_absence_matrix.txt",
        checkm_file = f"{config['output_dir']}/checkm_out.tsv"	
    output:
        cgt_output = f"{config['output_dir']}/cgt_output.txt"
    threads: 1
    params:
        exe = f"{config['cgt_exe']}",
        breaks = f"{config['cgt_breaks']}",
        error = f"{config['cgt_error']}"
    resources:
        mem_mb=5000
    shell:
        """
        {params.exe} --completeness-column 12 --breaks {params.breaks} --error {params.error} --output-file {output.cgt_output} {input.checkm_file} {input.matrix}
	    """
