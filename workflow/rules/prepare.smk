rule strip_chr_prefix_from_fasta:
    """Remove chr prefix from FASTA
    
    This is needed to match the MGP VCF file. 
    """
    input: config.get("genome"),
    output:
        fasta=f"{OUTPUT_DIR}/reference_assembly/genome.no_chr_prefix.fa",
        fasta_generic=f"{OUTPUT_DIR}/reference_assembly/genome.fa",
        dir=directory(f"{OUTPUT_DIR}/reference_assembly/")
    params:
        sed_expr="s/>.* />/g"
    message: "edit fasta headers"
    log:
        "logs/edit_fasta_header.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename fasta header' > {log}; "
        "mkdir -p {output.dir}; "
        "sed '{params.sed_expr}' {input} > {output.fasta} "
        "ln -s {output.fasta} {output.fasta_generic}"

rule add_chr_prefix_to_fasta:
    input: config.get("genome"),
    output:
        fasta=temp(f"{OUTPUT_DIR}/reference_assembly/genome.chr_prefix.fa"),
        fasta_generic=f"{OUTPUT_DIR}/reference_assembly/genome.fa",
        dir=directory(f"{OUTPUT_DIR}/reference_assembly/")
    params:
        sed_expr="s/>([0-9MXY]{1,2})/>chr\1/g"
    message: "edit fasta headers"
    log:
        "logs/edit_fasta_header.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename fasta header' > {log}; "
        "mkdir -p {output.dir}; "
        "sed -E '{params.sed_expr}' {input} > {output.fasta}; "
        "ln -s {output.fasta} {output.fasta_generic}"



rule strip_chr_prefix_from_gtf:
    input: config.get("annotation"),
    output:
        gtf=temp(f"{OUTPUT_DIR}/reference_assembly/annotation.no_chr_prefix.gtf")
    params:
        sed_expr="s/^chr//g"
    message: "edit GTF headers"
    log:
        "logs/edit_gtf_names.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename GTF header' > {log}; "
        "mkdir -p {output.dir}; "
        "sed '{params.sed_expr}' {input} > {output.gtf}"

rule create_genome_file:
    input:
        fasta=config["genome"]
    output:
        "{OUTPUT_DIR}/reference_assembly/genome_chromosome_lengths.txt"
    envmodules:
        "samtools"
    params:
        input_filename=lambda wildcards, input: os.path.basename(input.fasta)
    shell:
        "cp {input.fasta} {params.input_filename} && "
        "samtools faidx {params.input_filename} && "
        "cut -f1,2 {params.input_filename} | > {output}"
