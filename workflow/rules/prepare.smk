rule strip_chr_prefix_from_fasta:
    """Remove chr prefix from FASTA
    
    This is needed to match the MGP VCF file. 
    """
    input: config.get("genome"),
    output:
        fasta=temp(f"{OUTPUT_DIR}/mm10_fixed_headers/genome.fa"),
        dir=directory(f"{OUTPUT_DIR}/mm10_fixed_headers/")
    params:
        sed_expr="s/>.* />/g"
    message: "edit fasta headers"
    log:
        "logs/edit_fasta_header.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename fasta header' > {log}; "
        "sed '{params.sed_expr}' {input} > {output.fasta}"

rule add_chr_prefix_to_fasta:
    input: config.get("genome"),
    output:
        fasta=temp(f"{OUTPUT_DIR}/mm10_fixed_headers/genome_chr.fa"),
        dir=directory(f"{OUTPUT_DIR}/mm10_fixed_headers/")
    params:
        sed_expr="s/>([0-9MXY]{1,2})/>chr\1/g"
    message: "edit fasta headers"
    log:
        "logs/edit_fasta_header.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename fasta header' > {log}; "
        "sed -E '{params.sed_expr}' {input} > {output.fasta}"


rule strip_chr_prefix_from_gtf:
    input: config.get("annotation"),
    output:
        gtf=temp(f"{OUTPUT_DIR}/mm10_fixed_headers.gtf")
    params:
        sed_expr="s/^chr//g"
    message: "edit GTF headers"
    log:
        "logs/edit_gtf_names.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename GTF header' > {log}; "
        "sed '{params.sed_expr}' {input} > {output.gtf}"