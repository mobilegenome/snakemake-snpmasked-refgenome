rule strip_chr_prefix_from_fasta:
    input: config.get("genome"),
    output:
        fasta=temp("output/mm10_fixed_headers/genome.fa"),
        dir=directory("output/mm10_fixed_headers/")
    params:
        sed_expr="s/>.* />/g"
    message: "edit fasta headers"
    log:
        "logs/edit_fasta_header.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename fasta header' > {log}; "
        "sed '{params.sed_expr}' {input} > {output.fasta}"


rule strip_chr_prefix_from_gtf:
    input: config.get("annotation"),
    output:
        gtf=temp("output/mm10_fixed_headers.gtf"),
    params:
        sed_expr="s/^chr//g"
    message: "edit GTF headers"
    log:
        "logs/edit_gtf_names.log"
    shell:
        "echo 'running sed {params.sed_expr} to rename GTF header' > {log}; "
        "sed '{params.sed_expr}' {input} > {output.gtf}"