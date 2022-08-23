rule cellranger_rna_mkref_merged:
    """Create Cellranger Reference for a unified N-masked genome.
    """
    input:
        fasta=rules.maskfasta.output.fasta,
        annotation=config.get("annotation"),
    output:
        directory("output/GRCm38_masked_allStrains"),
    params:
        mem=300,
        output_root_dir=lambda wildcards, output: os.path.split(output[0])[0],
        cellranger_mkref_bin=config["executables"]["cellranger-rna"],
        genome="GRCm38_masked_allStrains",
    log:
        "cellranger_rna_mkref_GRCm38_masked_allStrains.log",
    shell:
        "cd output && "
        "{params.cellranger_mkref_bin} "
        "--genome={params.genome} "
        "--fasta=../{input.fasta} "
        "--genes=../{input.annotation} "
        "--memgb={params.mem} > {log}"

if PER_STRAIN_FASTA:
    rule cellranger_rna_mkref:
        input:
            fasta=rules.snp_split_concat.output,
            annotation=rules.strip_chr_prefix_from_gtf.output.gtf,
        output:
            directory("output/GRCm38_masked_{strain}"),
        wildcard_constraints:
            strain="|".join(config.get("strains")),
        params:
            mem=300,
            output_root_dir=lambda wildcards, output: os.path.split(output[0])[0],
            cellranger_mkref_bin=config["executables"]["cellranger-rna"],
            genome=lambda wildcards, input: f"GRCm38_masked_{wildcards.strain}",
        log:
            "cellranger_rna_mkref_GRCm38_masked_{strain}.log",
        shell:
            "cd output && "
            "{params.cellranger_mkref_bin} "
            "--genome={params.genome} "
            "--fasta=../{input.fasta[0]} "
            "--genes=../{input.annotation} "
            "--memgb={params.mem} > {log}"