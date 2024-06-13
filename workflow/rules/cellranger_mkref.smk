from pathlib import Path

if config["cellranger_filter_gtf"]:
    rule cellranger_rna_extract_version_suffix_from_annotation:
        input:
            config["annotation"] 
        output:
            f"{OUTPUT_DIR}/reference_assembly/genes.modified.gtf"
        log:
            "logs/cellranger_rna_extract_version_suffix_from_annotation.log",
        params:
            regex_pattern = "(gene|transcript|exon)_id \"(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)\"",
            repl_pattern = "\\1_id \\2; \\1_version \\4"
        script:
            "../scripts/cellranger_modify_gtf.py"

if config["cellranger_filter_gtf"]:
    rule cellranger_rna_filter_gtf:
        input:
            gtf=rules.cellranger_rna_extract_version_suffix_from_annotation.output
        output:
            f"{OUTPUT_DIR}/reference_assembly/genes.filtered.gtf"
        log:
            "logs/cellranger_rna_filter_gtf.log",

        params:
            accepted_gene_types=" ".join(
                map(lambda val: f"--attribute=gene_type:{val}",
                    config["cellranger_accepted_gene_transcript_types"])
            ),
            accepted_transcript_types=" ".join(
                map(lambda val: f"--attribute=transcript_type:{val}",
                    config["cellranger_accepted_gene_transcript_types"])
            ),
        envmodules:
            "cellranger/6.1.1"
        shell:
            "cellranger mkgtf "
            "{input.gtf} "
            "{output} "
            "{params.accepted_gene_types} "
            "{params.accepted_transcript_types} "
            "> {log}"


def cellranger_mkref_get_gtf_input():
    if config["cellranger_filter_gtf"]:
        return rules.cellranger_rna_filter_gtf.output
    else:
        return config["annotation"]

rule cellranger_rna_mkref_merged:
    """Create Cellranger Reference for a unified N-masked genome.
    """

    input:
        fasta=rules.maskfasta.output.fasta,
        annotation=cellranger_mkref_get_gtf_input()
    output:
        directory(f"{OUTPUT_DIR}/merged/GRCm38_masked_allStrains"),
    params:
        mem=300,
        genome=lambda wildcards, output: Path(output[0]).parts[-1] #"GRCm38_masked_allStrains", 
    envmodules:
        "cellranger/6.1.1"
    log:
        "cellranger_rna_mkref_GRCm38_masked_allStrains.log",
    script:
        "../scripts/cellranger_mkref.py"

# This rule in mode = "maskfasta" with condition per_strain_fasta = True
# only generates output for CAST and SPRET, not CARO

# In mode = "incorporate_snvs", output of this rule is not included in targets,
# so no cellranger reference is generated from the SNV-injected fastas
# it would also only apply to CAST and SPRET, not CARO
if PER_STRAIN_FASTA:
    rule cellranger_rna_mkref: 
        input:
            fasta=rules.snp_split_concat.output,
            annotation=rules.strip_chr_prefix_from_gtf.output.gtf,
        output:
            directory(f"{OUTPUT_DIR}/{{strain}}/GRCm38_masked_{{strain}}"),
        wildcard_constraints:
            strain="|".join(config.get("strains")),
        params:
            mem=300,
            output_root_dir=lambda wildcards, output: os.path.split(output[0])[0],
            cellranger_mkref_bin=config["executables"]["cellranger-rna"], 
            genome=lambda wildcards, input: f"GRCm38_masked_{wildcards.strain}",
        envmodules:
            "cellranger/6.1.1"
        log:
            "logs/cellranger_rna_mkref_GRCm38_masked_{strain}.log",
        shell:
            "cd output && "
            "cellranger mkref "
            "--genome={params.genome} "
            "--fasta=../{input.fasta[0]} "
            "--genes=../{input.annotation} "
            "--memgb={params.mem} > {log}"
