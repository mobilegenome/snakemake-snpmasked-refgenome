import os
import pathlib
import tempfile


checkpoint extract_snv_with_snp_split:
    """
    Extract SNVs from the Mouse Genome Project (MGP) VCF file 
    for strains listed in configured `strains` using SNPsplit's 
    genome_preparation tool
    
    SNPsplit: https://github.com/FelixKrueger/SNPsplit
    
    """
    input:
        vcf_file="data/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",
        ref_genome=rules.strip_chr_prefix_from_fasta.output.dir,
    output:
        dir_snp=temp(directory(f"{OUTPUT_DIR}/SNPsplit/{{strain}}/SNPs_{{strain}}/")),
        dir_fa=directory(f"{OUTPUT_DIR}/SNPsplit/{{strain}}/{{strain}}_N-masked/"),
        report=f"{OUTPUT_DIR}/SNPsplit/{{strain}}/{{strain}}_SNP_filtering_report.txt",
        archive=f"{OUTPUT_DIR}/SNPsplit/{{strain}}/all_SNPs_{{strain}}_GRCm38.txt.gz",
    message:
        "snpsplit prepare"
    conda:
        "../envs/snpsplit.yml"
    log:
        f"{OUTPUT_DIR}/SNPsplit/snp_split_{{strain}}.log",
    params:
        strain=lambda wildcards: wildcards.strain,
    script:
        "../scripts/snpsplit.py"


rule snp_split_create_sorted_bed:
    """
    Convert extracted SNV from from SNPsplit to a gzipped BED file.
    
    
    """
    input:
        lambda wildcards: checkpoints.extract_snv_with_snp_split.get(**wildcards).output.archive
    output:
        bed=f"{OUTPUT_DIR}/SNPsplit/{{strain}}/all_SNPs_{{strain}}_GRCm38.bed.gz"
    log:
        logs="logs/create_bed_{strain}.log"
    shell:
        "zcat {input} | sort -V -k2,2 -k3,3 | awk 'OFS=\"\t\" {{print $2, $3-1, $3, $5}}' | gzip -c > {output.bed} 2> {log}"


def add_additional_variants(*labels):
    add_bed_files = config["additional_variant_bed_files"]
    paths = []
    for label in labels:
        if label in add_bed_files:
            paths.append(add_bed_files[label])
        else:
            raise KeyError(f"Key {label} not found in config: {add_bed_files.keys()}")
    if paths:
        return paths
    else:
        raise ValueError

rule merge_bed_files:
    """Merge multipe BED files
    
    Add "chr" prefix
    """
    input:
        expand(f"{{output_dir}}/SNPsplit/{{strain}}/all_SNPs_{{strain}}_GRCm38.bed.gz",
            strain=STRAINS,
            output_dir=[OUTPUT_DIR]
        ) +
        add_additional_variants("Mus_caroli")
    output:
        bed=OUTPUT_DIR + "/merged/all_SNPs_all_strains_GRCm38.bed.gz",
    params:
        add_chr="| sed 's/^/chr/g' | "
    resources:
        mem_mb=64000
    envmodules:
        "bedtools/2.24.0"
    shell:
        "zcat {input} | "
        "sed 's#/#\t#g' | "
        "bedtools sort -i - | "
        "bedtools merge -i - -c 4,5 -o collapse "
        "{params.add_chr} "
        " gzip -c > {output.bed}"


rule intersection:
    input:
        bed_files=expand(f"{{output_dir}}/SNPsplit/{{strain}}/all_SNPs_{{strain}}_GRCm38.bed.gz",
            strain=STRAINS,
            output_dir=[OUTPUT_DIR]
        ) +
        add_additional_variants("Mus_caroli"),
        genome_file=rules.create_genome_file.output
    output:
        OUTPUT_DIR + "/merged/all_SNPs_all_strains_GRCm38.intersect.bed.gz",
    params:
        names="CAST SPRET CAROLI"
    envmodules:
        "bedtools/2.24.0"
    shell:
        "bedtools multiinter -i {input.bed_files} -empty -g {input.genome_file} -names {params.names} | "
        "gzip -c > {output}"

rule maskfasta:
    """Mask coordinates from BED file in FASTA"""
    input:
        bed=f"{OUTPUT_DIR}/merged/all_SNPs_all_strains_GRCm38.bed.gz",
        fasta=config["genome"]
    output:
        fasta=f"{OUTPUT_DIR}/merged/GRCm38_masked_all_strains.fa"
    resources:
        mem_mb=64000
    envmodules:
        "bedtools/2.24.0"
    shell:
        "bedtools maskfasta -bed {input.bed} -fi {input.fasta} -fo {output.fasta}"

if PER_STRAIN_FASTA:
    rule snp_split_concat:
        """Concatenate per-chromsome/scaffold FASTA files from SNPsplit 
        to a unified genome-wide FASTA.
        """
        input:
            lambda wildcards: sorted(
                glob.glob(
                    "{directory}/{filename_pattern}".format(
                        directory=checkpoints.extract_snv_with_snp_split.get(**wildcards).output.dir_fa,
                        filename_pattern="*.fa",
                    )
                )
            ),
        wildcard_constraints:
            strain="|".join(config.get("strains")),
        output:
            f"{OUTPUT_DIR}/GRCm38_masked_{{strain}}.fa",
        run:
            print("{input}")
            shell("cat {input} > {output}")
