# copy workflow/Makefile containing recipes to generate the required 
# CAROLI_EiJ SNP file into OUTPUT_DIR_MCAROLI and run it there 

rule copy_makefile:
    input:
        copy_from = "workflow/Makefile",
    output:
        copy_to = OUTPUT_DIR_MCAROLI + "/Makefile"
    shell:
        """
        cp {input.copy_from} {output.copy_to}
        """

rule run_makefile:
    input:
        rules.copy_makefile.output
    output:
        OUTPUT_DIR_MCAROLI + "/Caroli.snp",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.bed",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.mm10.bed",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.mm10.no_indels.bed",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.mm10.unmapped.bed",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.tar.gz",
        OUTPUT_DIR_MCAROLI + "/mm9ToMm10.over.chain.gz",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.mm10.no_indels.no_chr.bed",
        OUTPUT_DIR_MCAROLI + "/Caroli.snp.mm10.no_indels.no_chr.bed.sorted.gz",
        #OUTPUT_DIR_MCAROLI + "/Caroli.snp.no_indels.fa",
        #OUTPUT_DIR_MCAROLI + "/Caroli.snp.fa"
    conda:
        "../envs/makefile_mcaroli.yml"
    params:
        path = OUTPUT_DIR_MCAROLI
    shell:
        """
        make -C {params.path}
        """