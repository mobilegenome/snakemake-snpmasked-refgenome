# SNP-masked refererence genome

Create a reference genome sequence with "injected" SNPs/SNVs from other mouse strains.

![Rule graph](rulegraph.png)


## Input files

- mouse reference genome
- mouse annotation (GTF)
- SNP file

## Dependencies

The `snpsplit` tool is defined in a conda environment in `workflow/envs/snpslit.yml`.

*CellRanger* is proprietary software and cannot be installed via Conda. It needs to be manually downloaded and the path
to its executable stored in `workflow/config.yaml`.

On the DKFZ Compute Cluster Cellranger is installed as module and available at:
`/software/cellranger/6.0.0/bin/cellranger mkref`

## Configuration

The configuration in `workflow/config.yaml` contains:

- path to reference genome
- path to gene annotation
- path to `cellranger mkref`
- list of Strains from the mouse genome project (MGP)

### Sequence names

There is no consensus whether to use `chr` prefixes for sequence/chromosome names in 
FASTA and annotation files. 

In the results, we aim to use `chr` prefixes matching the convention used by the cellranger 
reference dataset. 

This matrix gives an overview what to expect and requires modification:

| Mod? | Dataset / File                    | prefix | example | contents                  |
|----- |-----------------------------------|--------|---------|--------------------------|
| Yes  | `mgp_REL2005_snps_indels.vcf.gz`  | None   | `1`     | 1-19,X missing: Y, MT |
| Yes  | `Caroli.snp`                      | None   | `1`     | 1-19
| No   | `Caroli.snp.bed`                  | `chr`  | `chr1`  | chr1-19,X,Y,M  + unlocalized | 
| No   | cellranger reference 20202-A      | `chr`  | `chr1`  | chr1-19,X,Y,M  + unlocalized | 
| No   | `gencode.vM25.annotation.gtf `    | `chr`  | `chr1`  | chr1-19,X,Y,M  |
| No   | `mm9ToMm10.over.chain.gz`         | `chr`  | `chr1`  | chr1-19,X,Y,M  |



## Run

The workflow can be run by executing:

```bash
snakemake -c1 --configfile config.yaml --use-conda
```

To use the LSF job scheduler on the DKFZ cluster run with

```bash
snakemake --cluster "bsub -n4 -q verylong -R rusage[mem=100GB]" -p -j2 -c4 --configfile config.yaml --use-conda
```

## Output files

Output files of the workflow are stored in the subdirectory `output/`

# TO DO

- Add a download function for the SNP data
  from `ftp.sanger.ac.uk/pub/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`
  This is currently not activated due to proxy issues in the DKFZ working.
