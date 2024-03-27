
FASTA_MM10 = /omics/odcf/analysis/OE0538_projects/DO-0008/data/snp_masked_genome/input_files/cellranger-ref/mm10/refdata-gex-mm10-2020-A/fasta/genome.fa
FASTA_CAROLI = Caroli.snp.fa
FASTA_CAROLI_NOINDELS = Caroli.snp.no_indels.fa

all: Caroli.snp.mm10.bed Caroli.snp.mm10.unmapped.bed $(FASTA_CAROLI) $(FASTA_CAROLI_NOINDELS) Caroli.snp.tar.gz mm9ToMm10.over.chain.gz

#all: Caroli.snp.tar.gz

#download: Caroli.snp.tar.gz mm9ToMm10.over.chain.gz
#	wget ftp://ftp.ebi.ac.uk/pub/databases/vertebrategenomics/FOG09/snp/Caroli.snp.tar.gz
#	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz

# download
Caroli.snp.tar.gz:
	wget ftp://ftp.ebi.ac.uk/pub/databases/vertebrategenomics/FOG09/snp/Caroli.snp.tar.gz

mm9ToMm10.over.chain.gz:
	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz


Caroli.snp: Caroli.snp.tar.gz
	tar xzf $<

Caroli.snp.bed: Caroli.snp
	awk -F "\t" '{print "chr"$$1, $$2-1, $$2, $$3, $$4}' $< > $@

Caroli.snp.mm10.bed: Caroli.snp.bed mm9ToMm10.over.chain.gz
	liftOver $^ $@ Caroli.snp.mm10.unmapped.bed

Caroli.snp.mm10.no_indels.bed: Caroli.snp.mm10.bed
	grep -v "\*" $< > $@

$(FASTA_CAROLI): Caroli.snp.mm10.bed
	module load bedtools;
	bedtools maskfasta -fi $(FASTA_MM10) -fo $@ -bed $<

$(FASTA_CAROLI_NOINDELS): Caroli.snp.mm10.no_indels.bed
	module load bedtools;
	bedtools maskfasta -fi $(FASTA_MM10) -fo $@ -bed $<


