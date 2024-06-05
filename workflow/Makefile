
#FASTA_MM10 = /omics/odcf/analysis/OE0538_projects/DO-0008/data/snp_masked_genome/input_files/cellranger-ref/mm10/refdata-gex-mm10-2020-A/fasta/genome.fa
#FASTA_CAROLI = Caroli.snp.fa
#FASTA_CAROLI_NOINDELS = Caroli.snp.no_indels.fa

all: Caroli.snp.mm10.bed Caroli.snp.mm10.unmapped.bed Caroli.snp.tar.gz mm9ToMm10.over.chain.gz Caroli.snp.mm10.no_indels.no_chr.bed.sorted.gz Caroli.snp.mm10.no_indels.no_chr.bed 
#$(FASTA_CAROLI) $(FASTA_CAROLI_NOINDELS)

# download Caroli SNP file
Caroli.snp.tar.gz:
	wget ftp://ftp.ebi.ac.uk/pub/databases/vertebrategenomics/FOG09/snp/Caroli.snp.tar.gz

# download liftover chain file for converting mm9 to mm10 coordinates
mm9ToMm10.over.chain.gz:
	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz

# extracting
Caroli.snp: Caroli.snp.tar.gz
	tar xzf $<

# adding "chr", turning start = pos-1, end = pos
# to adjust the format from VCF (1-based) to BED (0-based) 
Caroli.snp.bed: Caroli.snp
	awk -F "\t" '{print "chr"$$1, $$2-1, $$2, $$3, $$4}' $< > $@

# perform coordinate conversion using liftover chain file
# mm9ToMm10.over.chain.gz is from UCSC, uses 0-based coordinate system
Caroli.snp.mm10.bed: Caroli.snp.bed mm9ToMm10.over.chain.gz
	liftOver $^ $@ Caroli.snp.mm10.unmapped.bed

# remove indels (marked by "*")
Caroli.snp.mm10.no_indels.bed: Caroli.snp.mm10.bed
	grep -v "\*" $< > $@

# removing "chr" again 
Caroli.snp.mm10.no_indels.no_chr.bed: Caroli.snp.mm10.no_indels.bed
	sed 's/^chr\|%$$//g' $< > $@

# sorting by chromosome (first column) and start&end (second&third column)
# input file for downstream steps
Caroli.snp.mm10.no_indels.no_chr.bed.sorted.gz: Caroli.snp.mm10.no_indels.no_chr.bed
	sort -V -k1,1 -k2,2 -k3,3  $< | gzip -c > $@
	
	
	
# everything in FASTA_MM10 covered by Caroli.snp.mm10.bed is changed to N
# this file isn't required 
#$(FASTA_CAROLI): Caroli.snp.mm10.bed
#module load bedtools;
#bedtools maskfasta -fi $(FASTA_MM10) -fo $@ -bed $<

# same as above, but without indels
# this file isn't required 
#$(FASTA_CAROLI_NOINDELS): Caroli.snp.mm10.no_indels.bed
#module load bedtools;
#bedtools maskfasta -fi $(FASTA_MM10) -fo $@ -bed $<

