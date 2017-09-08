#STRAINS=Case1a Case1b Case2a Case2b Case3a Case3b Case3c Case3d Case3e Case3f Case4a Case5a Case5b Case5c Case6a Case6b Case6c Case7a Case7b Case8a Case8b Case8c
STRAINS=Case1a Case1b Case2a Case2b Case3a Case3b Case3c Case4a Case5a Case5b Case5c Case6a Case6b Case6c Case7a Case7b Case8a Case8b Case8c

CPUS=1
REF=NC_000962

## file locations
GENOMES    = genomes
ALIGNMENTS = alignments
PHYLOGENY  = phylogeny
DATA       = data
VCF        = vcf

## exe locations
EXE=./bin
TMAP=$(EXE)
SAMTOOLS=$(EXE)

## Workflow

indexed = $(addsuffix .tmap.bwt, $(GENOMES)/$(REF).fna )
bam     = $(addprefix $(ALIGNMENTS)/$(REF)_, $(addsuffix .bam, $(STRAINS) ) )
vcf     = $(addprefix $(VCF)/$(REF)_, $(addsuffix .vcf, $(STRAINS) ) )
sites   = $(addprefix $(VCF)/$(REF)_, $(addsuffix .all.vcf.gz, $(STRAINS) ) )

all: index alignments sites trees
index: $(indexed)
alignments: $(bam)
sites: $(sites)
snps: $(vcf)

versions:
	@samtools 2>&1 | grep Version | perl -p -e 's/Version/samtools/' > versions.txt
	@bwa 2>&1 | grep Version | perl -p -e 's/Version/bwa/' >> versions.txt
	@raxmlHPC-PTHREADS-SSE3 -v | grep "RAxML version" | perl -p -e 's/.+RAxML version (\d+.\d+.\d+).+/RAxML: \1/' >> versions.txt

$(GENOMES)/$(REF).fna.tmap.bwt:
	$(TMAP)/tmap index -f $(GENOMES)/$(REF).fna

$(VCF)/$(REF)_%.all.vcf.gz: $(ALIGNMENTS)/$(REF)_%.bam
	$(SAMTOOLS)/samtools mpileup -gf $(GENOMES)/$(REF).fna $(ALIGNMENTS)/$(REF)_$*.bam | $(SAMTOOLS)/bcftools view -cg - | $(EXE)/bgzip > $(VCF)/$(REF)_$*.all.vcf.gz
	$(EXE)/tabix -p vcf $(VCF)/$(REF)_$*.all.vcf.gz

#	$(EXE)/bgzip $(VCF)/$(REF)_$*.all.vcf

$(ALIGNMENTS)/$(REF)_%.bam: $(GENOMES)/$(REF).fna.tmap.bwt
	$(TMAP)/tmap mapall -f $(GENOMES)/$(REF).fna -r $(DATA)/$*.bam -n $(CPUS) -v -Y -u -o 1 stage1 map4 | $(SAMTOOLS)/samtools sort -o - $(ALIGNMENTS)/temp.$* | $(SAMTOOLS)/samtools rmdup - - > $@
	$(SAMTOOLS)/samtools index $@

$(PHYLOGENY)/NC_000962.b1.infile:
	perl scripts/snp_caller.pl --chrom NC_000962.3 --qual 30 --dp 2 --dp4 75 --dpmax 5000 --af 0 --mq 30 --noindels --noheader \
		--vcf Case5a,vcf/NC_000962_Case5a.all.vcf.gz \
		--vcf Case5b,vcf/NC_000962_Case5b.all.vcf.gz \
		--vcf Case5c,vcf/NC_000962_Case5c.all.vcf.gz \
		--vcf Case4a,vcf/NC_000962_Case4a.all.vcf.gz \
		--vcf Case2a,vcf/NC_000962_Case2a.all.vcf.gz \
		--vcf Case1a,vcf/NC_000962_Case1a.all.vcf.gz \
		--vcf Case1b,vcf/NC_000962_Case1b.all.vcf.gz \
		--vcf Case3a,vcf/NC_000962_Case3a.all.vcf.gz \
		--vcf Case3b,vcf/NC_000962_Case3b.all.vcf.gz \
		--vcf Case3c,vcf/NC_000962_Case3c.all.vcf.gz \
		--vcf Case3d,vcf/NC_000962_Case3d.all.vcf.gz \
		--vcf Case2b,vcf/NC_000962_Case2b.all.vcf.gz \
		--vcf Case6c,vcf/NC_000962_Case6c.all.vcf.gz \
		--vcf Case6b,vcf/NC_000962_Case6b.all.vcf.gz \
		--vcf Case6a,vcf/NC_000962_Case6a.all.vcf.gz \
 		\
		--vcf Case7b,vcf/NC_000962_Case7b.all.vcf.gz \
		--vcf Case7a,vcf/NC_000962_Case7a.all.vcf.gz \
		\
 		--vcf Case8a,vcf/NC_000962_Case8a.all.vcf.gz \
 		--vcf Case8b,vcf/NC_000962_Case8b.all.vcf.gz \
		--vcf Case8c,vcf/NC_000962_Case8c.all.vcf.gz \
	--dir $(PHYLOGENY)/ --phylip NC_000962.b1.infile --verbose 1 --refilter -b 1 --cpus $(CPUS)

$(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1: $(PHYLOGENY)/NC_000962.b1.infile
	$(EXE)/raxmlHPC-PTHREADS-SSE3 -T $(CPUS) -f a -s $(PHYLOGENY)/NC_000962.b1.infile -x 12345 -p 1234 -# 1000 -m GTRGAMMA -w $(shell pwd)/$(PHYLOGENY) -n NC_000962.b1 -o NC_000962.3
	perl -pi -e "s/NC_000962.3/H37Rv/" $(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1

trees: $(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1

tree.png:
	Rscript scripts/tree.R

