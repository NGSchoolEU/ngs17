

ll data/
tmap index -f genomes/NC_000962.fna
tmap-mapall genomes/NC_000962.fna data/Case1a.bam > Case1a.bam
samtools sort Case1a.bam Case1a.sort
samtools index Case1a.sort.bam
samtools flagstat Case1a.sort.bam
samtools view Case1a.sort.bam | head
samtools-mpileup genomes/NC_000962.fna Case1a.sort.bam > Case1a.all.vcf.gz
tabix -p vcf Case1a.all.vcf.gz

tabix vcf/Case1a.all.vcf.gz NC_000962.3:761094-761096
tabix vcf/Case1a.all.vcf.gz NC_000962.3:761154-761156
tabix vcf/Case1a.all.vcf.gz NC_000962.3:2155167-2155169
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247429-4247431
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247729-4247731
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7569-7571
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7581-7583

snp-filter
estimate-tree

