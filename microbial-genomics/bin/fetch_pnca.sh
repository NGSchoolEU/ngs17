for i in alignments/*.bam ; do 
j=`echo $i | sed 's/alignments\/NC_000962_//' | sed 's/.bam//'`;
./bin/samtools mpileup -uf genomes/NC_000962.fna \
  -r NC_000962.3:2288681-2289241 $i \
  | ./bin/bcftools view -cg - \
  | ./bin/vcfutils.pl vcf2fq - \
  | seqtk trimfq - \
  | seqtk seq -rA - \
  | transeq -filter \
  | perl -p -e "s/NC_000962.3_1/${j}/";
done

