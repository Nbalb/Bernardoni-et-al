# Quantify fastq using salmon in a for loop
# index_orig_cdna was created using (http://ftp.ensembl.org/pub/release-103/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz)
# which is the previous release, needed in order to use tximeta later. Always check at https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html#Pre-computed_checksums
# which is the latest tximeta-supported release of transcriptome and gtf.
cd $linux/mycmax/salmon_alignments
transcripts="Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"

salmon index -t $transcripts -i index_orig_cdna -k 31

# GTF (http://ftp.ensembl.org/pub/release-103/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.103.gtf.gz)
# Make sure to change the -o option to match the output folder
index="$linux/mycmax/salmon_alignments/index_orig_cdna"
cd $linux/mycmax/fastq_merged

for rname in *_R1.fastq.gz
do

base=`basename $rname _R1.fastq.gz`
echo "Doing $base"
time salmon quant -i $index \
-l A \
-p 5 \
-1 ${base}_R1.fastq.gz \
-2 ${base}_R2.fastq.gz \
--validateMappings \
--gcBias \
-o $linux/mycmax/salmon_alignments/${base}_quant_orig_index

done

