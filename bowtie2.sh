MN908947.3

wget -O sars2_ref.fasta https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?db\=nucleotide\&id\=MN908947.3\&rettype\=fasta\&retmode\=text

mamba env create -f alignment_consensus.yaml --name alignment_consensus_india

bowtie2-build sars2_ref.fasta sars2_ref

bowtie2 -x sars2_ref -U sampleName_singletons.fastq.gz -p 4 --sensitive
        --no-unal 2> stderr.log | samtools view -bS - > sampleName_single.bam


bowtie2 -x sars2_ref -1 sampleName_R1.fastp.fastq.gz -2 sampleName_R2.fastp.fastq.gz -p 4 --sensitive --no-unal | samtools view -bS - > sampleName_paired.bam

samtools merge sampleName_merged.bam sampleName_single.bam sampleName_paired.bam
samtools sort sampleName_merged.bam -o sampleName_sorted.bam
samtools index sampleName_sorted.bam


ivar trim -i sampleName_sorted.bam -b SARS-CoV-2.primer.bed -p sampleName_sorted_primer_removed

samtools sort sampleName_sorted_primer_removed.bam -o sampleName_sorted_primer_removed_sorted.bam

samtools mpileup -aa -A -d 0 -Q 0 sampleName_sorted_primer_removed_sorted.bam | ivar consensus -p sampleName_consensus


