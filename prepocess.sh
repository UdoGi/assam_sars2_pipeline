AdapterRemoval --file1 sampleName_R1_001.fastq.gz --file2 sampleName_R2_001.fastq.gz --threads 32 --basename test --gzip --trim5p 31 --trim3p 31 --minlength 0 --minquality 30 --output1 sampleName_R1_trimmed.fastq.gz --output2 sampleName_R2_trimmed.fastq.gz

fastp --in1 sampleName_R1_trimmed.fastq.gz --in2 sampleName_R2_trimmed.fastq.gz
    --out1 sampleName_R1.fastp.fastq.gz --out2 sampleName_R2.fastp.fastq.gz
    --unpaired1 sampleName_single.fastp.fastq.gz --unpaired2 sampleName_single.fastp.fastq.gz
    --merge --merged_out sampleName_merged.fastp.fastq.gz --dedup --length_required 10 -p
    --json sampleName.json --html sampleName.html --thread 4

touch {output.single} {output.merge};
cat {output.single} {output.merge} > {output.singletons}

cat sampleName_single.fastp.fastq.gz sampleName_merged.fastp.fastq.gz > sampleName_singletons.fastq.gz
