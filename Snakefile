from pathlib import Path
raw_data_dir=Path(config['path_to_samples']).resolve().absolute()
# print(raw_data_dir)
fq_paths=list(raw_data_dir.glob('*.fastq.gz'))
# print(fq_paths)
all_samples = [p.stem.split('_R1')[0] for p in fq_paths if "_R1" in p.stem]
all_samples = [e for e in all_samples if 'Undetermined' not in e]
# print(fq_paths)
print(len(all_samples))
print(all_samples[:2])

rule all:
    input:
        # fastp=[f'fastp/{sample}_singletons.fastp.fastq.gz' for sample in all_samples[:]],
        # fastp=[f'fastp/{sample}.fastp.json' for sample in all_samples[:]],
        #  "multiqc_report.html"
        # db="bowtie2/bowtie2_db/sars2_ref.1.bt2",
        # bam=['bowtie2/data/'  + f"{sample}_no_clip.bam" for sample in all_samples[:]]
        # [ f"{sample}_consensus.fasta" for sample in all_samples]
        [ f"{sample}.consensus.done" for sample in all_samples]

rule adapter_removal:
    input:
        R1=str(raw_data_dir) + "/{sample}_R1_001.fastq.gz",
        R2=str(raw_data_dir) + "/{sample}_R2_001.fastq.gz",
    output:
        R1='adapter_removal/' + '{sample}_R1_trimmed.fastq.gz',
        R2='adapter_removal/' + '{sample}_R2_trimmed.fastq.gz',
    conda:
        "preprocess.yaml"
    log:
        stdout="logs/adapter_removal/{sample}.stdout.log",
        stderr="logs/adapter_removal/{sample}.stderr.log"
    threads:
        2
    params:
        basename=lambda wildcards: wildcards.sample
    shell:
        "AdapterRemoval --file1 {input.R1} --file2 {input.R2} --threads {threads} --basename {params.basename} "
        "--gzip --trim5p 31 --trim3p 31 --minlength 30 --minquality 0 --output1 {output.R1} --output2 {output.R2} "
        "--adapter1 CTGTCTCTTATACACATCT --adapter2 CTGTCTCTTATACACATCT > {log.stdout} 2> {log.stderr};"


rule fastp:
    input:
        R1=str(raw_data_dir) + "/{sample}_R1_001.fastq.gz",
        R2=str(raw_data_dir) + "/{sample}_R2_001.fastq.gz",
        # R1='adapter_removal/' + '{sample}_R1_trimmed.fastq.gz',
        # R2='adapter_removal/' + '{sample}_R2_trimmed.fastq.gz',
    output:
        R1="fastp" + "/{sample}_R1.fastp.fastq.gz",
        R2='fastp' + "/{sample}_R2.fastp.fastq.gz",
        unpaired='fastp' + '/{sample}_unpaired.fastp.fastq.gz',
        merge='fastp' + '/{sample}_merged.fastp.fastq.gz',
        singletons='fastp' + '/{sample}_singletons.fastp.fastq.gz',
        json='fastp' + '/{sample}.fastp.json',
        html='fastp' + '/{sample}.html'
    conda:
        "preprocess.yaml"
    log:
        stdout="logs/fastp/{sample}.stdout.log",
        stderr="logs/fastp/{sample}.stderr.log"
    threads:
        2
    params:
        basename=lambda wildcards: wildcards.sample
    shell:
         "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
         "--unpaired1 {output.unpaired} --unpaired2 {output.unpaired} --merge "
         "--merged_out {output.merge} --dedup --length_required 10 -p "
         "--json {output.json} --html {output.html} --thread {threads} > {log.stdout} 2> {log.stderr};"
         "touch {output.unpaired} {output.merge};"
         "cat {output.unpaired} {output.merge} > {output.singletons}"


rule multiqc:
  input:
        expand("fastp/{sample}.fastp.json", sample=all_samples)
  output:
        "multiqc_report.html"
  log:
        "logs/multiqc/logfile.log",
  resources:
        memory=64000,
        time="00:59"
  conda:
        "preprocess.yaml"
  shell:
        "multiqc fastp"

rule download_ref:
    output:
        fasta='sars2_ref.fasta'
    params:
        acc='MN908947.3'
    shell:
        "wget -O {output.fasta} https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?db\=nucleotide\&id\={params.acc}\&rettype\=fasta\&retmode\=text"



rule make_bowtie2_db:
    input:
        '/Users/udo/Documents/code/pythonProject/sars2_ref.fasta'
    output:
        "bowtie2/bowtie2_db/sars2_ref.1.bt2"
    conda:
        # "/Users/udo/opt/miniconda3/envs/darkmatter/"
        'alignment_consensus.yaml'
    log:
        stdout="logs/make_bowtie2_db/stdout.log",
        stderr="logs/make_bowtie2_db/stderr.log",
    threads:
        2
    shell:
        "rm -rf bowtie2/bowtie2_db/ > {log.stdout} 2> {log.stderr} ; mkdir -p bowtie2/bowtie2_db/ >> {log.stdout} 2>> {log.stderr};"
        "cd bowtie2/bowtie2_db/; "
        "bowtie2-build {input} sars2_ref --threads {threads} >> ../../{log.stdout} 2>> ../../{log.stderr}; "


rule run_bowtie2_no_clip:
    input:
        db1="bowtie2/bowtie2_db/sars2_ref.1.bt2",
        R1="fastp" + "/{sample}_R1.fastp.fastq.gz",
        R2='fastp' + "/{sample}_R2.fastp.fastq.gz",
        singletons='fastp' + '/{sample}_singletons.fastp.fastq.gz',
    output:
        bam1=temp('bowtie2/data/' + "{sample}_1_no_clip.bam"),
        bam2=temp('bowtie2/data/' + "{sample}_2_no_clip.bam"),
        bam_merged=temp('bowtie2/data/' + "{sample}_merged_no_clip.bam"),
        bam='bowtie2/data/'  + "{sample}_no_clip.bam",
        bai='bowtie2/data/' + "{sample}_no_clip.bam.bai",
    conda:
        'alignment_consensus.yaml'
    log:
        stdout="logs/bowtie2/run_bowtie2_no_clip{sample}.stdout.log",
        stderr="logs/bowtie2/run_bowtie2_no_clip{sample}.stderr.log"
    threads:
        4
    shell:
        "export BOWTIE2_INDEXES=bowtie2/bowtie2_db/; "
        "bowtie2 -x sars2_ref -U {input.singletons} -p {threads} --sensitive --no-unal 2> {log.stderr} "
        "| samtools view -bS - > {output.bam1} 2>> {log.stderr}; "
        "bowtie2 -x sars2_ref -1 {input.R1} -2 {input.R2} -p {threads} --sensitive --no-unal 2> {log.stderr} "
        "| samtools view -bS - > {output.bam2} 2>> {log.stderr}; "
        "samtools merge {output.bam_merged} {output.bam1} {output.bam2} >> {log.stdout} 2>> {log.stderr}; "
        "samtools sort {output.bam_merged} -o {output.bam} --threads {threads} >> {log.stdout} 2>> {log.stderr}; "
        "samtools index {output.bam} >> {log.stdout} 2>> {log.stderr}; "


rule make_consensus:
    input:
        bam='bowtie2/data/'  + "{sample}_no_clip.bam",
        bai='bowtie2/data/' + "{sample}_no_clip.bam.bai",
    output:
        # fasta="{sample}_consensus.fasta"
        touch("{sample}.consensus.done")
    conda:
        'alignment_consensus.yaml'
    log:
        stdout="logs/make_consensus/{sample}.stdout.log",
        stderr="logs/make_consensus/{sample}.stderr.log"
    params:
        prefix= lambda wildcards: f"{wildcards.sample}_consensus"
    threads:
        4
    shell:
        "samtools mpileup -aa -A -d 0 -Q 0 {input.bam} | ivar consensus -p {params.prefix} > {log.stdout} 2> {log.stderr}; "


# rule pangolineage:
#     input:
#         fasta='consensus/' + "{sample}_consensus.fasta"
#     output:
#         csv='pangolineage' + "{sample}_pangolineage.csv"
#     # conda:
#     #     'alignment_consensus.yaml'
#     log:
#         stdout="logs/pangolineage/{sample}.stdout.log",
#         stderr="logs/pangolineage/{sample}.stderr.log"
    # params:
    #     prefix= lambda wildcards: f"{wildcards.sample}_consensus"
    # threads:
    #     4
    # shell:
    #     "samtools mpileup -aa -A -d 0 -Q 0 {input.bam} | ivar consensus -p {params.prefix} > {log.stdout} 2> {log.stderr} "

