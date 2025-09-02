import os
import pandas as pd
import Bio


configfile: "config.yaml"

innerbc_sheet = pd.read_csv(config["innerbc_sheet"])  # Generate inner_bc names by reading the innerbc_96.csv
INNER_BC = innerbc_sheet['barcodeRV_id'].tolist()  # Extract the list of inner bc
BATCHID = config["batch_id"]
OUTDIR = config["output_dir"]
LOGDIR = config["log_dir"]
BARCODE = config["input_fastq"].keys()


#rule all:
#    input:
#        OUTDIR + "/stats/"  + "readcounts_after_trimming.csv"

rule all:
    input:
        OUTDIR + "/stats/" + "readcount_events_gene.csv",
        OUTDIR + "/stats/" + "readcount_cas12a_cassette.csv",
        OUTDIR + "/stats/" + "readcount_offtarget.csv",
        OUTDIR + "/stats/" + "readcount_misAssignedM1M1.csv",
        OUTDIR + "/stats/" + "final_outcomes_meta.csv"


# pre-mapping --- demultiplexing and trimming
rule innerbc_demux_trim:
    input:
        innerbc_to_trim=config["innerbc_sheet"],    
        fastq=lambda wildcards: config["input_fastq"][wildcards.barcode],
        env="envs/trim_env"
    output:
        directory(OUTDIR + "/seq/{barcode}")   
    log:
        LOGDIR + "/{barcode}.log"
    threads: 1
    shell:
        "echo {input.fastq} && echo {input.innerbc_to_trim} && "
        "julia --project={input.env} julia/nano_read_demux_dl2.jl {input.fastq} {input.innerbc_to_trim} "
        "{output} 2>&1 | tee {log}"

rule dir_exist:
    input:
        OUTDIR + "/seq/{barcode}"
    output:
        OUTDIR + "/seq/{barcode}_bc_{inner_barcode}.fastq.gz"
    threads: 1
    shell:
        "mv "+ OUTDIR + "/seq/{wildcards.barcode}/{wildcards.inner_barcode}.fastq.gz {output}" 


#1 index reference genome 
rule indexgenome:
    input:
        genome = "{genome}"
    output:
        bwt =  "{genome}.bwt"
    threads: 1
    shell:
        "bwa index {input.genome}"

#3 mapping to reference sequences
rule bwa_mem_map:
    input:
        genome = config["ref_genome"],
        read = OUTDIR + "/seq/" + "{barcode}_bc_{innerbc}.fastq.gz",
        index = config["ref_genome"] + ".bwt"
        #trimmed_read = get_reads
    output:
        temp(OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.raw.bam")
    log:
        LOGDIR + "/{barcode}_{innerbc}.log"
    threads: 8
    shell:
        "(bwa mem -x ont2d -t {threads} {input.genome} {input.read}| "
        "samtools view -bS - > {output}) 2> {log}"


#3.1 sort bam file
rule sortbam:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.raw.bam"
    output:
        protected(OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.bam")
    threads: 1
    shell:
        "samtools sort {input} -o {output}"

rule filterMisassignedM1M1:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.bam"
    output:
        fil1_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.bam",
        rm_read = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.misassigedM1M1.tsv"
    threads: 1
    shell:
        "samtools view -h {input} > {input}.tmp.sam && "
        "python python/003_filter_bam_SA_v2.py -in {input}.tmp.sam -out {input}.tmpout.sam -read {output.rm_read} -samid \"{wildcards.barcode}\t{wildcards.innerbc}\" && "
        "samtools view -bS {input}.tmpout.sam > {output.fil1_bam} && "
        "rm {input}.tmp.sam {input}.tmpout.sam"

rule combine_SA_M1M1:
    input:
        expand(OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.misassigedM1M1.tsv", barcode = config["sample"].keys(), innerbc = INNER_BC)
    output:
        OUTDIR + "/stats/" + "misAssignedM1M1_summary.tsv"
    threads: 1
    shell:
        "echo \"barcode\tinner_bc\tread_id\tgene\tdescription\" > {output} && "
        "cat {input} >> {output}"

#3.2 index bam file
rule indexbam:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    threads: 1
    shell:
        "samtools index {input}"

#4.1 filter MAPQ score >=1 #alternative: >=10
rule filterMAPQ:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.bam"
    output:
        mapqhigh=OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.bam",
        mapqlow=OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.mapq0.bam"
    threads: 1
    shell:
        "samtools view -bSq 1 -U {output.mapqlow} {input} > {output.mapqhigh} "

rule filterHDRreads:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.bam"
    output:
        fil3=OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.bam",
        indel=OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.indel.bam"
    log:
        LOGDIR + "/{barcode}_{innerbc}.filterHDRreads.log"
    threads: 1
    shell:
        "python python/004_filter_HDR_reads_v2.py -in {input} -hdr {output.fil3} -indel {output.indel} > {log} 2>&1"


rule filterPrimarilyAlignedreads:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.bam"
    output:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.pri.bam"
    threads: 1
    shell:
        "samtools view -b -F 0x900 {input} > {output} "


rule filterXAtoFASTQ:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.mapq0.bam"
    output:
        OUTDIR + "/seq/" + "{barcode}_bc_{innerbc}.unmapped.fastq"
    threads: 1
    shell:
        "samtools view -h {input} > {input}.tmp.sam && "
        "python python/006_soft_clipping_to_fastq.py -in {input}.tmp.sam -out {output} && "
        "rm {input}.tmp.sam "

rule bwa_mem_map_offtar:
    input:
        genome = config["offtar_ref_genome"],
        read = OUTDIR + "/seq/" + "{barcode}_bc_{innerbc}.unmapped.fastq",
        index = config["offtar_ref_genome"] + ".bwt"
    output:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.unmapped.bam"
    log:
        LOGDIR + "/{barcode}_{innerbc}_unmapped.log"
    threads: 8
    shell:
        "(bwa mem -x ont2d -t {threads} {input.genome} {input.read}| "
        "samtools view -bS - > {output}) 2> {log}"

rule get_offtar_tsv:
    input:
        OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.unmapped.bam"
    output:
        OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.offtar.tsv"
    threads: 1
    shell:
        "samtools view -h {input} > {input}.tmp.sam && "
        "python python/007_parse_bam_to_table.py -in {input}.tmp.sam -out {output} -samid \"{wildcards.barcode}\t{wildcards.innerbc}\" && "
        "rm {input}.tmp.sam "

rule combine_tsv:
    input:
        expand(OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.offtar.tsv", barcode = config["sample"].keys(), innerbc = INNER_BC)
    output:
        OUTDIR + "/stats/" + "off_target_summary.tsv"
    threads: 1
    shell:
        "echo \"barcode\tinner_bc\tread_id\tgene\tchr\tstrand\tcoord\tlen_mapped\tlen_read\" > {output} && "
        "cat {input} >> {output}"
    

#5 generate statistcs files
rule generateMapStats:
    input:
        nofil_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.bam",
        nofil_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.bam.bai",
        fil1_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.bam",
        fil1_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.bam.bai",
        fil2_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.bam",
        fil2_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.bam.bai",
        indel_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.indel.bam",
        indel_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.indel.bam.bai",
        fil3_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.bam",
        fil3_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.bam.bai",
        pri_bam = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.pri.bam",
        pri_bai = OUTDIR + "/bwa/" + "{barcode}/" + "{innerbc}.sorted.fil1.fil2.fil3.pri.bam.bai"

    output:
        nofil_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.flagstat",
        nofil_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.idxstats",        
        fil1_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.flagstat",
        fil1_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.idxstats",
        fil2_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.flagstat",
        fil2_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.idxstats",
        indel_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.indel.flagstat",
        indel_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.indel.idxstats",
        fil3_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.fil3.flagstat",
        fil3_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.fil3.idxstats",
        pri_flagstat = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.fil3.pri.flagstat",
        pri_idxstats = OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.fil3.pri.idxstats"
        
    shell:
        "samtools flagstat {input.nofil_bam} > {output.nofil_flagstat} && "
        "samtools idxstats {input.nofil_bam} > {output.nofil_idxstats} && "
        "samtools flagstat {input.fil1_bam} > {output.fil1_flagstat} && "
        "samtools idxstats {input.fil1_bam} > {output.fil1_idxstats} && "
        "samtools flagstat {input.fil2_bam} > {output.fil2_flagstat} && "
        "samtools idxstats {input.fil2_bam} > {output.fil2_idxstats} && "
        "samtools flagstat {input.indel_bam} > {output.indel_flagstat} && "
        "samtools idxstats {input.indel_bam} > {output.indel_idxstats} && "
        "samtools flagstat {input.fil3_bam} > {output.fil3_flagstat} && "
        "samtools idxstats {input.fil3_bam} > {output.fil3_idxstats} && "
        "samtools flagstat {input.pri_bam} > {output.pri_flagstat} && "
        "samtools idxstats {input.pri_bam} > {output.pri_idxstats} "


#7 summarize read counts from fil1.fil2.fil3.pri/fil1.fil2.indel.idxstats with R
rule idxstats_summary:
    input:
        pri_idxstats = expand(OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.fil3.pri.idxstats", barcode = config["sample"].keys(), innerbc = INNER_BC),
        indel_idxstats = expand(OUTDIR + "/stats/" + "{barcode}/" + "{innerbc}.fil1.fil2.indel.idxstats", barcode = config["sample"].keys(), innerbc = INNER_BC)
    output:
        OUTDIR + "/stats/" + "fil3_pri_idxstats_summary.csv",
        OUTDIR + "/stats/" + "fil2_indel_idxstats_summary.csv",
    script:
        "R/001_idxstats_summary_v2.R"


#8 transform "idxstats_summary.csv" and "off_target_summary.tsv" tables into wide format
rule counts_events_summary:
    input:
        OUTDIR + "/stats/" + "fil3_pri_idxstats_summary.csv",
        OUTDIR + "/stats/" + "fil2_indel_idxstats_summary.csv",
        OUTDIR + "/stats/" + "off_target_summary.tsv",
        OUTDIR + "/stats/" + "misAssignedM1M1_summary.tsv",
        METADATA = config["metadata_sheet"]
    output:
        OUTDIR + "/stats/" + "readcount_events_gene.csv",
        OUTDIR + "/stats/" + "readcount_cas12a_cassette.csv",
        OUTDIR + "/stats/" + "readcount_offtarget.csv",
        OUTDIR + "/stats/" + "readcount_misAssignedM1M1.csv",
        OUTDIR + "/stats/" + "final_outcomes_meta.csv"
    script:
        "R/002_counts_events_summary_v2.R"



