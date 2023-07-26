## Calling SNPs with freebayes

Install with `conda install -c bioconda freebayes`

Want to do this analysis for a "pooled" sample because there is a pool of virus genomes in the samples I'm looking at.  
Because we don't know exactly how to analyze this data, I want it to give me the allele frequencies if I assume the number of alleles (average coverage) or if we don't know the number of alleles. 

Here are the two main flags I want to use: 

 -K --pooled-continuous:  
    Output all alleles which pass input filters, regardles of
    genotyping outcome or model.

-J --pooled-discrete:  
    Assume that samples result from pooled sequencing.
    Model pooled samples using discrete genotypes across pools.
    When using this flag, set --ploidy to the number of
    alleles in each sample or use the --cnv-map to define
    per-sample ploidy.

I will first start with pooled-continuous, which is to "Generate frequency-based calls for all variants passing input thresholds. You'd do this in the case that you didn't know the number of samples in the pool." Also it has been described as "FreeBayes can act as a frequency-based pooled caller and describe variants and haplotypes in terms of observation frequency rather than called genotypes. To do so, use --pooled-continuous and set input filters to a suitable level. Allele observation counts will be described by AO and RO fields in the VCF output."

For running freebayes, I need to add readgroups to my bam files. I could have done this when making the files, but I didn't know that at the time. Read groups basically tag each file and read with the sample information, that way when the program is going through the reads to find SNPs it knows which sample they are from. 

This can be done with a samtools command: samtools addreplacerg – adds or replaces read group tags. I'll use the -r flag which "Allows you to specify a read group line to append to the header and applies it to the reads specified by the -m option. If repeated it automatically adds in tabs between invocations." The default m option is to overwrite all existing read groups, which is fine because there are none. ID and SM (sample) are basically the same here because I only have two samples to separate. PL is for the sequencing type. If I had more than two samples, the SM is what would separate them out in the vcf file. 

`samtools addreplacerg -r ID:Dv-1 -r SM:VIR -r PL:ILLUMINA -o KM3-trim-h-mapped-re-mapped-DiNVrg.bam KM3-trim-h-mapped-re-mapped-DiNV.bam`  
`samtools addreplacerg -r ID:D-inn -r SM:INN -r PL:ILLUMINA -o 16cq-combo-mapped-re-mapped-DiNvrg.bam 16cq-combo-mapped-re-mapped-DiNv.bam`


After the read groups are added, the bam files are indexed 

`samtools index 16cq-combo-mapped-re-mapped-DiNvrg.bam`   
`samtools index KM3-trim-h-mapped-re-mapped-DiNVrg.bam`  

Then for --pooled-continuous there are some flags I want to use

 -F --min-alternate-fraction N :   
    Require at least this fraction of observations supporting
    an alternate allele within a single individual in the
    in order to evaluate the position.  default: 0.2

 -C --min-alternate-count N :  
    Require at least this count of observations supporting
    an alternate allele within a single individual in order
    to evaluate the position.  default: 2

For -C I want it to be 1, because with the pooled samples, even 1 read could represent a unique virus "individual".  
For -F I am not sure what the best number is here. If I use 1/20th (0.05), this might capture all the variants in the 16Cq sample, but miss some in the Dv-1 sample. I could use 1/50 (0.02) as well. Right now I think I will make it even lower, 0.01, because I can always filter out SNPs later that don't give any information or might just be sequencing error. 

I need to provide the path to the genome, and I can just list the two bam files I want it to use 

`freebayes -f /Users/maggieschedl/Desktop/KU/sequences/16Cq-DiNV-Test/mapping-combo/GCF_004132165.1_DiNV_CH01M_genomic.fna -F 0.01 -C 1 --pooled-continuous 16cq-combo-mapped-re-mapped-DiNvrg.bam KM3-trim-h-mapped-re-mapped-DiNVrg.bam > continuous-DiNV.vcf`



