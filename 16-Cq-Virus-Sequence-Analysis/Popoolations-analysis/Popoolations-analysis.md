# Popoolation Analysis to Compare Virus from Dv-1 Cells and DinnDiNV Cells 

Downloaded [Popoolation2](https://sourceforge.net/projects/popoolation2/) and put in my programs folder. "PoPoolation2 allows to compare allele frequencies for SNPs between two or more populations and to identify significant differences. PoPoolation2 requires next generation sequencing data of pooled genomic DNA (Pool-Seq)" Because each virus sample is likely a pool of viruses, this seems to be the best analysis. 

Made a new folder for these analyses and copy in the two BAM files that were mapped to the virus. 

Basing the analysis off of their [tutorial](https://sourceforge.net/p/popoolation2/wiki/Tutorial/). Because I already have BAM files for the samples, I will start with making the synconized files section:

"Synchronized files are the main input files for PoPoolation2. They basically contain the allele frequencies for every population at every base in the reference genome in a concise format. Note that the synchronized file format contains the allele frequencies after filtering for base quality."

Using first the samtools mpileup function to make an mpileup file. The -B flag says "disable BAQ (per-Base Alignment Quality)" which I'm not sure why this flag is needed, but it's in the tutorial. 

`samtools mpileup -B 16Cq-trim-h-mapped-DiNV.bam KM_3q-trim-h-mapped-DiNV.bam > 16_KM.mpileup`

Then create the synchronized file using a custom perl script from the program. 

`perl /Users/maggieschedl/Documents/programs/popoolation2_1201/mpileup2sync.pl --min-qual 20 --input 16_KM.mpileup --output 16_KM.sync`

The synchronized file has these columns:
- col1: reference contig
- col2: position within the refernce contig
- col3: reference character
- col4: allele frequencies of population number 1
- col5: allele frequencies of population number 2
- coln: allele frequencies of population number n

"The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A', count of bases 'T',... and deletion count in the end (character '*' in the mpileup)"

My sync file looks similar to the tutorial file, so this seems good. 
Might want to do some summary statistics with the synchronized file. 

Now to calculate allele frequency differences between the two samples. I think that min-count is the minimum number of alleles, because the coverage for these samples is kind of low for pool seq data, I am setting it to 2. Same with the min-coverage, 
my coverage is not very high on some places in the genome. I think the max coverage is fine. 

`perl /Users/maggieschedl/Documents/programs/popoolation2_1201/snp-frequency-diff.pl --input 16_KM.sync --output-prefix 16_KM --min-count 2 --min-coverage 5 --max-coverage 200`

This gives me a file with allele frequency differences for 419 sites. This is ~0.3% of the DiNV genome. 