# Mapping Trimmed Reads to the Various Reference Genomes 


**Downloading Genomes**

`GCF_004132165.1_DiNV_CH01M_genomic.fna` is the DiNV genome, downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2057187)

`GCA_007989325.2_ASM798932v2_genomic.fna` is the D. virilis genome, used for Dv-1 cells, downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7244), ASM798932v2 chromosome assembled.

`GCA_004354385.2_ASM435438v2_genomic.fna` is the D. innubila geneome, used for the DinnDiNV cells, downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=198719), ASM435438v2 chromosome assembled from 2020

Moved all the trimmed filed into a new folder for mapping 

**Map 16 Cq samples to DiNV genome**

- [BWA](https://bio-bwa.sourceforge.net/bwa.shtml) is an alignment tool that is widely used, and that I have [used previously](https://github.com/meschedl/DiNV-Dv1-Genome-Integration/blob/main/Mapping/mapping_Dvir_DiNV_combined.md)
- First "BWA first needs to construct the FM-index for the reference genome"

`bwa index GCF_004132165.1_DiNV_CH01M_genomic.fna`

- Output:
```
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.03 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.01 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index GCF_004132165.1_DiNV_CH01M_genomic.fna
[main] Real time: 0.045 sec; CPU: 0.048 sec
```
- Now to map the 16 Cq samples to the genome. BWA outputs files in SAM format, which are really large. I can use [samtools](http://www.htslib.org/doc/samtools.html) to convert the SAM file to a BAM file
    - The command samtools view prints all alignment files, and the options -h means to include the header in the output
    - The file also needs to be sorted so that it is in order (it will be random otherwise) and this is done with samtools sort which will sort the alignments by the leftmost coordinates
    - These commands can all be "piped" together because they use the output of one command as the input of the next. The "|" is the pipe
    - I can use both the paired reads in this mapping
`bwa mem GCF_004132165.1_DiNV_CH01M_genomic.fna 16Cq-DiNV-R1-trim.fastq.gz 16Cq-DiNV-R2-trim.fastq.gz | samtools view -h | samtools sort - > 16Cq-trim-h-mapped-DiNV.bam`

**Map 16 Cq sample to D. innubila genome**

- Index the innubila genome

`bwa index GCA_004354385.2_ASM435438v2_genomic.fna`

- Oupput: 
```
[bwa_index] 59.74 seconds elapse.
[bwa_index] Update BWT... 0.58 sec
[bwa_index] Pack forward-only FASTA... 0.45 sec
[bwa_index] Construct SA from BWT and Occ... 39.30 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index GCA_004354385.2_ASM435438v2_genomic.fna
[main] Real time: 100.979 sec; CPU: 100.913 sec
```
- Map the 16 Cq samples to the innubila genome, using similar code as above. However this time I will try multi-threading so it is a little quicker 

`bwa mem -t 5 GCA_004354385.2_ASM435438v2_genomic.fna 16Cq-DiNV-R1-trim.fastq.gz 16Cq-DiNV-R2-trim.fastq.gz | samtools view --threads 5 -h | samtools sort --threads 5 - > 16Cq-trim-h-mapped-inn.bam`

**Map Dv-1 sample to DiNV genome**

- Genome is already indexed, so that's good
- Map to the DiNV genome, using multi-threading

`bwa mem -t 5 GCF_004132165.1_DiNV_CH01M_genomic.fna KM_3_1_trim.fq.gz KM_3_2_trim.fq.gz | samtools view --threads 5 -h | samtools sort --threads 5 - > KM_3q-trim-h-mapped-DiNV.bam`

**Map Dv-1 sample to D. virilis genome**

- Need to index this genome first 

`bwa index GCA_007989325.2_ASM798932v2_genomic.fna`

- Output:

```
[bwa_index] Pack FASTA... 0.84 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=339543490, availableWord=35891056
[BWTIncConstructFromPacked] 10 iterations done. 59204242 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 109375826 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 153964450 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 193590914 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 228806962 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 260103058 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 287915090 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 312630434 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 334593490 characters processed.
[bwt_gen] Finished constructing BWT in 93 iterations.
[bwa_index] 59.84 seconds elapse.
[bwa_index] Update BWT... 0.59 sec
[bwa_index] Pack forward-only FASTA... 0.46 sec
[bwa_index] Construct SA from BWT and Occ... 36.62 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index GCA_007989325.2_ASM798932v2_genomic.fna
[main] Real time: 98.405 sec; CPU: 98.342 sec
```

- Map Dv-1 samples to the D. virilis genome, using multithreading 

`bwa mem -t 5 GCA_007989325.2_ASM798932v2_genomic.fna KM_3_1_trim.fq.gz KM_3_2_trim.fq.gz | samtools view --threads 5 -h | samtools sort --threads 5 - > KM_3q-trim-h-mapped-Dvir.bam`


## Determining mapping percentages 

Using [Qualimap] because it "Examines sequencing alignment data according to the features of the mapped reads and their genomic properties"

Download [Qualimap](http://qualimap.conesalab.org/doc_html/intro.html#installation)

Download [Java](https://www.java.com/en/)

Try to run in [terminal](http://qualimap.conesalab.org/doc_html/command_line.html) instead of using the GUI. Going to use BAM QC: "BAM QC reports information for the evaluation of the quality of the provided alignment data." This outputs all of the mapping statistics into a pdf file for each bam.  

Run as a loop to do all of the bam files, adding in -ip for "Activate this option to collect statistics of overlapping paired-end reads" 

```
for bam in *.bam 
do
/Users/maggieschedl/Documents/programs/qualimap_v2.2.1/qualimap bamqc -bam $bam -ip -outfile ${bam}.pdf
done
```



