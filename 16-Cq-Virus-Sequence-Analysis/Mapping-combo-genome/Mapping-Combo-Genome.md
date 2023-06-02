## Mapping Reads to Combined Virus and Fly Genomes 

Here I want to map the reads to combined genomes of either DiNV-D.innubila or DiNV-D.virilis. 

Using a new folder for mapping, first I used cat to combine genomes together 

DiNV-D. virilis  
`cat GCA_007989325.2_ASM798932v2_genomic.fna GCF_004132165.1_DiNV_CH01M_genomic.fna > DiNV-Dvir.fna`

DiNV-D. innubila  
`cat GCA_004354385.2_ASM435438v2_genomic.fna GCF_004132165.1_DiNV_CH01M_genomic.fna > DiNV-Dinn.fna`

Again, Ill be using [BWA](https://bio-bwa.sourceforge.net/bwa.shtml) mem to do the mapping, and I'll need to index the combo genomes before doing the real mapping. 

`bwa index DiNV-Dvir.fna`

Output:
```
[bwa_index] Pack FASTA... 0.85 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=339854600, availableWord=35912944
[BWTIncConstructFromPacked] 10 iterations done. 59240328 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 109442536 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 154058360 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 193708968 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 228946536 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 260261704 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 288090696 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 312821144 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 334797592 characters processed.
[bwt_gen] Finished constructing BWT in 93 iterations.
[bwa_index] 59.50 seconds elapse.
[bwa_index] Update BWT... 0.59 sec
[bwa_index] Pack forward-only FASTA... 0.45 sec
[bwa_index] Construct SA from BWT and Occ... 36.98 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index DiNV-Dvir.fna
[main] Real time: 98.430 sec; CPU: 98.372 sec
```
`bwa index DiNV-Dinn.fna`

Output:
```
[bwa_index] Pack FASTA... 0.82 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=332880780, availableWord=35422588
[BWTIncConstructFromPacked] 10 iterations done. 58431468 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 107948108 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 151954604 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 191063676 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 225819916 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 256707276 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 284156028 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 308548524 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 330224588 characters processed.
[bwt_gen] Finished constructing BWT in 92 iterations.
[bwa_index] 59.04 seconds elapse.
[bwa_index] Update BWT... 0.57 sec
[bwa_index] Pack forward-only FASTA... 0.45 sec
[bwa_index] Construct SA from BWT and Occ... 37.95 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index DiNV-Dinn.fna
[main] Real time: 98.866 sec; CPU: 98.840 sec
```

Now, map the 16Cq reads to the combined DiNV-D. innubila genome. I am going to use multi-threading, and convert it to a sorted bam file in the same piece of code. I have to use the full path to the trimmed files for this. 

`bwa mem -t 5 DiNV-Dinn.fna /Users/maggieschedl/Desktop/KU/sequences/16Cq-DiNV-Test/mapping/16Cq-DiNV-R1-trim.fastq.gz /Users/maggieschedl/Desktop/KU/sequences/16Cq-DiNV-Test/mapping/16Cq-DiNV-R2-trim.fastq.gz | samtools view --threads 5 -h | samtools sort --threads 5 - > 16Cq-trim-h-mapped.bam`

`Real time: 961.361 sec; CPU: 4831.324 sec`

Now, map the KM3 reads to the combined DiNV-D. virilis genome, again I'll use the same code as above. 

`bwa mem -t 5 DiNV-Dvir.fna /Users/maggieschedl/Desktop/KU/sequences/16Cq-DiNV-Test/mapping/KM_3_1_trim.fq.gz /Users/maggieschedl/Desktop/KU/sequences/16Cq-DiNV-Test/mapping/KM_3_2_trim.fq.gz | samtools view --threads 5 -h | samtools sort --threads 5 - > KM3-trim-h-mapped.bam`

`Real time: 613.574 sec; CPU: 3071.495 sec`

Now I am going to use [Qualimap](http://qualimap.conesalab.org/doc_html/intro.html#installation) again to look at the mapping statistics. 

I want to run this on the combined mapping bam files as well as on just the reads that map to DiNV from these mappings. So first I have to separate out the DiNV-only-mapping reads from these bam files. 

Samtools has a function for this, but I first have to index the bam file:

`samtools index 16Cq-trim-h-mapped.bam`

Then samtools view can separate out the lines mapping to DiNV, and put it into a new file. The -b flag says to output as a bam file. 

`samtools view 16Cq-trim-h-mapped.bam NC_040699.1 -b > 16Cq-trim-h-mapped-DiNV-only.bam`

Now I want to do the same thing for the D. virilis cell reads:

```
samtools index KM3-trim-h-mapped.bam  
samtools view KM3-trim-h-mapped.bam NC_040699.1 -b > KM3-trim-h-mapped-DiNV-only.bam
```

Now I'll use Qualimap to get mapping statistics and information on the DiNV only and the compo genome mappings. 

Run as a loop to do all of the bam files (even though there are only 2), adding in -ip for "Activate this option to collect statistics of overlapping paired-end reads"

```
for bam in *.bam 
do
/Users/maggieschedl/Documents/programs/qualimap_v2.2.1/qualimap bamqc -bam $bam -ip -outfile ${bam}.pdf
done
```

