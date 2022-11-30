
## Aligning Consensus Sequence from Plasmidsaurus Sequencing of pSL1142 pSPIN onto the Addgene Sequence 

This is to assess the accuracy of the sequencing/plasmid. I may also align the raw reads. 

- Saved the [addgene sequence](https://www.addgene.org/160730/sequences/) as pSL1141.fa
- Download bwa to my laptop from [bwa website](https://github.com/lh3/bwa), this was put in my user directory 

```
git clone https://github.com/lh3/bwa.git
cd bwa; make
```
- Add bwa to my path 
```
nano ~/.zshrc
export PATH=$PATH:/Users/maggieschedl/bwa #add this into the "document" and write out and exit
source ~/.zshrc
```
- Download samtools and put in my user directory 
```
cd samtools-1.16.1
make
sudo make install #I had to put in my password here
#then add it to the path
export PATH=/Users/maggieschedl/samtools-1.16.1:$PATH
```

- Index the addgene sequence as the "reference"
```
cd Desktop/
bwa index  pSL1142.fa.txt
```
- Align plasmidsaurus consensus to the reference from addgene, and make this into a bam file for safe storage
```
bwa mem pSL1142.fa.txt /Users/maggieschedl/Desktop/KU/Plasmidsaurus_results/Unckless_k5n_results/Unckless_k5n_1_pSL1142.fastq | samtools view -h -b | samtools sort - > pSL1142.bam
```
- Output from bwa mem
```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 1 sequences (14032 bp)...
[M::mem_process_seqs] Processed 1 reads in 0.008 CPU sec, 0.008 real sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa mem pSL1142.fa.txt /Users/maggieschedl/Desktop/KU/Plasmidsaurus_results/Unckless_k5n_results/Unckless_k5n_1_pSL1142.fastq
[main] Real time: 0.009 sec; CPU: 0.010 sec
```
- Look at the alignments 
```
samtools view pSL1142.bam
```
- This seems pretty good, there are only 2 alignments, but I want to look at the CIGAR strings to be sure 
    - 3604S3103M1I7324M
    - 3604M10428H  
- This seems to me like it's saying 3,103 bases match, then a 1 base insertion, then 7,324 bases match, then 3,604 bases match. This is a total of 14,031 bases matching with a 1 base insertion. This is the exact number of bases on the Addgene listing of the plasmid, so I would say we are good! 
