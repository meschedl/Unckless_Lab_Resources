# Various Methods for Assembling Plasmidsaurus Raw Data into Consensus

**Using Geneious to Assemble Raw Reads from Plasmidsaurus and Mapping to BAC "Genomes"**

- Loaded in the raw reads to Geneious and look at them under sequence view
- Looking at them in sequence view you can see the spread of the data. It actually looks like pBeloBAC11 sequenced more poorly than pBACe3.6. There are a lot of short reads for pBeloBAC11, and there are way less reads too ~300 compared to ~1000
- pBeloBAC11:
![](https://raw.githubusercontent.com/meschedl/Unckless_Lab_Resources/main/images/pBeloBAC11_raw.png)
- pBACe3.6
![](https://raw.githubusercontent.com/meschedl/Unckless_Lab_Resources/main/images/pBACe3.6_raw.png)
- This is another sign that we may not want to use pBeloBAC11 as out BAC for the project
- I'm using the [De novo Assembly Tutorial](https://www.geneious.com/tutorials/de-novo-assembly/) from Geneious as a guide for this process
- I started with pBACe3.6
- Installed BBDuk trimming program
- I tried trimming the reads based off of what the tutorial said, but that ended up shortening a lot of my sequences. This is probably because nanopore sequencing is pretty error prone, so the 20 base quality score is probably way higher than what you normally get with these sequences. So I decided just to "trim" the reads by removing any reads shorter than 20bp and lowered the base quality score to 5
- There weren't any reads that got cut out with this threshold, so it didn't really do anything
- Went through running assembly with default conditions and it gave me something 10X as long as it should be, I don't think this assembler thought the consensus could/should be so short like a plasmid. Even running with an option for "circularization of contigs" doesn't end up giving me anything close, this time it gave me a "consensus" that was +700,000bp long! It should be 11,000 for pBACe3.6.
- At this point I stopped and I'm going to try a different program

**Using [SPAdes](https://github.com/ablab/spades) to Assemble Raw Reads into "Genomes"**

- Using [SPAdes manual](https://github.com/ablab/spades#sec3) for help
- This program will take nanopore reads and also has an option for plasmids
- General code is: `spades.py [options] -o <output_dir>`
- Options I think I need to use for this:
  - `--nanopore`
  - `--careful` "Tries to reduce the number of mismatches and short indels, only for small genomes"
  - `-s <file_name>` "File with unpaired reads"
  - or `--nanopore <file_name>`
  - `--plasmid` "This flag is required when assembling only plasmids from WGS data sets" hopefully this can take nanopore data
- Default number of threads is 16
- Default amount of memory is 250GB (usually not reached)
- Raw reads from Plasmidsaurus are in the Linux at /Maggie/Plasmidsaurus/BACs/ Unckless_ed_2_pBeloBAC11.fastq and Unckless_ed_1_pBACe3.6.fastq
- There are other files in this directory, so I'm going to make a SPAdes directory for this work  
`mkdir SPAdes`
- Then I'm going to move those files into that directory, they aren't needed in the one above  
`mv Unckless_ed_2_pBeloBAC11.fastq SPAdes/`  
`mv Unckless_ed_1_pBACe3.6.fastq SPAdes/`  
`cd SPAdes/`
- I need to create an output directory for the program  
`mkdir 3.6_SPAdes_out`
- Try running program?  
`spades.py --plasmid --careful --nanopore Unckless_ed_1_pBACe3.6.fastq -o 3.6SPAdes_out`
- Nope, did not like "== Error ==  you should specify at least one unpaired, paired-end, or high-quality mate-pairs library!"
- So I have to try using -s  
`spades.py --plasmid --careful --nanopore -s Unckless_ed_1_pBACe3.6.fastq -o 3.6SPAdes_out`
- Also another error "spades.py: error: argument --nanopore: expected 1 argument(s)"
- I think this means it needs the name of the file of the nanopore reads after the --nanopore. Maybe I can get away with putting in the file name twice??  
`spades.py --plasmid --careful --nanopore Unckless_ed_1_pBACe3.6.fastq -s Unckless_ed_1_pBACe3.6.fastq -o 3.6SPAdes_out`
- This does not work either "== Error ==  file /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/Unckless_ed_1_pBACe3.6.fastq was specified at least twice"
- Ok I will try without specifying nanopore but this may not work right  
`nohup spades.py --plasmid --careful -s Unckless_ed_1_pBACe3.6.fastq -o 3.6SPAdes_out`
- This ended up working   

```
* Corrected reads are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/corrected/
 * Assembled contigs are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/contigs.fasta
 * Assembled scaffolds are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/scaffolds.paths
 * Assembly graph is in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/assembly_graph.fastg
 * Assembly graph in GFA format is in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/assembly_graph_with_scaffolds.gfa
 ```
 - I looked in the contigs.fasta file and it actually looks pretty short, this might have worked better? I want to know the number of bases
 - So there are 6 "contigs", but they don't add up quite right to the size that pBACe3.6 should be, which is 11,612 bases
 - I copy pasted the contig bases into a [character counter](https://wordcounter.net/character-count)
  - Contig 1 is 11,660 bases
  - Contig 2 is 1,206 bases
  - Contig 3 is 247 bases
  - Contig 4 is 259 bases
  - Contig 5 is 252 bases
  - Contig 6 is 249 bases
- I'm not for sure, but contig 1 is basically the length of what the BAC should be, so that might be as good as I can get for a consensus sequence
- I also tried running this for pBeloBAC11
`mv nohup.out 3.6nohup.out`  
`mkdir 11SPAdes_out`  
`nohup spades.py --plasmid --careful -s Unckless_ed_2_pBeloBAC11.fastq -o 11SPAdes_out`
- This ended up not working  

```
* Corrected reads are in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/11SPAdes_out/corrected/  
No plasmid contigs assembled!!  
 * Assembly graph is in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/11SPAdes_out/assembly_graph.fastg  
 * Assembly graph in GFA format is in /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/11SPAdes_out/assembly_graph_with_scaffolds.gfa
 ```
 - I'm not sure if this is because there are so many short reads in pBeloBAC11, there are so many less reads in pBeloBAC11, or because I couldn't specify nanopore it didn't work here but was somehow able to grind through it in pBACe3.6


At least for pBACe3.6, I am going to take contig 1 and map that to the pBACe3.6 genome from NCBI and see what I get.

- In 3.6SPAdes_out directory:  
`mkdir Map`  
`cd Map/`  
- Make a file with contig 1 from the 3.6 assembly  
`nano pBACe3.6_SPAdes_contig_1.txt`
```
> NODE_1_length_10408_cov_268.384226_component_0
CGGATCAATTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCA
CAGTTAAATTGCTAACGCAGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTC
.....
```
- Link in the sequence from NCBI   
`ln -s /home/runcklesslab/Maggie/Plasmidsaurus/BACs/CVU80929_pBACe3.6.txt /home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/Map/`
- Index the NCBI sequence  
`bwa index CVU80929_pBACe3.6.txt`
- Then map  
`bwa mem CVU80929_pBACe3.6.txt pBACe3.6_SPAdes_contig_1.txt | samtools view -h | samtools sort - > SPAdes3.6contig1.bam`

```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 1 sequences (10408 bp)...
[M::mem_process_seqs] Processed 1 reads in 0.008 CPU sec, 0.121 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem CVU80929_pBACe3.6.txt pBACe3.6_SPAdes_contig_1.txt
[main] Real time: 0.121 sec; CPU: 0.010 sec
[W::sam_parse1] empty query name
[W::sam_parse1] empty query name
[W::sam_parse1] empty query name
[W::sam_parse1] empty query name
```
- Look at the alignment summary  
`samtools view SPAdes3.6contig1.bam`
- There are 2 alignments:
  - 5591H773M1D310M1D3734M
  - 5591M4817S
- So, how much of this did map? 5591 + 773 + 310 + 3734 = 10,408 bases
- This isn't the whole sequence, but it is weird to me that the program says `[M::process] read 1 sequences (10408 bp)`, this makes it seem like it's saying that pBACe3.6_SPAdes_contig_1.txt is only that many bp long, which it isn't. I used 2 different character counters to check the number of characters and its 11,660
- If I add up the second alignment 5591 + 4817 = 10,408 bases
- I'm not sure what is going on here, why doesn't bwa see those other bases?  
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Plasmidsaurus/BACs/SPAdes/3.6SPAdes_out/Map/pBACe3.6_SPAdes_contig_1.txt /Users/m741s365/Desktop/Github/Unckless_Lab_Resources/Plasmidsaurus/BACs/ `
- Ok so this is not actually right, there are only 10,408 bases, the character counters I was using were detecting the "enters" as letters. So contig 1 is not actually a great consensus sequence, but it is close
- Maybe I should try mapping all the contigs to the NCBI reference?  
`cp contigs.fasta ./Map/`
- Map all contigs to NCBI reference  
`bwa mem CVU80929_pBACe3.6.txt contigs.fasta | samtools view -h | samtools sort - > SPAdes3.6.bam`

```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 6 sequences (12621 bp)...
[M::mem_process_seqs] Processed 6 reads in 0.009 CPU sec, 0.009 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem CVU80929_pBACe3.6.txt contigs.fasta
[main] Real time: 0.010 sec; CPU: 0.011 sec
```
- Now there are 7 alignments
- The same 2 alignments for contig 1
  - 5591H773M1D310M1D3734M
  - 5591M4817S
- As well as ones for the other contigs
  - 49M198S contig 3
  - 9S84M166S contig 4
  - 1120M86S contig 2
  - 44M208S contig 5
  - 122M1D127M contig 6
- Contig 1 (10,408 bases) is missing 1,254 bases to get to the NCBI sequence at 11,662 bases
- I'm not entirely sure what to make of the alignment, but they align up pretty good, but it's not exact. It's hard to say because I know the assembly didn't go "the right way" so should I be expecting this to align up perfectly anyways?
- I also cannot do anything with pBeloBAC11 because it din't assemble...
