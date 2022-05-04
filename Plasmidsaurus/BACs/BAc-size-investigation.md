## Investigating Plasmidsarus Reads Further Based on Size


I want to re-check the read sizes from Plasmidsarus to see if it's possible the second band in pBeloBAC11 is likely to be a doublet/triplet or a second construct.  
They gave us little graphs of the read lengths:  

pBACe3.6:
![](https://raw.githubusercontent.com/meschedl/Unckless_Lab_Resources/main/images/Unckless_ed_1_pBACe3.6.png)  

pBeloBAC11:
![](https://raw.githubusercontent.com/meschedl/Unckless_Lab_Resources/main/images/Unckless_ed_2_pBeloBAC11.png)

**Get list of too long reads**

Sizes are supposed to be 11612bp for pBACe3.6 and 7,507bp for pBeloBAC11. Based on the graphs, there are larger reads in both sequencing runs, but especially for pBeloBAC11. Based on the x-axis of the graphs, not all the larger reads look like they are doublets or triplets just on exact size, however it's hard to know if if the DNA was broken or not. In Geneious I'm able to look at the exact length of each of the reads.

pBACe3.6 reads longer than expected (11612bp)
- 28,408bp read 10
- 20,493bp read 111
- 16,297bp read 197
- 58,936bp read 329
- 21,080bp read 380
- 39,623bp read 490
- 30,718bp read 609
- 22,490bp read 708
- 23,099bp read 1001

pBeloBAC11 reads longer than expected (7,507bp)
- 8,820bp read 8
- 11,033bp read 74
- 10,851bp read 100
- 25,321bp read 108
- 12,108bp read 113
- 17,443bp read 123
- 9,731bp read 162
- 31,523bp read 174
- 12,739bp read 285

Pretty much none of these seem like doublets except for read 1001 for pBACe3.6, because that's close to exactly twice the expected size 11612. Both these samples have the same number of longer reads, but there are 3 times as many 3.6 reads as there are 11 reads, so it's a much larger proportion of reads for the 11 sample. I think the best thing to do is try to map the reads to the reference sequences, and then check these reads for how they map - basically look for overlapping mapping or look to see if there are any sections that don't map. If they don't map then maybe I BLAST them?

**Map long reads to NCBI reference**

In the Linux:  
`cd /home/runcklesslab/Maggie/Plasmidsaurus/BACs`  
`mkdir read_mapping`  
I have the reference sequences from NCBI already on the Linux, I want to copy them to this folder, as well as the associated files that were made as a bwa index.  
`cp CVU51113* /home/runcklesslab/Maggie/Plasmidsaurus/BACs/read_mapping/`  
`cp CVU80929* /home/runcklesslab/Maggie/Plasmidsaurus/BACs/read_mapping/`  
Then I need to copy over the raw reads that are in the SPAdes directory.  
`cd SPAdes`  
`cp *fastq /home/runcklesslab/Maggie/Plasmidsaurus/BACs/read_mapping/`  
Ok now I should just be able to map the raw reads to the NCBI reference.
`bwa mem CVU80929_pBACe3.6.txt Unckless_ed_1_pBACe3.6.fastq | samtools view -h -b | samtools sort - > pBACe3.6.raw.bam`  
`bwa mem CVU51113_pBeloBAC11.txt Unckless_ed_2_pBeloBAC11.fastq | samtools view -h -b | samtools sort - > pBeloBAC11.raw.bam`

There are a lot of lines in this file so I am just going to pull out the reads that have the long lengths and look at them.   
28,408bp read 10  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_10`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
20,493bp read 111  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_111`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string   
16,297bp read 197  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_197`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
58,936bp read 329  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_329`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
21,080bp read 380  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_380`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
39,623bp read 490  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_490`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
30,718bp read 609  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_609`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
22,490bp read 708  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_708`  
This sequence maps! Actually it maps 5 times. But the CIGAR strings are almost unreadable, this is one example:  
```
7M1D125M1D44M1I31M4D286M1D32M1D38M1D36M1I23M1D21M2D18M1D100M1I149M1I120M1I71M1I110M1I29M1D63M1I8M1I126M1D11M1D73M3D62M4D73M1I12M1D121M2D37M2I84M1D124M5I35M1D75M2I2M1D28M1D10M1D49M1D43M1D15M1I66M1D121M1D78M3I13M2D127M2D54M1D190M1D75M1I85M3I137M1D5M2D13M1I46M2D105M1I225M1D119M3I44M1I356M5D14M1I48M1D43M1D9M3D59M9D3M8I7M2I13M1D15M1D106M1I144M1D182M1D25M2D147M3D38M1D155M6D5M1D4M1D14M1D2M3I35M2D5M1D87M1D22M2D146M4D3M2I67M2D60M2D8M1D26M2D48M1D73M2I40M2D150M1D16M2D83M1I116M1D138M1D251M1I3M1I22M1D43M1D48M1I153M1D355M2D48M1I3M1D320M1I64M1D27M1I159M1D167M2D71M1I39M3D36M4D77M3D27M2I151M1I62M1I26M3I34M1D11M1D6M1D32M1I107M1D70M1D62M1I75M1I16M1I360M1I217M2D28M1I340M1D100M1I159M1D83M2D84M1I262M1I123M1D69M2D212M11506S
```
A lot of the reads that map only once have CIGAR strings that look like this. My guess is that this one might be a doublet, it's close the the right size to be one. And if it maps multiple times that seems to make sense for a doublet.    
23,099bp read 1001  
`samtools view pBACe3.6.raw.bam | grep Unckless_ed_1_pBACe3.6_1001`  
This sequence also maps multiple times, 4 times. Again, the CIGAR strings are hard to make anything out of:  
```
82M1D242M5D69M1D40M1D72M2I3M1D38M2I99M1I197M3D155M1D140M1I191M1D7M1I182M1D92M1I30M1D135M1I115M1D227M1I19M1I68M1I23M1I4M3I100M5D134M1D9M2D8M1D57M2I16M1D17M2D6M1I7M1D33M1I32M1I142M6D113M3D203M5I3M3D9M1I4M1D71M4D46M1D20M2I11M3D7M3D1M1I45M11D104M2D29M2D52M1D83M1D25M1I83M2I10M3D1M1I75M3D2M1I6M2D194M1D5M2D43M3I1M2D10M2D196M1I109M1D21M1I80M2I42M4D66M1D20M1D81M3I19M1D124M2I2M2D16M1D8M2D10M1I35M4D131M1D6M2D146M2D37M2D4M2D44M2I41M2I2M1D87M2D25M1D99M1D5M1I8M2D26M1D3M1I111M1D10M1I37M1I4M3D9M1D8M6D10M3I37M5D10M4I24M3D45M4I1M4D20M1D133M2D21M2D12M1I14M2D293M1I1M2D15M1I55M1D17M3I4M3D144M1D43M4D65M2D57M1D5M2D14M1D253M1D11M7D6M3I4M1D73M1I6M4D155M2I81M2D20M1D5M1I29M1I60M1D6M2I14M1D36M5D26M1D16M1I52M2D10M1D45M1D115M1I11M1I1M3D9M1D168M5D121M1I18M1D39M1D197M1I117M2D156M2D56M3D7M1I60M1I73M2I163M2D81M13948H
```
Again, this is one where the size of the read is in doublet range, and having multiple mappings it seems like it probably has sequence all along it that map.

Check the pBeloBAC11 long reads and their mapping:   
8,820bp read 8  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_8`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
11,033bp read 74  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_74`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
10,851bp read 100  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_100`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
25,321bp read 108  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_108`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
12,108bp read 113  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_113`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
17,443bp read 123  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_123`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
9,731bp read 162  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_162`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
31,523bp read 174  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_174`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  
12,739bp read 285  
`samtools view pBeloBAC11.raw.bam | grep Unckless_ed_2_pBeloBAC11_285`  
This sequence didn't map (bit flag 4), and therefore has no CIGAR string  


This is interesting, I was thinking that it was more likely that pBeloBAC11 had a doublet but it looks like all of the too long reads are not mapping at all to the BAC sequence. I wonder if this is _E. coli_ DNA? I did have to grow up A LOT of bacteria to do these extractions. I'm going to BLAST the long sequences and see what comes up as the top hit.

**BLAST non mapping reads**

Doing the grep for the reads in samtools (above) prints out the read for me. I am not putting any parameters in BLAST other than I'm looking for a highly similar sequence.  
Start with pBACe3.6.  
- 28,408bp read 10  - all top BLAST hits are to _E. coli_ strains
- 20,493bp read 111 - all top BLAST hits are to _E. coli_ strains and one Salmonella
- 16,297bp read 197 - all top BLAST hits are to _E. coli_ strains
- 58,936bp read 329 - all top BLAST hits are to _E. coli_ strains
- 21,080bp read 380 - all top BLAST hits are to _E. coli_ strains and one Salmonella
- 39,623bp read 490 - all top BLAST hits are to _E. coli_ strains
- 30,718bp read 609 - all top BLAST hits are to _E. coli_ strains
- 22,490bp read 708 - top BLAST hit is pBACe3.6!
- 23,099bp read 1001 - top BLAST hit is pBACe3.6!      

pBeloBAC11
- 8,820bp read 8 - all top BLAST hits are to _E. coli_ strains
- 11,033bp read 74 - all top BLAST hits are to _E. coli_ strains
- 10,851bp read 100 - all top BLAST hits are to _E. coli_ strains
- 25,321bp read 108 - all top BLAST hits are to _E. coli_ strains
- 12,108bp read 113 - all top BLAST hits are to _E. coli_ strains and one shigella (bacteria)
- 17,443bp read 123 - all top BLAST hits are to _E. coli_ strains
- 9,731bp read 162 - all top BLAST hits are to _E. coli_ strains and one Salmonella
- 31,523bp read 174 - all top BLAST hits are to _E. coli_ strains
- 12,739bp read 285 - all top BLAST hits are to _E. coli_ strains and one shigella (bacteria)

Very surprising! Based on the gel I thought that it was pBeloBAC11 that had the doublet, but it just looks like it is all contamination in the extraction from the bacteria. And that is in both of them. And it looks like there are two reads that hint at a possible doublet in pBACe3.6.
