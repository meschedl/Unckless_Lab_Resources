
## Sequence QC

I have 4 files for the 16Cq reads, presumably the data was from two lanes, so there is a lane 1 and lane 2 file for each forward and reverse. I can concatenate the forward and reverse files from each lane together. 

`cat 16Cq-DiNV_S232_L001_R1_001.fastq.gz 16Cq-DiNV_S232_L002_R1_001.fastq.gz > 16Cq-DiNV-R1.fastq.gz`
`cat 16Cq-DiNV_S232_L001_R2_001.fastq.gz 16Cq-DiNV_S232_L002_R2_001.fastq.gz > 16Cq-DiNV-R2.fastq.gz`

Use [fastp](https://github.com/OpenGene/fastp) to do read trimming and quality control of reads. 

download fastp to my local computer 

`conda install -c bioconda fastp`

Start with trimming the Dv-1 reads the way I had done them [previously](https://github.com/meschedl/DiNV-Dv1-Genome-Integration/blob/main/Data_QC/copy-data-and-fastp.md):

`fastp -i KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz -I KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz --trim_front1 15 --trim_front2 15 --html trim-fastp.html -o KM_3_1_trim.fq.gz -O KM_3_2_trim.fq.gz`

Output:

```
Read1 before filtering:
total reads: 7150325
total bases: 1072548750
Q20 bases: 1052709214(98.1502%)
Q30 bases: 1030796684(96.1072%)

Read2 before filtering:
total reads: 7150325
total bases: 1072548750
Q20 bases: 1032920148(96.3052%)
Q30 bases: 994062589(92.6823%)

Read1 after filtering:
total reads: 7026662
total bases: 812090153
Q20 bases: 799654639(98.4687%)
Q30 bases: 784993787(96.6634%)

Read2 after filtering:
total reads: 7026662
total bases: 812090153
Q20 bases: 787620681(96.9869%)
Q30 bases: 758952060(93.4566%)

Filtering result:
reads passed filter: 14053324
reads failed due to low quality: 247324
reads failed due to too many N: 2
reads failed due to too short: 0
reads with adapter trimmed: 7718326
bases trimmed due to adapters: 273565416

Duplication rate: 9.25084%

Insert size peak (evaluated by paired-end reads): 170

JSON report: fastp.json
HTML report: trim-fastp.html
```

Run default on the 16 Cq samples to see what they look like:

`fastp -i 16Cq-DiNV-R1.fastq.gz -I 16Cq-DiNV-R2.fastq.gz -o 16Cq-DiNV-R1-trim.fastq.gz -O 16Cq-DiNV-R2-trim.fastq.gz`

Output:

```
Read1 before filtering:
total reads: 12724286
total bases: 1921367186
Q20 bases: 1875546106(97.6152%)
Q30 bases: 1797984762(93.5784%)

Read2 before filtering:
total reads: 12724286
total bases: 1921367186
Q20 bases: 1866484094(97.1435%)
Q30 bases: 1772666316(92.2607%)

Read1 after filtering:
total reads: 12652278
total bases: 1841646916
Q20 bases: 1803328217(97.9193%)
Q30 bases: 1731845913(94.0379%)

Read2 after filtering:
total reads: 12652278
total bases: 1841634816
Q20 bases: 1795722864(97.507%)
Q30 bases: 1708962398(92.7959%)

Filtering result:
reads passed filter: 25304556
reads failed due to low quality: 137576
reads failed due to too many N: 5592
reads failed due to too short: 848
reads with adapter trimmed: 3827470
bases trimmed due to adapters: 135304308

Duplication rate: 5.21651%

Insert size peak (evaluated by paired-end reads): 237

JSON report: fastp.json
HTML report: fastp.html
```
Looked at the html file  INSERT LINK

Looks like the reads are really good, quality scores between 37 and 35. There is some variability in the QC content until about 12 bases in, so I think I will trim that off. Otherwise these should be good reads. 

Redo the fastp, renaming the html file to trim-fastp.html:

`fastp -i 16Cq-DiNV-R1.fastq.gz -I 16Cq-DiNV-R2.fastq.gz --trim_front1 15 --trim_front2 15 --html trim-fastp.html -o 16Cq-DiNV-R1-trim.fastq.gz -O 16Cq-DiNV-R2-trim.fastq.gz`

INSERT LINK TO FASTP HERE




