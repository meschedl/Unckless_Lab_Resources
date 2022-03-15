## Confirming BAC Sequences from Plasmidsaurus

- Make new directories in Linux for Plasmidsaurus and BACs
- Copy in data from Plasmidsaurus to Linux machine
- Copy the two complete fasta sequences:  
`scp /Users/m741s365/Desktop/Lab/Plasmidsaurus/Unckless_ed_results/*.fasta runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Plasmidsaurus/BACs`
- Copy in the raw reads:  
`scp /Users/m741s365/Desktop/Lab/Plasmidsaurus/Unckless_ed_raw_reads/*.fastq runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Plasmidsaurus/BACs`
- Then I took the pBeloBAC11 and pBACe3.6 sequences from Geneious (downloaded from NCBI) and exported them to my desktop in txt format. These are in a sort of strange format, there is an enter/extra blank line between each line, and a space every 10 bases. I can easily remove the lines by quickly deleting those, but it will be too tedious by hand to remove the spaces. I'll have to do that computationally
- So I need to copy the txt files to the Linux:   
`scp /Users/m741s365/Desktop/Geneious/*.txt runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Plasmidsaurus/BACs`  
- Then I want to get rid of the spaces by using sed. Sed -i "updates" the file, and whatever is between the // get's found in the first one, and replaced by the second one   
`sed -i 's/ //g' CVU51113_pBeloBAC11.txt`  
`sed -i 's/ //g' CVU80929_pBACe3.6.txt`
- I used `wc -m` to get the number of characters in each file
  - Both the NCBI and the Plasmidsaurus pBeloBAC11s had 7555 characters
  - The NCBI pBACe3.6 has 11657 characters
  - And the Plasmidsaurus pBACEe3.6 has 11686 characters
  - Characters should also mean bases
  - These are either exactly the same or very close
- Because these files don't necessarily start and stop in the same part of the plasmid, I can't use a simple command like `cmp` or `diff` to compare the differences between the NCBI and the Plasmidsaurus files
- I wonder if I can try mapping the Plasmidsaurs file to the NCBI one, and hopefully get only 1 alignment?
- First try adding on a header line to these text files  
`nano CVU80929_pBACe3.6.txt`  
and added `> NCBI pBACe3.6 CVU80929`  
`nano CVU51113_pBeloBAC11.txt`  
and added `> NCBI pBeloBAC11 CVU51113`
- Let's see if I can get bwa to index my NCBI files   
`bwa index CVU80929_pBACe3.6.txt`    
```
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index CVU80929_pBACe3.6.txt
[main] Real time: 4.395 sec; CPU: 0.006 sec
```  
`bwa index CVU51113_pBeloBAC11.txt`  
```
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index CVU51113_pBeloBAC11.txt
[main] Real time: 3.572 sec; CPU: 0.006 sec
```   
- Now try mapping the full sequence from Plasmidsarus to the NCBI file   
`bwa mem CVU80929_pBACe3.6.txt Unckless_ed_1_pBACe3.6.fasta | samtools view -h -b | samtools sort - > pBACe3.6.bam`  
```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 1 sequences (11612 bp)...
[M::mem_process_seqs] Processed 1 reads in 0.012 CPU sec, 0.041 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem CVU80929_pBACe3.6.txt Unckless_ed_1_pBACe3.6.fasta
[main] Real time: 0.056 sec; CPU: 0.014 sec
```  
`bwa mem CVU51113_pBeloBAC11.txt Unckless_ed_2_pBeloBAC11.fasta | samtools view -h -b | samtools sort - > pBeloBAC11.bam`   
```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 1 sequences (7507 bp)...
[M::mem_process_seqs] Processed 1 reads in 0.001 CPU sec, 0.001 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem CVU51113_pBeloBAC11.txt Unckless_ed_2_pBeloBAC11.fasta
[main] Real time: 0.002 sec; CPU: 0.002 sec
```
- Ok, now to look at these  
`samtools view pBeloBAC11.bam | less`
- There are just 2 alignments for pBeloBAC11  
`samtools view pBACe3.6.bam | less`
- There are also only 2 alignments for pBACe3.6
- **Investigating pBeloBAC11 further:** looking at the CIGAR string on the alignments to see how many bases are a match and how many don't match. This is with H for hard clipped (not present in mapping) and M for matched
- The 2 alignments for pBeloBAC11 have reciprocal CIGAR strings:
  - 4743H2764M
  - 4743M2764S
- I am pretty sure this means that the entire consensus sequence from Plasmidsaurus maps exactly to the sequence from NCBI, it just has different start and stop points
- 4743 + 2764 = 7507 bases
- This is actually less than the number of bases that I got with the `wc -m` above, which is 7555. Seems like I am missing ~50 bases
- I also might have counted the header lines when doing `wc -m`, so I copied and pasted just the bases into an [online character counter](https://lettercounter.github.io/)
  - Unckless_ed_2_pBeloBAC11.fasta 5507 characters, no spaces
  - CVU51113_pBeloBAC11.txt 5507 characters without spaces, 5553 characters **with** spaces
- So, even though it doesn't look like there are spaces in CVU51113_pBeloBAC11.txt, there still are some (although I really thought I got rid of all of them). So, this looks like a perfect mapping for pBeloBAC11!
- **Investigating pBACe3.6 further:** looking at the CIGAR string
- The 2 alignments for pBACe3.6 are also reciprocal:
  - 2452S9160M
  - 2452M9160H
- 2452 + 9160 = 11,612 bases
- I want to check the actual bases character count for pBACe3.6, same way as above
  - CVU80929_pBACe3.6.txt 11,612 characters without spaces, 11,684 characters **with** spaces
  - Unckless_ed_1_pBACe3.6.fasta 11,612 characters
- So I would say pBACe3.6 is also an exact map!   
