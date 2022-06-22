# INTEGRATE Guide RNA Tool Design

Using [INTEGRATE Guide RNA Tool](https://github.com/sternberglab/INTEGRATE-guide-RNA-tool) and scripts to find Type I-F Cas9 guide RNAs for cloning DiNV


- Forked repository then cloned it to my desktop computer
- Going to edit it on the desktop, then copy it to the linux to run
- Pretty much the only thing I did was edit the spacer_eval.py
  - I set region type to custom
  - And generated a csv file that has the name of the genome (GCF_004132165.1_ASM413216v1_genomic.fna)
  - And the start and stop points: 76741 - 78478 bases
  - This should be the 77kb region of the DiNV genome
- I also said to output the regions in the same folder, generate 1000 spacers per region (suggested number), and have GC content between 35-65%
- That should be all I want to modify right now, so I will copy it all to the Linux
`scp /Users/m741s365/Desktop/Github/INTEGRATE-guide-RNA-tool/* runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/INTEGRATE-gRNA`
- That did not bring along the src directory, so I have to also use the -r
`scp -r /Users/m741s365/Desktop/Github/INTEGRATE-guide-RNA-tool/src runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/INTEGRATE-gRNA`
- Now I need to copy over the DiNV genome
`cp /home/runcklesslab/Maggie/DiNV-DV-1/GCF_004132165.1_ASM413216v1_genomic.fna /home/runcklesslab/Maggie/INTEGRATE-gRNA`
- I created a new conda environment with python version 3.9 to be able to run these python scripts (specified on the Github page, needs to be 3.8 or higher)
  - `conda create -n IntegrateEnv python=3.9`
  - `conda activate IntegrateEnv`
- Also need to copy in the custom_regions.csv from the desktop
`scp /Users/m741s365/Desktop/custom_regions.csv runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/INTEGRATE-gRNA`
- Tool Github says that I can run this command to make sure dependencies are installed
`pip install -r requirements.txt`
- This seemed to have installed things, however it gave a warning at the beginning about the depreciation of python 2.7...
- Ok, I am just going to run the spacer_gen.py wiht python and see what happens
`python spacer_gen.py`
- Ok this immediately gave me an error saying that it didn't like the blank space after email = on line 43
  - I had gone and left every section blank for the genebank ID and genebank file sections because I am only doing the custom setting
  - Seems like that wasn't right
- I have gone back and made every bank into something = []
- I will remove the spacer_gen.py, then copy and paste what is in my desktop version of the file to the Linux when I nano spacer_gen.py (new file)

**Some more errors occured, until we figured out that installing the requirements with the text method does not actually install them, and that to run the python scripts you have to specify python3**


### code to run this program locally on my laptop
- All code below is done in a new environment  
`conda create -n INTEGRATEenv python=3.9`  
`conda activate INTEGRATEenv`
- Make sure bowtie2 is installed in the PATH  
`conda install -c bioconda/label/broken bowtie2` #works
- List installed programs (make sure bowtie2 is included)
`conda list`
- Show the location of installed programs in the current environment
`ls -lh /opt/anaconda3/envs/INTEGRATEenv/`
- Add bowtie2 to your path
`sudo nano /etc/paths`
- (need to use password to computer here)
- Add the following line:
`/opt/anaconda3/envs/INTEGRATEenv/`
- Check that it worked
`echo $PATH`
- Open up a new terminal window to allow the changes to take effect
- Now you should be able to run the following without an error:
`bowtie2 -h`
- Dependencies
- Install biopython individually  
`pip install biopython==1.76`
- Install simplesam individually   
`pip install simplesam==0.1.3`
- There is another program called boto3 that needs to be installed (version unknown)
`pip install boto3`
- You have to rename your reference genome as a single letter (no idea why, this program is garbage)
`cp GCF_004132165.1_ASM413216v1_genomic.fna d`
- Make new directory to hold indexed reference genome
`mkdir /Users/maggieschedl/Desktop/KU/INTEGRATE-gRNA/assets/bowtie/DiNV`
- Run this line to index the reference genome specifying 'd' as input reference genome and a path to output the index files since the program still won't recognize the reference genome
`bowtie2-build d /Users/maggieschedl/Desktop/KU/INTEGRATE-gRNA/assets/bowtie/DiNV/index`
- Now finally run the program
`python spacer_gen.py`
- Look in the output folder for your identified spacers with minimal off-target binding potential
