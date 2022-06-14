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
  - `conda create -n integrateEnv python=3.9`
  - `conda activate integrateEnv`
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
- Now it is having an issue with the first line, it cannot find a module named src.main
  - There is a folder called src, but it's not called src.main
  - Should I change it to say just src?
  - No, this does not work. It also says no module named src
- Ok now I checked the version of python now and it says 2.7.17! This is probably why it doesn't know what src.main is
- This happened when I did the pip install code... which it said I needed to do.
- I uninstalled the dependencies, then reinstalled them with conda
- My issue might be me checking the python version with "python" not "python3" and also not running the code that way
- So now if I run `python3 spacer_gen.py` I get an error that says there is no module named Bio in the main.py script
  - I think this might be one of the dependencies Biopython
  - Maybe using conda to install them was not right, I will retry the pip install code
  - But first I will uninstall them with conda
    - This ended up not working, it says they aren't there
    - So I just reran `pip install -r requirements.txt`
- Now I will try `python3 spacer_gen.py` again
- Same issue again, it just says No module named "Bio", I don't know what to do about this 
