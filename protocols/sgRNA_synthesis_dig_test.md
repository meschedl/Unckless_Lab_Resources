# sgRNA Synthesis and Digestion Test Protocol


**This protocol is for testing sgRNAs for use in a CRISPR experiment. There are multiple steps:**
1. 50ul PCR of DNA region to be test digested
2. Bead cleanup and quantification of PCR product
3. sgRNA synthesis, purification, and quantification
4. Digestion test with Cas9 and gel

### Step 1 50ul PCR reaction of DNA region
- This will be based off your own optimized protocol for the PCR you are doing
- To make it a 50ul reaction, multiply every component of the reaction by 5 (assuming you start with a 10ul reaction). This includes the amount of DNA you put in
- **It is very important that you also change the protocol on the PCR machine to say the volume is 50ul when you run it**
- After the PCR, take 5ul of the reactions and run on a gel to confirm amplification as normal, but make sure to save the remaining 45ul in the fridge until you are ready to go to the next step

### Step 2 Bead clean and quantify PCR product
**Bead Clean**
- Before using in the digestion test, the enzymes and primers need to be cleaned out of you PCR product tubes
- Take the Ampure XP beads out of the fridge and place in a drawer for ~30 minutes to let come to room temperature. They are slightly light sensitive so that is why you keep them in a drawer
- The beads (brown) settle out of solution in the fridge, and thus need to be resuspended. **Swirl** the bottle of beads gently until the whole liquid is brown and you cannot see concentrated beads on any of the sides of the bottle. **Do not vortex the beads**
- Take PCR products out of the fridge and put on bench at room temp (you can remove the negative and/or positive controls, they aren't needed)
- Once beads are at room temp, add either 1X or 2X volume to the PCR product tube(s)
  - For PCR products larger than 100bp, use 1X (or 45ul)
  - For PCR products smaller than 100bp, use 2X (or 90ul)
- Pipette mix the beads into the liquid gently by pipetting up and down at least 10 times. You will see the liquid become homogeneously brown tinted (although this can be hard to see because the green in the PCR mix)
- Place the bead-PCR tubes in the drawer for 15 minutes
  - This allows the DNA time to adhere to the beads. Under certain pH DNA will bind to the beads but all else will not. These beads can be washed, and then the DNA can be eluted in a purified form.
- While you are waiting, make fresh 80% ethanol using 100% ethanol and molecular grade water. You will need ~400ul per tube.
- Get the 96 well magnet plate out of the drawer where the qubit supplies are
- After the 15 minutes, place the tubes on the magnet plate and wait ~2 minutes
  - You should be able to see that the beads have gone to the side of the tube against the magnet. There should be a clear brown spot where the beads are, and the liquid should be clear (but green). It is very important to make sure all the beads have gone to the magnet
- Keep tubes on magnet. Pipette the liquid off without disturbing the beads: to do this set your pipette to ~7ul **less** than the volume in the tube. Do not worry about the extra volume, you just want to avoid removing any beads
- Keep tubes on magnet. Add 200ul of fresh 80% ethanol to each tube
- Wait ~30 seconds
- Keep tubes on magnet. Remove 200ul liquid from the tubes
- Keep tubes on magnet. Add 200ul of fresh 80% ethanol to each tube
- Keep tubes on magnet. Remove 200ul liquid from the tubes
- Keep tubes on magnet. Use a small pipette to remove any excess liquid from the tubes, try to remove any droplets
- Let the tubes sit open for ~2 minutes to let any excess ethanol dry, but not too long because over-dried beads will not release DNA
- Take the tubes off the magnet
- Resuspend the beads in 30ul of molecular grade water. Pipette them off the side of the tube and pipette mix the liquid multiple times gently to make the solution homogenously brown
- Let the tubes sit in the drawer for 5 minutes
- Label a new set of tubes
- After the 5 minutes, place the tubes on the magnet again and wait until the solution goes clear
- Remove 29ul of liquid from the tubes to their new tubes
- Keep those tubes on ice while quantifying

**Quantify**
- Use the Broad Range dsDNA Qubit Kit (includes a large bottle and a small tube)
- Take the BR dsRNA standards out of the fridge and let come to room temperature
- Determine _n_ number: this is (# of samples to qubit + 2) * 1.1
  - Example: (5 samples + 2) * 1.1 = 7.7
- Multiply 199 by your n number
  - Example: 199 * 7.7 = 1532.3
- Add that amount in ul (ex 1532.3ul) of buffer (big bottle) to a new tube
- Add _n_ ul of quant IT reagent (small tube) to the same tube
- Vortex and spin down, this is your working solution
- Set up qubit tubes (in the same drawer as the qubit kit and machine)
  - You will need 2 tubes for the standards, and 1 tube per sample
- Label standard tubes S1 and S2, and sample tubes by sample names
- Add 190ul of working solution to each standard tube
- Add 199ul of working solution to each sample tube
- Add 10ul of standard 1 to the S1 tube
- Add 10ul of standard 2 to the S2 tube
- Add 1ul of cleaned DNA sample to their respective tubes
- Vortex all tubes and quick spin them down
- Place them in a drawer for 2 minutes
- Plug in the Qubit and select dsDNA broad range
- Tell it to read new standards
- When prompted, add the S1 tube, close lid, and follow instructions
- When prompted, add the S2 tube, close lid, and follow instructions
- When prompted, add sample tube, close lid, and follow instructions
- The reading will come up, but you will need to toggle to "calculate sample concentration"
- Tell the Qubit you used 1ul input DNA
- Then it will say the concentration of your sample in ug/mL which is the same as ng/ul
- Add in your next sample tube and repeat for all samples
- When done freeze DNA if not using that day

### Step 3 sgRNA synthesis, purification, and quantification

**sgRNA synthesis**
- Use [EnGen® sgRNA Synthesis Kit, _S. pyogenes_](https://www.neb.com/products/e3322-engen-sgrna-synthesis-kit-s-pyogenes#Product%20Information) for synthesis reaction and [Zymo RNA Clean & Concentrator](https://www.zymoresearch.com/collections/rna-clean-concentrator-kits-rcc/products/rna-clean-concentrator-25) for purification
- **IMPORTANT** When working with RNA (like in this protocol) wipe the surface and all materials with RNase Away before using. Use filtered tips if possible, and if non are available make sure to use freshly autoclaved tips. Always wear gloves and be extra careful about contamination
- Synthesis kit should be in the -20 door, and purification kit is in the cabinet with the extraction kits
- Set a heat block to 37 degrees C
- Thaw all kit components on ice, as well as the guide RNA primers
- Guide RNA primers need to be diluted to 1uM, if you have them resuspended at 100uM, they need to be diluted 1:100 in molecular grade water to a concentration of 1uM before use
- Flick to mix kit tubes and spin them all down, don't vortex any tube
- Primers can be vortexed and spun down
- Assemble **separate** reactions in 1.5mL tubes **on the bench,** not on ice
- You will need 1 synthesis reaction per sgRNA (you may only have 1 to make)
- Each reaction contains (in order):
  - 2ul molecular grade water
  - 10ul 2X sgRNA reaction mix
  - 5ul of 1uM primer
  - 1ul 0.1M DTT
  - 2ul sgRNA enzyme mix
- Pipette mix tubes with 10ul
- Spin down tubes
- Place tubes in the heat block at 37 for 30 minutes
- Afterwards, transfer the tubes to ice
- Add 30ul molecular grade water to each tube
- Add 2ul of DNaseI and pipette mix
- Spin down tubes
- Incubate tubes in the heat block at 37 for 15 minutes

**RNA clean and concentrator kit purification**
- Take tubes out of the heat block and let come to room temp
- Add 2X the volume of sample (104ul) of RNA binding buffer to each tube
- Pipette mix
- Add 1.5X the volume (231ul) 100% ethanol to each tube (the kit says to do 1.5X for RNA expected to be between 17-200nt in size, the NEB kit says they should be ~100nt)
- Pipette mix
- Transfer liquid from each tube into their own spin column
  - Spin column set up is one column slid into a collection tube
- Centrifuge columns at 15,800rcf for 30 seconds
- Discard flow through
- Add 400ul RNA prep buffer to each column
- Centrifuge columns at 15,800rcf for 30 seconds
- Discard flow through
- Add 700ul RNA wash buffer to each column
- Centrifuge columns at 15,800rcf for 30 seconds
- Discard flow through
- Add 400ul RNA wash buffer to each column
- Centrifuge columns at 15,800rcf for 1 minute
- Discard flow through
- Centrifuge columns "dry" at 15,800rcf for 30 seconds
- Transfer columns to new labeled 1.5mL tubes
- Added 40ul nuclease free water to each column and let sit for ~2 minutes
- Centrifuge columns at 15,800rcf for 30 seconds
- Placed tubes on ice

**Quantify**
- Follow the exact protocol as above for the dsDNA quantification but instead use the broad range RNA quantification kit, reagents, standards, and settings on the Qubit

### Step 4 Digestion Test

Following [this protocol](https://www.neb.com/protocols/2018/01/30/in-vitro-digestion-of-dna-with-engen-cas9-nls-s-pyogenes-m0646) from NEB

**Preparing sgRNAs**
- sgRNAs need to be diluted to a 300nM concentration for this protocol
- To determine nM concentration of your sgRNAs, it is a little complicated
- First, calculate the total amount of RNA synthesized:
  - ng/ul of sgRNA * ul = total ng
  - Example: 45ng/ul * 39 = 1,800ng or 1.8ug RNA
- Based off the NEB synthesis kit information, the length should be ~100nt for each sgRNA
- Use the NEB [mass to mole converter](https://nebiocalculator.neb.com/#!/ssrnaamt) to find the pmol of your sgRNA(s)
  - Example: 1.8ug at 100nt is 55.96pmol
- This can now be converted to nM in a few steps
- 1 pmol/ul = 1000nM
- Need the total volume of sgRNA(s) = should be 39ul
- Use conversion equation:
  - pmol / volume = pmol/ul * 1000 = nM
  - Example: 55.96/39 = 1.434 * 1000 = 1434nM
- Now you can dilute your sgRNA to 300nM with molecular grade water
  - You may want to determine how much volume sgRNA you will need to make for the dilution
  - For example, if you need 10ul of sgRNA at 300nM, you will need 3000nmols total
      - 3000nmols/1434nM = 2.09ul of sgRNA
      - 10ul total volume - 2.09 = 7.91ul of molecular grade water

**Determining nM concentration of PCR product**
- The PCR product needs to be added to the reaction at 30nM (although I have had success at 12nM)
- You will need to know the concentration of the DNA and the size of the PCR product
- Equation for DNA ng/ul to nM:
    - nM = ((ng/ul)/(660g/mol * size in bp)) * 1000000
- Example
  - PCR product is 34ng/ul and 500bp in size
  - ((34ng/ul)/(660g/mol * 500) * 1000000) = 103nM
- Dilute your PCR product with molecular grade water, you will need 3ul per reaction that uses that PCR product (number of reactions may vary, see below)

**Reaction setup**
- How many reactions will you do? You will need to do some controls as well as the experimental reaction. This is how many you would need to do if you had 1 PCR product and 1 sgRNA to test:

|tube number | PCR product|sgRNA|Cas9|reaction type|
|---|---|---|---|---|
|1|product 1|+ sgRNA 1|No Cas9| control|
|2|product 1|No sgNRA|+ Cas9|control|
|3|product 1|No sgRNA|No Cas9|control|
|4|product 1|+ sgRNA 1|+ Cas9|digest test|

- Dilute Cas9 to 1uM, it is at 20uM in the stock tube. You will need 1ul per reaction with Cas9, however if you need less than 20ul, it is easiest to make a 20ul diltuion:
  - 19ul diluent B
  - 1ul 20mM Cas9
  - This tube should be flick mixed, spun down, and kept on ice
- Set the thermocycler program:
  - 25 decree C hold
  - 25 degree C 10 minutes
  - 37 degree C hold
  - 37 degree C 15 minutes
- Program should have a 33ul volume
- Start the program so it is on the hold while you make your reactions
- Thaw all reagents and components on ice and make the mixes and reactions on ice
- _n_ is number of reactions * 1.1
- Make master mix:
  - 20ul molec grade water * _n_ = 140ul
  - 3ul NEB buffer r3.1 * _n_ = 21ul
- Pipette mix and spun down
- Add 23ul master mix to each tube
- **Note** here, the additions can get confusing, especially if you are working with multiple PCR products and/or sgRNAs. It can be helpful to number your tubes and when writing the protocol in your notebook to specify which number tube gets what
- Add 3ul 300nM sgRNA to tubes + for sgRNA
- Add 3ul water to tubes with no sgRNA
- Add 1ul 1nM Cas9 to tubes + for Cas9
- Add 1ul water to tubes with no Cas9
- Pipette mix tubes and spin them down
- Incubate the tubes in the thermocycler for 10 minutes at 25 degrees C (place tubes in and press skip step in the program. Set a timer for 10 minutes)
- Take out the tubes and add 3ul of PCR product to each tube that gets that product (do at room temp)
- Pipette mix and spin down tubes
- Incubate tubes in the thermocycler for 15 minutes at 37 degrees C (place tubes in and press skip step in the program, set a timer for 15 minutes)
- Take tubes out and add 1ul of Qiagen proteinase K to each
- Pipette mix and spin down tubes
- Incubate the tubes at room temp for 10 minutes

**Gel**
- Make at least a 2% gel to run the digestion to see separation of the bands. Also make it very thick with large combs if possible. It will run best if you use **15ul** of the digestion test liquid as input into the gel
- Use 3ul of 6X loading dye with 15ul of liquid, pipette mix, and add to the gel
