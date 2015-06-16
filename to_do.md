- [x] Rebundle the viral files for downloading
- [ ] Check to make sure the BMG files are correct
- [ ] start running code for BMG metagenomes as well
- [ ] take out cloning repository step? might need for figure examples? or put in interm directory?
- [ ] add trimmomatic and cd-hit-454 results to interm
- [ ] have the interm file as a compressed file for DL, change in the markdown
- [ ] add introduction about data - total vs. viral (VMG / BMG), research questions, hypothesis, etc
- [ ] add PBS cluster information to commands?
- [ ] add portion to explain why sample dep't trimming -- each library prep a bit diff with transposon
- [x] add email to markdown
- [ ] explain why the QC for total is different than viral
- [ ] add viral to outputs, like viral_cd_hit_454 and change to code to reflect this
- [ ] set home (~) varaible at beginning to make it easier to get around
- [ ] after QC, change the name of files so short, qc.fastq or something
- [ ] once done with QC, check all steps by running in Tusker
    - [ ] have shell script to run certain portions? like shell script to run all viral QC steps
- [ ] put warning about generating lot of data, make sure have room for files
- [ ] make sure interm file names match waht they should actually be - cdhit454 is wrong now
- [ ] look at and list modules needed for perl (python too?) for prinseq
- [ ] change all scripts to remove unnecessary cd .. steps and just call direct path of script in home
- [ ] make total trimmomatic step (first) a loop since not sample dep't
- [ ] re-compress the total metagenome samples. rename?
- [ ] if you make trimmomatic be underscore trimm then naming will work