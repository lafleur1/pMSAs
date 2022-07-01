# paired multiple sequence alignments (pMSAs)

Summary: 

Attempting simple paired MSA generation for eukaryotes using OMA DB to feed into RoseTTAFold for complex prediction 

(Based off of the RoseTTAfold pMSA process for yeast, modified to use OMA DB, and trying to make a isoform-specific pMSA generation process ) 

1st, prepare each member of the heterocomplex: 


Run for each complex member 

1) Download both orthologs and paralogs from OMA for the target (TODO: update this to pull from OMA with their python API)
  a) Keep difference of set of ortholog and paralog sequences 
  b) IF DOING A NON-CANNONICAL ISOFORM RUN THOSE STEPS NOW
3) Select a single ortholog to use per species
  a) Remove multi-orthologs due to genome duplication issues 
  b) Remove multi-orthologs due to overlapping proteins (TODO: not a problem for my test case septin & S100 complexes) 
  c) Remove multi-orthologs due to ancient gene duplicates (using decreasingly stringent min delta  Seq90 cutoffs until at least one ortholog is selected per species) 
  
4) Create MSA: 
  a) Transform orthologs into a fasta DB w/ query as first sequence
  b) Run phmmer to make the seed alignment
  c) Run hmmbuild
  d) Run hmmlaign 
  e) Convert sto -> a3m 

 For non-cannonical isoforms: (IN PROGRESS) 

 Route 1: (Isoform re-BLAST and reselection)
 
 1) For every ortholog to the cannonical sequence, download any isoforms present 
  a) For every ortholog with isoforms, pBLAST the isoform against these and see if there's a better match than the cannonical sequnece 
  b) Select the best match to the isoform (TODO: find a better criteria for a match.....) and select that per species 
  
 Route 2: (Complete proteome re-BLAST for each species which had an ortholog to the cannonical sequence) 
 
 1) Get proteomes for every species with an ortholog to the cannonical isoform:
  a) NCBI proteome downloads with entrex (TOO SLOW for more of species in test case)
  b) OMA DB proteome downloads (with API) (TOO SLOW for amount of species in test case)
  c) Download and parse OMA DB proteomes for all species (TODO)


 MSA pairing:

 1) For each species: 
  a) Concatenate the wt sequences as the first line
  b) For all species following this: 
    i) If each of the partners in the complex have orthologs in the sequence, concatenate the alignments 
    ii) If one of them is in the species, concatenate it with gaps the same length as the other sequence 
    
 The pMSA is now ready to use for RoseTTAFold complex prediction (?) (Note that if using Robetta, you need to upload the MSA first, delete the autofill of the WT-WT sequence, and then paste in the WT-WT sequence with the chain break symbol) 
 
 TODO: figure out how to do the structural inputs.....
 
 #Septin complexes 
 
 Computing pMSAs for septin-12 interactions with septin 1, septin1 splice variant, and septin 5. 
 

 
| % Seq ID| S1 | S1DS | S5 | S12
| ------------- | ------------- | -------- |-------- |-------- |
|S1 | - | 100 | 63 | 45 |
| S1DS| 100 | - | 63 | 45|
| S5| 64 | 64 | - | 46|
| S12| 45 | 45 | 46 | -|
  

Filtering out paralogs from orthologs from OMA db for all proteins: 

| Protein       |# orthologs | # paralogs| # filtered orthologs |
| ------------- | ------------- | -------- |-------- |
| S1   | 1003   | 2123 | 1003 | 
| S5   | 1143  | 2060  | 1066  | 
| S12  | 1034  | 2126  | 995 | 

