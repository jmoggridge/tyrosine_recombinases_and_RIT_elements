## Tyrosine recombination classification for RIT element data-mining

*Master's thesis work on tyrosine recombinases and recombinase-in-trio elements.*

--------------------------------------------------------------------

#### Contents
1. [Objectives](#obj)
2. [Part 1: Classifier for tyrosine recombinases](#p1)
   - [Overview](#p1o)
   - [Workflow](#p1w)
3. [Part 2: Pipeline for identification of RIT elements](#p2)
   - [Overview](#p2o)
   - [Workflow](#p2w)
4. [Daily Progress Tracking](#dpt)

----------------------------------------------------------------------

### Objectives <a name="obj"></a>

**1 **. Create a tool to **classify 'phage integrase' proteins** by their catalytic domains, in line with the **20 subfamilies** described by Smyshlaev *et al.* (2021).


**2**. Apply classifier to known integrases in sequenced genomes and mobile element databases. For **RitA, RitB, and RitC** domains, map their genomic coordinates and check whether there are multiple tyrosine recombinases constituting a **RIT element**. Characterize the arrangement of subfamilies, flanking inverted repeats, and their insertion sites in terms of flanking genes.


------------------------------------------------------------------------
 

### Part 1: Classifier for tyrosine recombinases <a name="p1"></a>


#### Overview  <a name="p1o"></a>
  
- [x] Obtain SMART reference data & other non-integrase sequences as negative examples.
- [x] Downsample & Split data: Create hold-out sets for final testing.
- [x] Modelling Pipeline - to take a data split, return model performance measures. 
   - create 20 alignments from YR catalytic domains of training seqs
   - create 20 HMMs from the alignments
   - score train and test seqs against HMMs
   - normalize scores and use SMOTE upsampling
   - train classifier with prepared training data
   - apply classifier to testing set. 
   - evaluate predictions and return performance metrics.
     
-   [x] 3-fold Nested-CV (rep 3 times)
    -  Apply above pipeline to with variety of models & hyperparameters for nested CV resamples.    
    -  Get an estimate of performance for *the model-fitting procedure* (not actually selecting the model hyperparameters for the final model).    
-   [x] 3-fold CV (rep 3 times).    
    -  Selection of models with best mean MCC across.
-   [x] Last validation    
    - Fit 'best' models from CV using pipeline, with all training data (75 %)  
    - Get point estimates of performance from prediction of holdout data (25%) 
-   [x] Fit final models. 
    - Align all seqs, build HMMs, score all seqs, normalize & SMOTE, train classifier.

*Then...* Apply classifier to MGE proteins and proteomes from assembled genomes (objective 2).


------------------------------------------------------------------------


#### Workflow  <a name="p1w"></a>

##### 1. Obtain data:

**SMART**

`1a_tidy_smart_data.R` organizes SMART data from fasta files.  
   \
   \
-   Has domain sequences and protein sequences for 20 integrase subfamilies.    
-   Total number of sequences is \~120k, but with large class imbalance, like Xer has \~34k members and Int_Des has only 62.    
-   Problematically, some sequences are annotated with \>1 subfamily, for the same portion of sequence (ambiguous classifications). I removed these double-labelled sequences, entirely.    
-   After filtering, 114,848 sequences remain.     
-   The dataset for modelling was downsampled by subfamily to a max of 10k, with the removed sequences retained for classifier validation.    
<!-- *show bar chart here* -->
  
**Non-integrases**


`1b_get_refseq_non_integrases.R`  
   - Uses keyword searches to Entrez API to get a variety of decoys.    
  

  -   Combination of:  
      -   From various searches of Refseq for recombinases and nucleases,    
      -   Several Pfam families of transposases (domain sequences only)    
      -   Uniprot proteomes of yeast, Arabidopsis, and human proteomes.    
  -   Large groups down-sampled to 2000 sequences.    
  -   A total of 39,691 sequences were retained after downsampling.    

  
`1c_tidy_non_integrase.R`  
   - Cleans up data gathered by step 1b


#### 2. Searching for best models

2a. Data splitting for classifier
  
  - 75-25 split for train (building & tuning classifiers) / test (final validation)
  - 2 files in `./data/`:
    - `classif_train_set.rds` (training 75 %)
    - `classif_test_set.rds` (testing 25 %)

2b. Nested cross-validation   

  - Nested cv: `2b_nested_CV.R`
  - Plots:     `2b_nested_cv_tables_and_plots.R`
    
2c. Regular CV

  - `2c_regular_CV.R`
  
See also `Classifier_results_plots.R` -- TODO combine all scripts with these various plots into one place for consistency....

#### 3. Final validation and model fitting

- `3a_final_test_set.R` - fit/eval of best models on the final validation set
- `3a_final_test_set.R` - visualize results from 3a.
- `3c_fit_final_model.R` - fits one best model of each type using the full, downsampled dataset





------------------------------------------------------------------------


### Part 2: A pipeline to identify new RIT elements  <a name="p2"></a>



#### Overview  <a name="p2o"></a>

- Obtain data from NCBI by linking CDD ->> proteins ->> nucleotides -> taxonomy
- Translate nucleotide sequences -> proteins with locations
- Classify protein sequences
- Figure out which nucleotides have RIT arrangement


#### Workflow  <a name="p2w"></a>



#### Other random scraps of work done so far

Various scripts for obtaining & tidying mobile genetic element data to search: Iceberg, Aclame, pVOG,

<!-- Other possible sources for MGE sequences.... PHAST (phaster), ISfinder, (others from Smyshlaev)? -->

-----------------------------------------------------------------------------


## Daily Progress Tracking <a name='dpt'></a>

#### July 8th

Restarted nested CV - runs ~40 hrs.  

- [x]  Created a new set of model specifications, including larger search grid for glmnet and rpart models, added random forest models with 5 mtry parameters.    
- [x]  Implemented new recipe with SMOTE and normalization (within CV).   
- [x]  Executed nested 3-fold CV, 3-repeats.   
- [x]  Had to fix many bugs in functions arising from adding SMOTE and different recipes...  

#### July 9th  

(*this wouldn't work on the full dataset, but only on the small testing set*)    
- Created script that gets all ids for CDD families Rit- A, B, C.     
- Gets cdd ids, gets protein ids linked to each cdd id.     
- Then gets nuccore ids linked to each protein id.     
- Then gets taxonomy id linked to each nuccore id.      

  
Downloading the sequences is more tricky. Keep getting HTTP errors. Probably sending too many requests - \> maybe use post request then set get request using the webhistory token. Alternately add sleep between blocks of x ids sent in batches.  
  
#### July 10th  
  
Nested CV finished today.  
 - saved `./results/07-08-3x3_nested....`.   
 - Results look awesome again (too awesome).  
  
Still trying to rewrite script to retrieve NCBI data linked to cdd specific proteins...  
 - trying to use webhistories  
 - trying to chunk requests into smaller than 50 for protein,  
 - trying to chunk into smaller than 10 for nucleotides... (not necessary, actually) 50 or 100 is ok.).
 

Restarting regular CV with new model set.   
Needed to change:  
 -  models  
 - {out_path} issues: needs to be complete path to output dirs for aligns, hmms....  
 - otherwise everything worked.  
 - Saved to `results/regular_cv_07-10/`   
 - *could* do stupid t-tests. to check for significance....  
 
#### July 11th

Running `3_final_test_set.R`  

   - had to change recipe to include smote and normalize   
   - changed file paths and names to make more sense    
   - Issue: system runs out of memory...   
      - I think this was before when trying to save unnested workflows/preds/data    
   - Decided: just use the few  specifications that performed best in regular CV    
   - Works fine for smaller amount (<10) models.    
   
Downloading `P2_A1_ncbi data.R`:      

- [x]  Finally managed to get all the sequences!  
- [x]  **TRICK** some genomes are really, large, need to get sequences in large chunks of 1000; do requests within in mini chunks of 100. Parse sequences and save compressed data. **Issue** even if managing to get full download complete before, ran out of memory trying to parse or save data.
 
- [ ]  still need the taxonomy but no rush  

    
#### July 12th (monday)
 
**Classifiers**

  - Re-ran final classifiers - one of each glmnet, knn, rf.
    - `3c_fit_final_model.R` script updated -
    - Now creates one model for each of rf, knn, glmnet; and two sets of thresholds: normalized+smoted or from raw hmm scores...  
    - Still missing the parts to generate the full alignments from all smart domains, and then to create the HMMs from these, but can copy from where those were actually done or rename the files...  
 
 
**Part 2 data**

 - Continuing `P2_A1_ncbi_downloads.R`
   - linking taxonomy ids to nucleotide ids. (`link_nuccore_taxonomy`)
   - downloading taxonomy summaries (old `fetch_taxonomy`)
   - weird error when trying to get nucleotide summaries: seems random
      ` x No esummary records found in file `
   - splitting into half and using `fetch_summary` to get nested data worked on the first half.
   - some cols with the same name have different types (list vs char or num); if we make all columns into character, then they can be unnested together. Not if the co
   - the dataframes have different numbers of columns.
   - **Solved** by making columns as.character in each dataframe, then unnesting... 
   - **Finished** Downloaded ids, sequences, summaries for cdd, prot, nuc, tax.   
   
 - Working on `P1_A2_repair_nucleotide_data.R`:
   - [x] fixed 2 superceded sequences and replaced ids, got new taxonomy and seq. *still need to replace seq*
   - [x] fixed 2 updated seqs... got data + linked tax
   - still need to remove 205 seqs 
   - get updated taxonomy for remaining that are missing it from the nuc_summary$taxid column
      
     
#### Tues, July 13th
 
**Part 2**

- Simplified `P2_A2_filter_nuc_ids_and_taxonomy` - first filtering is done now
- Wrote `P2_A3_fix_nucleotide_datasets.R` 
  - fixes the nucleotide datasets -> to *./data/CDD/nuc_data_fixed/*
  - tries to get the missing data for remaining few sequences
  - 9 ids are from assemblies - mulitple contigs
  - saved gi # (nuc_ids) from browser by searching nuc_id in nucleotide, then finding the assembly, going to that, and finding the 'nucleotide Refseq' or 'nucleotide INSD' links, then downloading the list to *./data/CDD/assemblies_gi_lists/*
  
- Updating clustering:
  - Added Dunn index and some other indexes using nb_clust
  - Not sure if wanting to use PCA figures, would need to try to explain, requires many dimensions but also have many groups.
  
    
#### Weds, July 14th
  
  - meeting: prioritize RIT finding + description
  - wrote genbank fetcher and parser functions: 
  

#### Thurs, July 15th

  - started: `P2_A4_get_genbank_data_for_probable_RITs.R`
    - filter ids down to nucleotide ids with three linked RIT proteins.
    - get genbank records for those - split the downloads into chunks of 100. Saved as `./data/CDD/RIT_gbk_<1-6>.rds`
  - still having some issues with iceberg data - HMM scores not being created in right folder...
  
#### Friday, July 16th

Iceberg data to Nicole - *done*.     
     
Working on P2_A5 - Rit identification pipeline. Documentation:     
     
**The main wrapper `rit_finder()` **:     

   - uses set of probable Rits (`ids_w_three_integrases.rds`)     
     - these already have genbank data from `P2_A3` + `A4`.     
   - parses genbank record with `open_genbank()`     
   - extracts feature table with ` extract_features_table() `     
   - classifies all translations of the CDSs with `classify_proteins()`     
     - uses consensus of knn and glmnet to create `consensus_pred` column.     
   - tidies up start, stop, orientation with `edit_feature_table()`.     
   - take lag + lead row data with `pivot_feature_table`     
     - now each CDS has linked up- and downstream CDSs for testing.     
   - Test trios with `rit_tester` to check whether trio of proteins are:     
     - `rit_length_check`: p1 start to p3 stop < 4000 bp     
     - `rit_distance_check`: CDS not separated by gap > 250 bp      
     - `rit_all`: all Rit subfamily members,     
     - `rit_ABC`: trio is Rit<A, B, C> or Rit<C, B, A>     
     - `rit_orientation`: all forward or reverse? or not all same?     
   - Select which are RITs with `rit_selector()` based on `rit_all`, `_distance_check`, and `_length_check`   
   - Return value is a tibble with nuc_id and nested tibble of Rit trios!    
   
**`rit_finder` works on nucs 1-4, not 5, 6-8...    **

  - **Issue** with HMM searches coming back empty for some domains.    
  - **Solution** create filler dataset `./data/hmmsearch_filler.rds` to concatenate with query sequences, such that there is always at least one hit and the hmmsearch output is not an empty table. Empty table causes an error because read_hmmsearches is expecting specific columns to be present in the table.    
  
  - **Issue** Some genbank files are full of partial sequences and pseudogenes but the only info I have is the absence of protein_id and translation values...    
  - **work-around** label these as pseudo_or_partial = TRUE.    
  
  - **Issue** not sure what data to collect for later analysis    
  - **action** show Nicole what I have so far...    
     
  - **Issue** no CDS....    
   - e.g.: 1817592545 (5),     
  - **Possible** download their files manually and parse the xml.    
      - **Issue** parse_genbank() 'Error: object 'GBSeq' not found'    
    
 
#### Saturday, July 17th   

Today:
1. Find out which Genbank files in three_ints don't have the necessary CDS.
2. Figure out if it is possible to create a different parser for genbank records
downloaded via browser for those currently missing CDS.
    
Tried to extract RIT sequences by locating protein with translation of ICEs for all reading frames, many protein simply not present in entire set of translations. Protein sequences do not seem to match the DNA sequences.

Instead, just blasted proteins against strain to find location.
Went to genbank and expanded range to find trio. Took start and stop from the range containing the trio. 
Aligned mesorhizobium rits and they were the same.
Only get one blast hit despite triplicate of double RIT segment.

Trying to finish RIT finder script. Still finding weird issues:
- "47118329" -> "Column `locus_tag` doesn't exist."
- "339284117", "737980678", '1980667557' ->  "Error in gzfile(file, "rb") : invalid 'description' argument"


#### Sunday, July 18th   

Trying to process all genbank files with `rit_finder` (P2_A5)  

- up to 475

Trying to figure out how to answer Nicole's qs about ICEs' RITs.

- are ICEs distinct 
  - 
Bordetella RITs:

RIT1 
  - rows (1,2,3) 
  - identical - except for insertion in RIT 3
  - one on each ICE...
RIT2
  - different from RIT1

Mesorhizobium
 - ICEberg data needs labels a,b,c for the three ICEs  
 - each ICE is not uniquely identified in the fasta, sequences are different
 - this is causing wrong proteins to be joined by id....
 - TODO need to fix



##### TO DO

- Find all RITs
- ICEberg -> locate, check if duplicated...
- Fix mesorhizobium in ICEberg data in P2_B3 before predictions...
- Start writing some parts...


### Writing

- METHODS
   - part1 data
   - part1 classifier overview
   - 

**Part 1**

  - [ ] redo plots for regular cv - compare to nested ; incorp final valid resuls and plot
  - [ ] redo k-means with letting k be variable...
 
**Part 2**

  - [ ] 



  
---------

<!-- USE THESE OLD SCRIPTS TO REUSE EXISTING CODE
- [x] `2b_align_all_domains.R`
  - Same as 2a but aligns all domain sequences to create the final HMMs.
  - Uses all domain sequences from *./data/SMART/smart_df.rds* created in `1a`
  - Saves them to *./data/SMART/domain_alignments/*.
  
- [x] `3b_hmmbuild_and_hmmsearch.sh`
  - Bash script to run from the project directory.
  - Builds the training and final HMMs from the alignments in step 2.
  - Saves training HMMs to _./data/SMART/domain_hmm_training/_, and the final HMMs to _./data/SMART/domain_hmm/_.
  - Run hmmsearch for sequences against 20 HMMs made from training sequences.
  - Save hmmsearch tables to _./data/hmmsearch_res/_

-->



