# Tyrosine recombination classification for RIT element data-mining

Master's thesis work on tyrosine recombinases and recombinase-in-trio elements.

------------------------------------------------------------------------

## Objective 1

**Create a tool to classify 'phage integrase' proteins by their catalytic domains, in line with the 20 subfamilies described by Smyshlaev *et al.* (2021).**

------------------------------------------------------------------------

#### Datasets

-   **SMART**:

  -   Has domain sequences and protein sequences for 20 integrase subfamilies.

  -   Total number of sequences is \~120k, but with large class imbalance, like Xer has \~34k members and Int_Des has only 62.

  -   Problematically, some sequences are annotated with \>1 subfamily, for the same portion of sequence (ambiguous classifications). I removed these double-labelled sequences, entirely.

  -   After filtering, 114,848 sequences remain.  

  -   The dataset for modelling was downsampled by subfamily to a max of 10k, with the removed sequences retained for classifier validation.

*show bar chart here*

-   **Non-integrase**:

  -   Combination of:

      -   From various searches of Refseq for recombinases and nucleases,
      -   Several Pfam families of transposases (domain sequences only)
      -   Uniprot proteomes of yeast, Arabidopsis, and human proteomes.

  -   Large groups down-sampled to 2000 sequences.

  -   A total of 39,691 sequences were retained after downsampling.

------------------------------------------------------------------------

#### Workflow

-   [x] Obtain SMART reference data & other non-integrase sequences as negative examples.
-   [x] Create hold-out set for final testing.
-   [x] Align SMART domain sequences (for assessment and final models).
-   [x] Build HMMs from training data.
-   [x] Score 20 integrase HMM alignments for each sequence.
-   [x] Gather HMM scores (& generate k-mer features).
-   [ ] Get baseline accuracy from HMM scoring approach.
-   [ ] Tune & assess classifiers. Possibly: attempt stacking classifiers.

*Then...* Apply classifier to MGE proteins and proteomes from  
assembled genomes (objective 2).

------------------------------------------------------------------------

### Scripts

All scripts are in the project directory `./code/`

------------------------------------------------------------------------

## Objective 2

**Identify integrases in assembled genomes** and mobile element databases and **classify by subfamily**. For **RitA, RitB, and RitC** domains, map their genomic coordinates and check whether these constitute a **RIT element**.

#### Search MGE databases

#### Search assembled genomes

#### Search for flanking repeats

-   Extract flanking genomic sequence

-   Algorithm for identifying inverted repeat:

    -   subsequence's reverse complement is the same as the reverse of the subsequence, for perfect base-pairing

------------------------------------------------------------------------

### Other work done so far

Various scripts for obtaining & tidying mobile genetic element data to search: Iceberg, Aclame, pVOG,

<!-- Other possible sources for MGE sequences.... PHAST (phaster), ISfinder, (others from Smyshlaev)? -->




------------------------------------------------------------------------

## New workflow

-   combine data and split data for CV
-   setup pipeline to align, create hmms for training data
-   score train and test sequences for each fold
-   do training on scored train sequences & evaluation on scored test sequences

\*\*

#### June 9th

Restarted nested CV - runs ~40 hrs.
 - Created a new set of model specifications, including larger search grid for glmnet and rpart models, added random forest models with 5 mtry parameters.
 - Implemented new recipe with SMOTE and normalization (within CV)
 - Executed nested 3-fold CV, 3-repeats. Results look awesome again (too awesome).
 - Had to fix many bugs in functions arising from adding SMOTE and different recipes...

#### June 10th

Created script that gets all ids for CDD families Rit- A, B, C. Gets cdd ids, gets protein ids linked to each cdd id. Then gets nuccore ids linked to each protein id. Then gets taxonomy id linked to each nuccore id.

Downloading the sequences is more tricky. Keep getting HTTP errors. Probably sending too many requests - \> maybe use post request then set get request using the webhistory token. Alternately add sleep between blocks of x ids sent in batches.

#### June 11th

Nested CV finished - -saved ./results/07-08

Still trying to rewrite script to retrieve NCBI data linked to cdd specific proteins...

Restarting regular CV with new model set. 
Needed to change:
 - models
 - {out_path} issues: needs to be complete path to output dirs for aligns, hmms....
 


<!--
**1. Data acquisition, tidying, joining**

- [x] `1a_tidy_smart_data.R`
  - Reads domain and protein fasta sequences for 20 subfamilies from *./data/SMART/domain_fasta/* and *full_protein_fasta/*.  
  - Joins domain and protein datasets into *./data/SMART/smart_df.rds*. 
  - Removes sequences found in more than 1 subfamily.
  - Splits ref integrases into test and training datasets. 

- [x] `1b_get_refseq_non_integrases.R`
  - Downloads various groups of non-integrases from NCBI entrez api.
  - Saves to *./data/non_integrase_seqs/refseq_non_integrases_raw.rds*.  

- [x] `1c_tidy_non_integrases.R`
  - Combines all non-integrase sequences in *./data/non_integrase_seqs/*.
  - Tidies data up to match integrase dataset.
  - Data sanity checks & filtering.
  - Train/test split 
  - Saves *nonint_train_df.rds* and *nonint_test_df.rds* dfs to *./data/non_integrase_seqs/*.

------------------------------------------------------------------------

**2. Alignment, HMM building**

- [ ] `2a_align_training_domains.R` (rerun)
  - Alignment of training domains for each of the subfamilies. 
  - Uses training domain sequences from *./data/SMART/smart_train.rds* created in `1a`
  - Saves them to _./data/SMART/domain_align_training/*.train.aln_.

- [x] `2b_align_all_domains.R`
  - Same as 2a but aligns all domain sequences to create the final HMMs.
  - Uses all domain sequences from *./data/SMART/smart_df.rds* created in `1a`
  - Saves them to *./data/SMART/domain_alignments/*.


------------------------------------------------------------------------

**3. Consolidate data, score sequences, join scores, prepare for classifier**
- [x] `3a_join_data.R`
  - Splits test/train from non_integrase data.
  - Consolidates integrases (SMART) & non-integrases data into training and test dataframes for the classifier. These are saved in ./data/ as *train_df.rds* and *test_df.rds*.
  - Consolidates fasta files for hmmsearch scores: *train_seq.fa* & *test_seq.fa*

- [x] `3b_hmmbuild_and_hmmsearch.sh`
  - Bash script to run from the project directory.
  - Builds the training and final HMMs from the alignments in step 2.
  - Saves training HMMs to _./data/SMART/domain_hmm_training/_, and the final HMMs to _./data/SMART/domain_hmm/_.
  - Run hmmsearch for sequences against 20 HMMs made from training sequences.
  - Save hmmsearch tables to _./data/hmmsearch_res/_

- [x] `3c_prep_for_classifier.R` (needs refactoring)
  - Processes hmmsearch output and joins to train and test data
  - Add kmer counts for each sequence
  - Add longer Dayhoff alphabet kmers? (5mers = 7,700 cols)
  
  
------------------------------------------------------------------------

**4. Classification: model selection and assessment**
- `4a_eda.R` (to do)
  - [x] plot sequence lengths
  - [ ] plot sequence composition profiles
  - [x] plot hmmsearch scores
  
- `4_classifier.R` (to do)
  - [ ] Check accuracy of classification by hmmsearch best score.
  - [ ] Nested cross-validation for assessment
  - [ ] Model tuning CV
  - [ ] Model stacking CV
  - [ ] Final model selection and training
 
  -->

<!-- TODO Continue code documentation here. -->

<!-- - `./code/classifier1.R`: adds kmer profiles and splits data, does resampling for tuning and assessment. Trains final model and saves it.... -->
