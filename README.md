# Tyrosine recombination classification for RIT element data-mining

Master's thesis work on tyrosine recombinases and recombinase-in-trio elements.

------------------------------------------------------------------------

## Objective 1

**Create a tool to classify 'phage integrase' proteins by their catalytic domains, in line with the 20 subfamilies described by Smyshlaev *et al.* (2021).**

------------------------------------------------------------------------

#### Datasets

- **SMART**: domain sequences and protein sequences for 20 subfamilies. Total number of sequences is ~120k, but with large class imbalance, where Xer has ~34k members and Int_Des has only 62. Notably some sequences are annotated with >1 subfamily, for the same portion of sequence (ambiguous classifications); any sequence belonging to more than one subfamily was removed. After filtering, 114,848 sequences remain.

- **non-integrase**: various searches of Refseq, some Pfam families of transposases, Uniprot proteomes of yeast, arabidopsis, and human proteomes. Large groups downsampled to 2000 sequences. A total of 39,691 sequences were retained after downsampling.

------------------------------------------------------------------------

#### Workflow

- [x]  Obtain SMART reference data & other non-integrase sequences as negative examples.
- [x] Create hold-out set for final testing
- [x] Align SMART domain sequences (for assessment and final models)
- [x] Build HMMs from training data
- [x] Score 20 integrase HMM alignments for each sequence    
- [ ] Gather HMM scores (& generate k-mer features?)  
- [ ] Tune & assess classifiers. Attempt stacking classifiers...

*Then...*  Apply classifier to MGE proteins and proteomes from    
assembled genomes (objective 2).
  
------------------------------------------------------------------------
  

### Scripts

All scripts are in the project directory `./code/`

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
  
  <!-- TODO Continue code documentation here. -->
  
<!-- - `./code/tidy_hmmsearch_res.R`: cleans up classifier data from hmmsearch results for both integrases and non-integrases. Saves files to *./data/smart_refseqs_hmm_scores.rds* and *./data/non_integrases_hmm_scores.rds*.   -->
<!-- -  `./code/add_kmer_features.R`: joins integrase and non-integrase datasets and computes kmer proportions, saves *./data/full_classifier_data.rds* for modelling. -->

------------------------------------------------------------------------

**4. Check accuracy of classification by hmmsearch best score**

- `4_classifier.R` (to do)
  - [ ] Check accuracy of classification by hmmsearch best score.
  - [ ] Nested cross-validation for assessment
  - [ ] Model tuning CV
  - [ ] Model stacking CV
  - [ ] Final model selection and training
  
  

<!-- - `./code/classifier1.R`: adds kmer profiles and splits data, does resampling for tuning and assessment. Trains final model and saves it.... -->

------------------------------------------------------------------------

 



## Objective 2

**Identify integrases in assembled genomes** and mobile element databases and **classify by subfamily**. For **RitA, RitB, and RitC** domains, map their genomic coordinates and check whether these constitute a **RIT element**.

#### Search MGE databases

#### Search assembled genomes

#### Search for flanking repeats

- Extract flanking genomic sequence
- Algorithm for identifying inverted repeat:
    - subsequence's reverse complement is the same as the reverse of the subsequence, for perfect base-pairing



---


### Other work done so far

Various scripts for obtaining & tidying mobile genetic element data to search: Iceberg, Aclame, pVOG, 

<!-- Other possible sources for MGE sequences.... PHAST (phaster), ISfinder, (others from Smyshlaev)? -->



\ 
    
    
    
----


New questions:

- what to do about ambiguous characters in sequences?? one sequence contains an XXXXX segment? **need to remove these seqs from kmer counting**? kmer::kcount doesn't work with them ...





