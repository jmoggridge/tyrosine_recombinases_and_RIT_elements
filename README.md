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


**2**. Apply classifier to sequenced genomes and mobile element databases. For **RitA, RitB, and RitC** domains, map their genomic coordinates and check whether there are multiple tyrosine recombinases constituting a **RIT element**. Characterize the arrangement of subfamilies, flanking inverted repeats, and their insertion sites in terms of flanking genes.


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

