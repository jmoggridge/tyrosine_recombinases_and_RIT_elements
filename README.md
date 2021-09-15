
## Tyrosine recombination classification and mining for Recombinase-in-Trio elements  

**BINF6999 Thesis Project**    

*Documentation in progress!*

--------------------------------------------------------------------

#### Contents
1. [Project Objectives](#obj)
2. [Part 1: Classifier for tyrosine recombinases](#p1)
   - [Overview](#p1o)
   - [Workflow](#p1w)
3. [Part 2: Pipeline for identification of RIT elements](#p2)
   - [Overview](#p2o)
   - [Workflow](#p2w)


----------------------------------------------------------------------

### Project Objectives <a name="obj"></a>
    
    
**1 **. Create an amino acid sequenece classifier to assign site-specific tyrosine recombinases (aka. phage integrase) proteins to the 20 subfamilies described by [Smyshlaev *et al.* (2021)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138268/), or 'Other' if the protein is not a tyrosine recombinase.


**2**. Create a pipeline for the identification of Recombinase-in-Trio (RIT) mobile elements from GenBank files, using the tyrosine recombinase classifier to identify RitA, RitB, and RitC proteins.


**3**. Classify tyrosine recombinases in the [ICEberg2.0 database](https://bioinfo-mml.sjtu.edu.cn/ICEberg2/index.php).
to sequenced genomes and mobile element databases. 


------------------------------------------------------------------------
 

### Part 1: Classifier for tyrosine recombinases (YRs)<a name="p1"></a>


#### Overview  <a name="p1o"></a>
  
- **Datasets**
  - SMART tyrosine recombinases reference data & other non-integrase sequences as negative examples.
  - These are combined, downsampled, and split to create a training/validation dataset for classifier development and a holdout dataset for final testing.
  
  
- **Modelling Pipeline**
   - For a given data split (eg. training/validation samples):
   - Create 20 alignments of YR catalytic domains, one for each subfamily.
   - Create 20 HMMs from the alignments, one for each subfamily
   - Score train and test seqs against HMMs using HMMER's `hmmsearch`
   - Normalize scores (using training data only to create normalizer)
   - Use SMOTE upsampling to create synthetic examples, bringing minority classes up to 25% of the largest class.
   - Fit classifiers to the prepared training sample.
   - Apply classifiers to validation sample. 
   - Evaluate predictions and return performance metrics.
     
- **3-fold Nested-CV (repeated 3 times)**
   - Apply the above pipeline, with a variety of models & hyperparameter settings, in nested CV.    
  -  Get an estimate of performance for the model-fitting procedure (not actually selecting the model hyperparameters at this stage).    
  
- **3-fold CV (rep 3 times)**    
  - For the selection of models with best mean MCC across.
  - Best models are kept for the final holdout test.
  
- **Holdout test**
  - Fit 'best' models from CV using pipeline, with all training data.
  - Get point estimates of performance from prediction of holdout data.

- **Fit final models** 
  - Align all catalytic domains for each subfamily, build 20 HMMs, score all sequences by `hmmsearch`, normalize scores, upsampling by SMOTE, fit classifiers.

- **Produce visualizations** 
  - To show evaluation of classifiers in nested CV, CV, and final test.
 

------------------------------------------------------------------------

#### Workflow  <a name="p1w"></a>


##### 1. Dataset acquisition:

**SMART database**

  - The script `1a_tidy_smart_data.R`:
  - Organizes data from the set of fasta files downloaded from the [SMART database](http://smart.embl-heidelberg.de/).  
  - Creates a dataframe with the names and sequences for tyrosine recombinase proteins and their catalytic domains. 
  - The total number of sequences is \~120k, but with large class imbalance, *e.g.*, Xer has \~34k members and Int_Des has only 62.    
  - Problematically, some sequences are annotated with \>1 subfamily, for the same portion of sequence (ambiguous classifications). These double-labelled sequences are removed, leaving 114,848 sequences.     
  - The dataset for modelling was downsampled to a max of 10k sequences per subfamily, with the removed sequences set aside for further classifier testing.   
  - The cleaned data was saved as `smart_df.rds`
  
**Non-integrases**

  - The script `1b_get_refseq_non_integrases.R` sent keyword searches to Entrez API to get a variety of non-YR protein sequences to serve as negative examples. 
  - I also downloaded:
      -   Several Pfam families of transposases (domain sequences only)    
      -   Uniprot proteomes of yeast, Arabidopsis, and human proteomes.    
  -   Any large groups were down-sampled to 2000 sequences.    
  -   A total of 39,691 sequences were retained after downsampling.    

  - The script `1c_tidy_non_integrase.R` cleans up data gathered by the steps above to create the dataset `nonint_df.rds`


#### 2. Searching for best models

The script `2a_data_splitting.R`:
  - Combines the SMART integrases dataset and the non-integrase dataset.
  - Splits data 75:25 for training/validation (building & tuning classifiers) and testing
  - Creates 2 files in `./data/`:
    - `classif_train_set.rds` (training 75 %)
    - `classif_test_set.rds` (testing 25 %)

The script `2b_nested_CV.R`:
  - Performs 3 x 3-fold nested CV 
  - Uses functions from the script `00_functions.R` for the model fitting steps
  - Uses the set of model specifications created in the script `00_get_model_specs.R` (models: `unfitted_parsnip_model_set.rds`)
  - Results are saved to `results/3x3-fold_07-08_nest_cv_summary.rds`

The script `2c_regular_CV.R`


#### 3. Final validation and model fitting

- `3a_final_test_set.R` - fits best models from CV using the full training data, and evaluates predictions of the test set
- `3b_fit_final_model.R` - fits the best model of each type (Elastic Net, Random Forest, k-Nearest Neighbors) using the full, down-sampled dataset


#### Visualizations

`4_Classifier_results_plots.R` generates the figures for the manuscript and presentation

------------------------------------------------------------------------


### Part 2: A pipeline to identify new RIT elements  <a name="p2"></a>

#### Overview  <a name="p2o"></a>

- Obtain data from NCBI by linking CDD ->> proteins ->> nucleotides -> taxonomy
- Translate nucleotide sequences -> proteins with locations
- Classify protein sequences
- Figure out which nucleotides have RIT arrangement


#### Workflow  <a name="p2w"></a>



-----------------------------------------------------------------------------



<!--
### Other random scraps of work done so far

Various scripts for obtaining & tidying mobile genetic element data to search: Iceberg, Aclame, pVOG,

 Other possible sources for MGE sequences.... PHAST (phaster), ISfinder, (others from Smyshlaev)? -->