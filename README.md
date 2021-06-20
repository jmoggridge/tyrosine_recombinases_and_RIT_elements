# tyrosine_recombinases_and_RIT_elements
Master's thesis work on tyrosine recombinases and recombinase-in-trio elements



## Objective 1

**Create a tool to classify 'phage integrase' proteins by their catalytic domains, in line with the 20 subfamilies described by Smyshlaev *et al.* (2021).**

---

#### Datasets

- **SMART**: domain sequences and protein sequences for 20 subfamilies  

- **non-integrase**: various searches of Refseq, some Pfam families of transposases, Uniprot proteomes of yeast, arabidopsis, and human proteomes

---

#### Workflow

1.  Obtain SMART reference data & other non-integrase sequences as negative examples    
2.  Align SMART domain sequences
3.  Create hold-out set for final testing
4.  Build HMMs from training data

4.  Score 20 integrase HMM alignments for each sequence    

5.  Generate k-mer features    

6.  Tune & assess classifiers    

7.  Attempt stacking classifiers    

8.  *next* Apply classifier to MGE proteins and proteomes from    
assembled genomes.
  
  ---
  

#### Scripts

- `./code/1_tidy_smart_data.R`: reads domain and protein fasta sequences for 20 subfamilies from *./data/SMART/domain_fasta/* and */full_protein_fasta/*.  Joins domain and protein datasets into *./data/smart_df.rds*. Splits refences integrases into test and training datasets. Creates a set of fasta files of training domains for alignment
in *./data/SMART/training_domain_fasta*

- `./code/1b_tidy_non_integrases.R`: combines all non-integrase sequences in *./data/non_integrase_seqs/* and saves a fasta file for hmmsearches *./data/non_integrases.fa* and dataframe for the classifier *./data/non_integrases_df.rds*.  

- `./code/2_align_smart_domains.R`: creates 20 alignments for each of the subfamilies using 'decipher'. Saves them to *./data/SMART/domain_alignments/*.  

- `./code/3_train_test_split.R`: joins domains to their parent protein and then samples rows to create a hold-out set for final classifier testing. Outputs aligned fasta for each family to *./data/training_aligns/*.

- `./code/4_hmms_build_and_score_SMART.txt`: bash loop to create 20 HMMs with `hmmbuild`, saved to *./data/SMART/domain_hmm/* another to get scores for all SMART full proteins against each model with `hmmsearch.` and saves to *./data/SMART/hmmsearch_res/*.  

- `./code/get_refseq_non_integrases.R`: downloads non-integrases from ncbi entrez. Saves raw data for next script *.data//non_integrase_seqs/refseq_non_integrases_raw.rds*.  

- `./code/score_non_integrases.txt`: bash loop to score all non-integrases against the 20 hmms. Saves results to *./data/non_integrase_seqs/hmmsearc_res/*.  

- `./code/tidy_hmmsearch_res.R`: cleans up classifier data from hmmsearch results for both integrases and non-integrases. Saves files to *./data/smart_refseqs_hmm_scores.rds* and *./data/non_integrases_hmm_scores.rds*.  
-  `./code/add_kmer_features.R`: joins integrase and non-integrase datasets and computes kmer proportions, saves *./data/full_classifier_data.rds* for modelling.

- `./code/classifier1.R`: adds kmer profiles and splits data, does resampling for tuning and assessment. Trains final model and saves it....

Markdown version:   
  -`./thesis_objective1.Rmd`

----


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



