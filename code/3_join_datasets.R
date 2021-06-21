## join integrase and non-integrase data


# smart data
smart_train <- read_rds('./data/SMART/smart_train.rds')
smart_test <- read_rds('./data/SMART/smart_test.rds')



try <- Biostrings::readAAStringSet('./data/SMART/smart_train.fa')
try

# non_integrases


## create fasta file for HMM scoring
# 
# # create fasta files for scoring train and test_data
# dest <- './data/train_seqs.fa'
# fasta <- Biostrings::AAStringSet(seq)
# names(fasta) <- acc
# writeXStringSet(fasta, filepath = dest)

## score HMMs





