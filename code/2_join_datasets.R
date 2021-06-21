## join integrase and non-integrase data




## then create training and testing fasta files



# create fasta files for scoring train and test_data
dest <- './data/train_seqs.fa'
fasta <- Biostrings::AAStringSet(seq)
names(fasta) <- acc
writeXStringSet(fasta, filepath = dest)




