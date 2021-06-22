## ./code/2c_build_HMMs.sh

## Make training and final HMMs; score train and test sequences against 'training' hmms 
## created from training sequences only.

# First: install HMMER with conda: conda install -c bioconda hmmer
# in bash: cd into project home directory
# takes ~1hr

mkdir ./data/SMART/domain_hmm_training
mkdir ./data/SMART/domain_hmm_final
mkdir ./data/hmmsearch_res/

# loop to build training HMM for each subfamily
for file in ./data/SMART/domain_align_training/*.aln
  do
  outfile=`echo $file | sed 's_align_hmm_g;s_\\.aln_.hmm_g'`
  echo $file
  echo $outfile
  hmmbuild $outfile $file
done

# loop to build final HMM for each subfamily
for file in ./data/SMART/domain_alignments/*.aln
  do
  outfile=`echo $file | sed 's/alignments/hmm_final/g;s_\\.aln_.hmm_g'`
  echo $outfile
  hmmbuild $outfile $file
done

# hmmsearch using the training models
for hmm in ./data/SMART/domain_hmm_training/*;
  do
  outfile=`echo $hmm | sed 's|SMART/domain_hmm_training|hmmsearch_res|g;s|\\.train\\.hmm||g'`
  echo $hmm
  hmmsearch --noali -o temp --tblout $outfile.train.tbl $hmm ./data/train_seq.fa
  hmmsearch --noali -o temp --tblout $outfile.test.tbl $hmm ./data/test_seq.fa
done

rm temp