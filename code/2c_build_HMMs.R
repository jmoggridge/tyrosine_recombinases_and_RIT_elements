### Make training and final HMMs

# install HMMER with conda: conda install -c bioconda hmmer
# bash: cd into project home directory


## training domain HMMs ----

system('mkdir ./data/SMART/domain_hmm_training')

## loop to build training HMM for each subfamily
## system call isn't working; paste into terminal:
"
for file in ./data/SMART/domain_align_training/*.aln
  do
  outfile=`echo $file | sed 's_align_hmm_g;s_\\.aln_.hmm_g'`
  echo $file
  echo $outfile
  hmmbuild $outfile $file
done
"


## final domain HMMs ----

# loop to build final HMM for each subfamily
"
for file in ./domain_alignments/*.aln
  do
  outfile=`echo $file | sed 's/alignments/hmm/g;s/\\.aln/.hmm/g'`
  echo $outfile
  hmmbuild $outfile $file
done
"


# ### HMM domain searches on full protein sequences of 20 SMART groups
# 
system('mkdir hmmsearch_res')
# 
# for hmm in ./data/SMART/domain_hmm_train/*;
# do
# outfile=`echo $hmm | sed 's|domain_hmm/||g;s|.hmm|.search|g'`
# echo $outfile
# hmmsearch -o hmmsearch_res/$outfile.raw --tblout hmmsearch_res/$outfile.tbl --domtblout hmmsearch_res/$outfile.domtbl $hmm all_SMART_integrases.fa
# done
