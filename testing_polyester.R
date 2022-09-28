#So brave, my first R script!
#Let's try out the Polyester package for simulating reads 

library(polyester)
library(Biostrings)

#generate a fold change matrix, all 1s since I don't want diff expression 
fold_changes=matrix(c(1), nrow=20, ncol=2)

#load fasta file
#fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
#fasta = readDNAStringSet(fasta_file)

fasta_file='C:/Users/chand/Box Sync/Krasileva_Lab/Research/chandler/Krasileva Lab/E14/NLR_transcripts.fa.txt'
fasta=readDNAStringSet(fasta_file)

small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'NLRs_small.fa')
readspertx = round(20 * width(small_fasta) / 100)

simulate_experiment('NLRs_small.fa', reads_per_transcript=readspertx, 
                    num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads') 
