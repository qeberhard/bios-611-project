#First we load in the unspliced and spliced data
unspliced_txptome <- read_csv(derived_data/unspliced_txptome.csv)
spliced_txptome <- read_csv(derived_data/spliced_txptome.csv)

#Next, let's find the number of RNAs that are spliced/unspliced
unspliced_rna_count <- length(unspliced_txptome$RNA)
spliced_rna_count <- length(spliced_txptome$RNA)

#As you cacn see, there are more spliced version than their are unspliced.
