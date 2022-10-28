library(ggplot2)
library(stringr)


#read in data
txptome <- read.csv("source_data/5mers_expressed_txpts.csv.gz")
colnames(txptome)[1] <- "RNA"
head(txptome)

#Initial visualization of the ACACA motif freq
png("figures/ACACA_motif_txptome.png")
ggplot(txptome, aes(x=ACACA)) + geom_histogram(binwidth=0.05)
dev.off()

#Visualizing highly expressed ACACA motif freq
expressed_ACACA <- subset(txptome, txptome$ACACA >= 3)
png("figures/ACACA_motif_expressed.png")
ggplot(expressed_ACACA, aes(x=ACACA)) + geom_histogram(binwidth=0.05)
dev.off()
ACACA_rnas <- expressed_ACACA$RNA

#Identify the percentage of unspliced txpts that highly express ACACA
unspliced_ACACA_ids <- str_extract(expressed_ACACA$RNA, "unspliced")
unspliced_ACACA <- subset(expressed_ACACA, unspliced_ACACA_ids == "unspliced")
length(unspliced_ACACA$RNA)/length(expressed_ACACA$RNA)

#Identify total unspliced transcriptome
unspliced_txptome_ids <- str_extract(txptome$RNA, "unspliced")
unspliced_txptome <- subset(txptome, unspliced_txptome_ids == "unspliced")
unspliced_txptome_freq <- length(unspliced_txptome$RNA)/length(txptome$RNA)
unspliced_txptome_freq
write.csv(unspliced_txptome,"derived_data/unspliced_txptome.csv", row.names = FALSE)

#Identify total spliced transcriptome
spliced_txptome <- subset(txptome, !(txptome$RNA %in% unspliced_txptome_ids))
write.csv(spliced_txptome, "derived_data/spliced_txptome.csv", row.names = FALSE)

#Identify the transcripts of XIST
XIST_rnas <- str_extract(txptome$RNA, "XIST")
XIST_transcripts <- subset(txptome, XIST_rnas == "XIST")
write.csv(XIST_transcripts, "derived_data/xist_transcripts.csv", row.names = FALSE)
