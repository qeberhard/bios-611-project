library(ggplot2)
library(stringr)
library(readr)
library(tidyr)

#Read in transcriptome data and expression from TPM derived expression
expressed.tpms <- read_csv("derived_data/expressed_tpm_txptome.csv")
txptome1 <- read_csv("source_data/5mers_expressed_txpts_1.csv.gz")
txptome2 <- read_csv("source_data/5mers_expressed_txpts_2.csv.gz")
txptome <- rbind(txptome1, txptome2)
colnames(txptome)[1] <- "RNA"
txptome$RNA <- sub('>>', '', txptome$RNA)
txptome$RNA <- sub('>', '', txptome$RNA)

#Identify total unspliced transcriptome
unspliced_txptome_ids <- str_extract(txptome$RNA, "unspliced")
unspliced_txptome <- subset(txptome, unspliced_txptome_ids == "unspliced")
unspliced_txptome_freq <- length(unspliced_txptome$RNA)/length(txptome$RNA)
#unspliced_txptome_freq
write.csv(unspliced_txptome,"derived_data/unspliced_txptome.csv", row.names = FALSE)

expressed.kmers.txptome <- subset(txptome, (txptome$RNA %in% expressed.tpms$RNA))
write.csv(expressed.kmers.txptome, file="derived_data/expressed_kmers_txptome.csv")




#Identify the transcripts of XIST
XIST_rnas <- str_extract(expressed.kmers.txptome$RNA, "XIST")
XIST_tpms <- subset(expressed.tpms, XIST_rnas == "XIST")

png("figures/xist_expressed_transcripts.png")
ggplot(data = XIST_tpms, aes(x= RNA, y = K562_averaged_total, fill = RNA)) +
  geom_bar(stat = "identity")
dev.off()

XIST_motifs <- subset(expressed.kmers.txptome, XIST_rnas == "XIST")


#Identify the transcript of XACT
XACT_rnas <- str_extract(expressed.kmers.txptome$RNA, "XACT")
XACT_tpms <- subset(expressed.tpms, XACT_rnas == "XACT")

png("figures/xact_expressed_transcripts.png")
ggplot(data = XACT_tpms, aes(x= RNA, y = K562_averaged_total, fill = RNA)) +
  geom_bar(stat = "identity")
dev.off()

XACT_motifs <- subset(expressed.kmers.txptome, XACT_rnas == "XACT")

xist <- subset(expressed.kmers.txptome, (expressed.kmers.txptome$RNA == "XIST_chrX_73820656_73852723_-_ENST00000429829.6"))
xact <- subset(expressed.kmers.txptome, (expressed.kmers.txptome$RNA == "XACT_chrX_113616300_114059121_-_ENST00000468762.3"))

xact.sort <- sort(xact[2:1025], decreasing = TRUE)
xist.sort <- sort(xist[2:1025], decreasing = TRUE)

write.csv(xact.sort[1:10], "derived_data/top_10_xact_5mers.csv", row.names = FALSE)
write.csv(xist.sort[1:10], "derived_data/top_10_xist_5mers.csv", row.names = FALSE)

xist.plot <- t(xist.sort[1:10])
xist.plot <- as.data.frame(xist.plot)
xist.plot$motifs <- rownames(xist.plot)
png("figures/xist_5mer_expression.png")
ggplot(data = xist.plot, aes(x= motifs, y = V1, fill = motifs)) +
  geom_bar(stat = "identity")
dev.off()

xact.plot <- t(xact.sort[1:10])
xact.plot <- as.data.frame(xact.plot)
xact.plot$motifs <- rownames(xact.plot)
png("figures/xact_5mer_expression.png")
ggplot(data = xact.plot, aes(x= motifs, y = V1, fill = motifs)) +
  geom_bar(stat = "identity")
dev.off()



#Identify the transcript of KCNQ1OT1
KCNQ1OT1_rnas <- str_extract(expressed.kmers.txptome$RNA, "KCNQ1OT1")
KCNQ1OT1_tpms <- subset(expressed.tpms, KCNQ1OT1_rnas == "KCNQ1OT1")

png("figures/kcnq1ot1_expressed_transcripts.png")
ggplot(data = KCNQ1OT1_tpms, aes(x= RNA, y = K562_averaged_total, fill = RNA)) +
  geom_bar(stat = "identity")
dev.off()

KCNQ1OT1_motifs <- subset(expressed.kmers.txptome, KCNQ1OT1_rnas == "KCNQ1OT1")

kcnq1ot1 <- subset(expressed.kmers.txptome, (expressed.kmers.txptome$RNA == "KCNQ1OT1_chr11_2608328_2699994_-_ENSG00000269821.1.unspliced"))
kcnq1ot1.sort <- sort(kcnq1ot1[2:1025], decreasing = TRUE)
write.csv(kcnq1ot1.sort[1:10], "derived_data/top_10_kcnq1ot1_5mers.csv", row.names = FALSE)

kcnq1ot1.plot <- t(kcnq1ot1.sort[1:10])
kcnq1ot1.plot <- as.data.frame(kcnq1ot1.plot)
kcnq1ot1.plot$motifs <- rownames(kcnq1ot1.plot)
png("figures/kcnq1ot1_5mer_expression.png")
ggplot(data = kcnq1ot1.plot, aes(x= motifs, y = V1, fill = motifs)) +
  geom_bar(stat = "identity")
dev.off()

#NEAT1 and MALAT1
NEAT1_rnas <- str_extract(expressed.kmers.txptome$RNA, "NEAT1")
NEAT1_motifs <- subset(expressed.kmers.txptome, NEAT1_rnas == "NEAT1")

MALAT1_rnas <- str_extract(expressed.kmers.txptome$RNA, "MALAT1")
MALAT1_motifs <- subset(expressed.kmers.txptome, MALAT1_rnas == "MALAT1")

neat1 <- subset(expressed.kmers.txptome, (expressed.kmers.txptome$RNA == "NEAT1_chr11_65423125_65426499_+_ENST00000645023.1"))
malat1 <- subset(expressed.kmers.txptome, (expressed.kmers.txptome$RNA == "MALAT1_chr11_65504519_65506468_+_ENST00000610851.1"))

neat1.sort <- sort(neat1[2:1025], decreasing = TRUE)
malat1.sort <- sort(malat1[2:1025], decreasing = TRUE)

neat1.plot <- t(neat1.sort[1:10])
neat1.plot <- as.data.frame(neat1.plot)
neat1.plot$motifs <- rownames(neat1.plot)
png("figures/neat1_5mer_expression.png")
ggplot(data = neat1.plot, aes(x= motifs, y = V1, fill = motifs)) +
  geom_bar(stat = "identity")
dev.off()

malat1.plot <- t(malat1.sort[1:10])
malat1.plot <- as.data.frame(malat1.plot)
malat1.plot$motifs <- rownames(malat1.plot)
png("figures/malat1_5mer_expression.png")
ggplot(data = malat1.plot, aes(x= motifs, y = V1, fill = motifs)) +
  geom_bar(stat = "identity")
dev.off()

write.csv(neat1.sort[1:10], "derived_data/top_10_neat1_5mers.csv", row.names = FALSE)
write.csv(malat1.sort[1:10], "derived_data/top_10_malat1_5mers.csv", row.names = FALSE)

#Visualize differences in nucleotide content
nucleotide.df <- as.data.frame(t(matrix(c(sum(str_count(colnames(xist.sort[1:10]), "A")),
  sum(str_count(colnames(xist.sort[1:10]), "T")),
  sum(str_count(colnames(xist.sort[1:10]), "G")),
  sum(str_count(colnames(xist.sort[1:10]), "C")),
  sum(str_count(colnames(xact.sort[1:10]), "A")),
  sum(str_count(colnames(xact.sort[1:10]), "T")),
  sum(str_count(colnames(xact.sort[1:10]), "G")),
  sum(str_count(colnames(xact.sort[1:10]), "C")),
  sum(str_count(colnames(kcnq1ot1.sort[1:10]), "A")),
  sum(str_count(colnames(kcnq1ot1.sort[1:10]), "T")),
  sum(str_count(colnames(kcnq1ot1.sort[1:10]), "G")),
  sum(str_count(colnames(kcnq1ot1.sort[1:10]), "C")),
  sum(str_count(colnames(neat1.sort[1:10]), "A")),
  sum(str_count(colnames(neat1.sort[1:10]), "T")),
  sum(str_count(colnames(neat1.sort[1:10]), "G")),
  sum(str_count(colnames(neat1.sort[1:10]), "C")),
  sum(str_count(colnames(malat1.sort[1:10]), "A")),
  sum(str_count(colnames(malat1.sort[1:10]), "T")),
  sum(str_count(colnames(malat1.sort[1:10]), "G")),
  sum(str_count(colnames(malat1.sort[1:10]), "C"))), nrow=4, ncol=5)))
colnames(nucleotide.df) <- c("A", "T", "G", "C")
rownames(nucleotide.df) <- c("XIST", "XACT", "KCNQ1OT1", "NEAT1", "MALAT1")
nucleotide.df$RNA <- rownames(nucleotide.df)
nuc.long.df <- pivot_longer(nucleotide.df, c("A", "T", "G", "C"), names_to = "nucleotide", values_to = "count" )

png("figures/nucleotide_components_kmers.png")
ggplot(data = nuc.long.df, aes(x= RNA, y = count, fill = nucleotide)) +
  geom_bar(stat = "identity")
dev.off()