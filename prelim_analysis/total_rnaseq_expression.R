library(ggplot2)
library(stringr)
library(readr)

#First: determine which genes are expressed!
#read in data
averaged.tpm <- read.csv("source_data/K562_averaged_table.tsv", sep="\t")

#visualize the expression of the averaged.counts
png("figures/total_transcriptome_expression_initial.png")
ggplot(averaged.tpm, aes(x=K562_averaged_total)) + 
  geom_histogram(color = "black", fill="lightblue", bins = 10000) + ylim(0,15) +
  xlab("Gene Expression (TPM)") + ylab("Counts")
dev.off()

minimum.expressed <- subset(averaged.tpm, (averaged.tpm$K562_averaged_total > 0))

#A common cutoff is a TPM of 0.25
log.threshold <- log(1.25)

minimum.expressed$log.K562_averaged_total <- log(minimum.expressed$K562_averaged_total + 1)

png("figures/log_transcript_threshold_expression.png")
ggplot(minimum.expressed, aes(x=log.K562_averaged_total)) + 
  geom_histogram(color = "black", fill="lightblue", bins = 700, binwidth=0.5) +
  geom_vline(xintercept=log.threshold, linetype="dotted") +
  xlab("log(TPM Expression)") + ylab("Counts")
dev.off()

expressed.tpm <- subset(minimum.expressed, (minimum.expressed$K562_averaged_total >= 0.25))
not.expressed.tpm <- subset(averaged.tpm, !(averaged.tpm$K562_averaged_total >= 0.25))
write.csv(expressed.tpm, "derived_data/expressed_tpm_txptome.csv")


# Which genes are not expressed vs. which are?
# with regards to spliced vs. unspliced:
unspliced.express.sum <- sum(grepl("unspliced", expressed.tpm$RNA, fixed=TRUE))
unspliced.not.express.sum <- sum(grepl("unspliced", not.expressed.tpm$RNA, fixed=TRUE))
spliced.express.sum <- sum(!(grepl("unspliced", expressed.tpm$RNA, fixed=TRUE)))
spliced.not.express.sum <- sum(!(grepl("unspliced", not.expressed.tpm$RNA, fixed=TRUE)))
unspliced.total.sum <- sum(grepl("unspliced", averaged.tpm$RNA, fixed=TRUE))
spliced.total.sum <- sum(!(grepl("unspliced", averaged.tpm$RNA, fixed=TRUE)))

splice.express.df <- as.data.frame(matrix(
  c(spliced.express.sum, unspliced.express.sum, spliced.not.express.sum, unspliced.not.express.sum, spliced.total.sum, unspliced.total.sum,
    "spliced", "unspliced", "spliced", "unspliced", "spliced", "unspliced",
    "expressed", "expressed", "unexpressed", "unexpressed", "total", "total",
    "lightblue", "coral", "lightblue", "coral", "lightblue", "coral"),
              nrow = 6, ncol = 4))

colnames(splice.express.df) <- c("Value", "Spliced", "Expression", "Color")

splice.express.df$Color
splice.express.df$Spliced
splice.express.df$Expression
splice.express.df$Value

png("figures/splicing_diff_expressed.png")
p <- ggplot(data=splice.express.df, aes(x=Expression,y=as.numeric(Value),fill=Spliced)) +
  geom_bar(position="dodge",stat="identity") + 
  xlab("Expression") + ylab("Number of Transcripts")
p
dev.off()