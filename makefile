.PHONY: clean

clean:
	rm -rf figures
	rm -rf derived_data
	rm -rf .created-dirs

.created-dirs:
	mkdir -p figures
	mkdir -p derived_data
	touch .created-dirs

#Expression exploration
figures/total_transcriptome_expression_initial.png\
 figures/log_transcript_threshold_expression.png\
 derived_data/expressed_tpm_txptome.csv\
 figures/splicing_diff_expressed.png:\
  source_data/K562_averaged_table.tsv\
  .created-dirs
  	Rscript prelim_analysis/total_rnaseq_expression.R


#k-mer analysis
derived_data/unspliced_txptome.csv\
 derived_data/expressed_kmers_txptome.csv\
 derived_data/top_10_xist_5mers.csv\
 derived_data/top_10_xact_5mers.csv\
 derived_data/top_10_kcnq1ot1_5mers.csv\
 derived_data/top_10_neat1_5mers.csv\
 derived_data/top_10_malat1_5mers.csv:\
  derived_data/expressed_tpm_txptome.csv\
  source_data/5mers_expressed_txpts_1.csv.gz\
  source_data/5mers_expressed_txpts_2.csv.gz\
  source_data/gencode.v39.basic.annotation.gtf.gz\
  .created-dirs
  	Rscript prelim_analysis/initial_kmer_analysis.R

report.pdf: source_data/K562_averaged_table.tsv\
 figures/total_transcriptome_expression_initial.png\
 figures/log_transcript_threshold_expression.png\
 figures/splicing_diff_expressed.png\
 derived_data/expressed_kmers_txptome.csv\
 derived_data/expressed_tpm_txptome.csv\
 derived_data/top_10_xist_5mers.csv\
 derived_data/top_10_xact_5mers.csv
 derived_data/top_10_kcnq1ot1_5mers.csv\
 derived_data/top_10_malat1_5mers.csv\
 derived_data/top_10_neat1_5mers.csv\
 figures/nucleotide_components_kmers.png
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"pdf_document\")"