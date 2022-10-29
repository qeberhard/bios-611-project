.PHONY: clean

clean:
	rm -rf figures
	rm -rf derived_data
	rm -rf .created-dirs

.created-dirs:
	mkdir -p figures
	mkdir -p derived_data
	touch .created-dirs

<<<<<<< HEAD
=======
#Develop the preliminary visualization
>>>>>>> 09fd905779f3dac18ae1b97f21b50ed8fdb6ced3
figures/ACACA_motif_txptome.png\
 figures/ACACA_motif_expressed.png\
 derived_data/unspliced_txptome.csv\
 derived_data/xist_motifs.csv:\
<<<<<<< HEAD
  source_data/5mers_expressed_txpts.csv.gz\
=======
  source_data/5mers_expressed_txpts_1.csv.gz\
  source_data/5mers_expressed_txpts_2.csv.gz\
>>>>>>> 09fd905779f3dac18ae1b97f21b50ed8fdb6ced3
  .created-dirs
	Rscript prelim_analysis/initial_data_exploration.R

report.pdf: figures/ACACA_motif_txptome.png figures/ACACA_motif_expressed.png
<<<<<<< HEAD
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"pdf_document\")"
=======
	R -e "rmarkdown::render(\"report.Rmd\", output_format=\"pdf_document\")
>>>>>>> 09fd905779f3dac18ae1b97f21b50ed8fdb6ced3
