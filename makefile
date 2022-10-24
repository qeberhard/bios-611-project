.PHONY: clean
.PHONY: visualization

clean:
	rm -rf figures
	rm -rf derived_data
	rm -rf .created-dirs

.created-dirs:
	mkdir -p figures
	mkdir -p derived_data
	touch .created-dirs

#Develop the preliminary visualization
#save figs as something and change .Rmd to .R
#save written part into a write up document
figures/ACACA_motif_txptome.png\
 figures/ACACA_motif_expressed.png\
 derived_data/unspliced_txptome.csv\
 derived_data/xist_motifs.csv:\
  source_data/5mers_expressed_txpts.csv.gz\
  .created-dirs
	Rscript prelim_analysis/initial_data_exploration.R

#figures/reduced-demo-gender_female-roc.png\
 figures/reduced-demo-ethnicity_white-roc.png\
 figures/reduced-demo-married-roc.png\
 figures/reduced_demographic_projection.png:\
  models/reduced-demographics-enc\
  models/reduced-demographics-ae\
  derived_data/reduced-demographics-one-hot.csv
#  .created-dirs
#	Rscript reduced-demographics-plots.R

writeup.pdf: figures/ACACA_motif_txptome.png figures/ACACA_motif_expressed.png
	pdflatex writeup.tex
