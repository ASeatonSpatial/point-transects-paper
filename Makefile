default:

run-all:
				Rscript R/fit_model.R
				Rscript R/eval_model.R
				Rscript R/eval_spde.R
				Rscript R/posterior_N.R
				Rscript R/excursions.R
				Rscript R/make_study_area_fig.R
