all: simple_figures data summary_figures supplemental_figures

figures: simple_figures summary_figures supplemental_figures

# make figures that require few simulations
simple_figures: 1BCDE.eps 1FGH_gamma.eps 1FGH_uniform.eps 1FGH_trend.eps 2ABC.eps \
	4B.eps 4G.eps 4B_trend.eps 4CD.eps 4HI.eps 4CD_trend.eps

# make data for figures that require many rounds of simulations
data: 3B.csv 3D.csv 3F.csv 4E.csv 4J.csv 4E_trend.csv 5B.csv 5D.csv \
	5B_trend.csv 5D_trend.csv 6B.csv S1B.csv
# Note: this is better than calling all functions in python because
#       memory allocated for python is cleared after each rule

# make figures that require many rounds of simulations
# Note: data must already be made when calling summary_figures
summary_figures: 3B.eps 3D.eps 3F.eps 4E.eps 4J.eps 4E_trend.eps \
	5B.eps 5D.eps 5B_trend.eps 5D_trend.eps 6B.eps 6B_dists.eps

# make supplementary figures which require many rounds of simulations
# Note: data must already be made when calling summary_figures
supplemental_figures: S1B.eps S2A.eps S2ABC.eps S2A_trend.eps S2ABC_trend.eps \
	S2C.eps S2D.eps S3A.eps S3B.eps S3A_trend.eps S3F.eps S3C.eps \
	S3D.eps S4A.eps S4B.eps



3B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3B()"

3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3D()"

3F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3F()"

4E.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4E()"

4J.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4J()"

4E_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4E_trend()"

5B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5B()"

5D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5D()"

5B_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5B_trend()"

5D_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5D_trend()"

6B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_6B()"

6B_dists.eps:
	cd scripts && python -c "from figures import fig_6B_dists; fig_6B_dists()"



S1B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S1B()"

S2A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2A()"

S2A_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2A_trend()"

S2ABC.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2ABC()"

S2ABC_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2ABC_trend()"

S2C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2C()"

S2D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2D()"

S3A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3A()"

S3A_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3A_trend()"

S3B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3B()"

S3F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3F()"

S3C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3C()"

S3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3D()"

S4A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S4A()"

S4B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S4B()"



3B.csv:
	cd scripts && python -c "from make_data import data_3B; data_3B()"

3D.csv:
	cd scripts && python -c "from make_data import data_3D; data_3D()"

3F.csv:
	cd scripts && python -c "from make_data import data_3F; data_3F()"

4E.csv:
	cd scripts && python -c "from make_data import data_4E; data_4E()"

4J.csv:
	cd scripts && python -c "from make_data import data_4J; data_4J()"

4E_trend.csv:
	cd scripts && python -c "from make_data import data_4E_trend; data_4E_trend()"

5B.csv:
	cd scripts && python -c "from make_data import data_5B; data_5B()"

5D.csv:
	cd scripts && python -c "from make_data import data_5D; data_5D()"

5B_trend.csv:
	cd scripts && python -c "from make_data import data_5B_trend; data_5B_trend()"

5D_trend.csv:
	cd scripts && python -c "from make_data import data_5D_trend; data_5D_trend()"

6B.csv:
	cd scripts && python -c "from make_data import data_6B; data_6B()"

S1B.csv:
	cd scripts && python -c "from make_data import data_S1B; data_S1B()"



1BCDE.eps:
	cd scripts && python -c "from figures import fig_1BCDE; fig_1BCDE()"

1FGH_gamma.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('gamma')"

1FGH_uniform.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('uniform')"

1FGH_trend.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('trend')"

2ABC.eps:
	cd scripts && python -c "from figures import fig_2ABC; fig_2ABC()"

4B.eps:
	cd scripts && python -c "from figures import fig_4B; fig_4B()"

4G.eps:
	cd scripts && python -c "from figures import fig_4G; fig_4G()"

4B_trend.eps:
	cd scripts && python -c "from figures import fig_4B_trend; fig_4B_trend()"

4CD.eps:
	cd scripts && python -c "from figures import fig_4CD; fig_4CD()"

4HI.eps:
	cd scripts && python -c "from figures import fig_4HI; fig_4HI()"

4CD_trend.eps:
	cd scripts && python -c "from figures import fig_4CD_trend; fig_4CD_trend()"



# delete simulated data and figures
clean: clean_data clean_figures

clean_data:
	cd data_simulated && rm -f *

clean_figures:
	cd figures && rm -rf *

.PHONY: clean clean_data clean_figures
