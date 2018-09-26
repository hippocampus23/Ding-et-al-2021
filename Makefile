all: simple_figures data summary_figures supplemental_figures

figures: simple_figures summary_figures supplemental_figures

# make figures that require few simulations
simple_figures: 1ABCD.eps 1FGH_gamma.eps 1FGH_uniform.eps 1FGH_trend.eps 2BDF.eps \
	4B.eps 4C.eps 4F.eps 4G.eps 4J.eps 4K.eps

# make data for figures that require many rounds of simulations
data: 3B.csv 3D.csv 3F.csv 4D.csv 4H.csv 4L.csv 5B.csv 5D.csv \
	5B_trend.csv 5D_trend.csv 6C.csv S1B.csv
# Note: this is better than calling all functions in python because
#       memory allocated for python is cleared after each rule

# make figures that require many rounds of simulations
# Note: data must already be made when calling summary_figures
summary_figures: 3B.eps 3D.eps 3F.eps 4D.eps 4H.eps 4L.eps \
	5B.eps 5D.eps 5B_trend.eps 5D_trend.eps 6C.eps 6C_dists.eps

# make supplementary figures which require many rounds of simulations
# Note: data must already be made when calling summary_figures
supplemental_figures: S1B.eps S2A.eps S2B.eps S2E.eps S2F.eps \
	S2C.eps S2D.eps S3A.eps S3B.eps S3A_trend.eps S3B_trend.eps S3C.eps \
	S3D.eps S3C_trend.eps S3D_trend.eps S4A.eps S4B.eps



3B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3B()"

3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3D()"

3F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3F()"

4D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4D()"

4H.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4H()"

4L.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4L()"

5B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5B()"

5D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5D()"

5B_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5B_trend()"

5D_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5D_trend()"

6C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_6C()"

6C_dists.eps:
	cd scripts && python -c "from figures import fig_6C_dists; fig_6C_dists()"



S1B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S1B()"

S2A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2A()"

S2E.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2E()"

S2B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2B()"

S2F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2F()"

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

S3B_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3B_trend()"

S3C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3C()"

S3C_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3C_trend()"

S3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3D()"

S3D_trend.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3D_trend()"

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

4D.csv:
	cd scripts && python -c "from make_data import data_4D; data_4D()"

4H.csv:
	cd scripts && python -c "from make_data import data_4H; data_4H()"

4L.csv:
	cd scripts && python -c "from make_data import data_4L; data_4L()"

5B.csv:
	cd scripts && python -c "from make_data import data_5B; data_5B()"

5D.csv:
	cd scripts && python -c "from make_data import data_5D; data_5D()"

5B_trend.csv:
	cd scripts && python -c "from make_data import data_5B_trend; data_5B_trend()"

5D_trend.csv:
	cd scripts && python -c "from make_data import data_5D_trend; data_5D_trend()"

6C.csv:
	cd scripts && python -c "from make_data import data_6C; data_6C()"

S1B.csv:
	cd scripts && python -c "from make_data import data_S1B; data_S1B()"



1ABCD.eps:
	cd scripts && python -c "from figures import fig_1ABCD; fig_1ABCD()"

1FGH_gamma.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('gamma')"

1FGH_uniform.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('uniform')"

1FGH_trend.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('trend')"

2BDF.eps:
	cd scripts && python -c "from figures import fig_2BDF; fig_2BDF()"

4B.eps:
	cd scripts && python -c "from figures import fig_4B; fig_4B()"

4F.eps:
	cd scripts && python -c "from figures import fig_4F; fig_4F()"

4J.eps:
	cd scripts && python -c "from figures import fig_4J; fig_4J()"

4C.eps:
	cd scripts && python -c "from figures import fig_4C; fig_4C()"

4G.eps:
	cd scripts && python -c "from figures import fig_4G; fig_4G()"

4K.eps:
	cd scripts && python -c "from figures import fig_4K; fig_4K()"



# delete simulated data and figures
clean: clean_data clean_figures

clean_data:
	cd data_simulated && rm -f *

clean_figures:
	cd figures && rm -rf *

.PHONY: clean clean_data clean_figures
