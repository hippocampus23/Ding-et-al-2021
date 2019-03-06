all: simple_figures data summary_figures supplemental_figures

figures: simple_figures summary_figures supplemental_figures

# make figures that require few simulations
simple_figures: 1ABCD.eps 1FGH_gamma.eps 1FGH_uniform.eps 1FGH_trend.eps 2BDF.eps \
	4B.eps 4C.eps 4F.eps 4G.eps 4J.eps 4K.eps

# make data for figures that require many rounds of simulations
data: 3B.csv 3D.csv 3F.csv 4D.csv 4H.csv 4L.csv 5B.csv 5D.csv \
	5B_trend.csv 5D_trend.csv 6B.csv 7C.csv S1B.csv
# Note: this is better than calling all functions in python because
#       memory allocated for python is cleared after each rule

# make figures that require many rounds of simulations
# Note: data must already be made when calling summary_figures
summary_figures: 3B.eps 3D.eps 3F.eps 4D.eps 4H.eps 4L.eps \
	5B.eps 5D.eps 6B.eps 7C.eps 7D.eps

# make supplementary figures which require many rounds of simulations
# Note: data must already be made when calling summary_figures
supplemental_figures: S1B.eps S2A.eps S2B.eps S2E.eps S2F.eps \
	S2C.eps S2D.eps S3A.eps S3B.eps S3C.eps \
	S3D.eps S4A.eps S4B.eps

3B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3B(); quit()"

3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3D(); quit()"

3F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_3F(); quit()"

4D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4D(); quit()"

4H.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4H(); quit()"

4L.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_4L(); quit()"

5B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5B(); quit()"

5D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_5D(); quit()"

6B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_6B(); quit()"

7C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_7C(); quit()"

7D.eps:
	cd scripts && python -c "from figures import fig_7D; fig_7D(); quit()"



S1B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S1B(); quit()"

S2A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2A(); quit()"

S2E.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2E(); quit()"

S2B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2B(); quit()"

S2F.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2F(); quit()"

S2C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2C(); quit()"

S2D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S2D(); quit()"

S3A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3A(); quit()"

S3B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3B(); quit()"

S3C.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3C(); quit()"

S3D.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S3D(); quit()"

S4A.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S4A(); quit()"

S4B.eps:
	cd scripts && Rscript -e "source('figures.R'); fig_S4B(); quit()"



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

6B.csv:
	cd scripts && python -c "from make_data import data_6B; data_6B()"

7C.csv:
	cd scripts && python -c "from make_data import data_7C; data_7C()"

S1B.csv:
	cd scripts && python -c "from make_data import data_S1B; data_S1B()"



1ABCD.eps:
	cd scripts && python -c "from figures import fig_1ABCD; fig_1ABCD(); quit()"

1FGH_gamma.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('gamma'); quit()"

1FGH_uniform.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('uniform'); quit()"

1FGH_trend.eps:
	cd scripts && python -c "from figures import fig_1FGH; fig_1FGH('trend'); quit()"

2BDF.eps:
	cd scripts && python -c "from figures import fig_2BDF; fig_2BDF(); quit()"

4B.eps:
	cd scripts && python -c "from figures import fig_4B; fig_4B(); quit()"

4F.eps:
	cd scripts && python -c "from figures import fig_4F; fig_4F(); quit()"

4J.eps:
	cd scripts && python -c "from figures import fig_4J; fig_4J(); quit()"

4C.eps:
	cd scripts && python -c "from figures import fig_4C; fig_4C(); quit()"

4G.eps:
	cd scripts && python -c "from figures import fig_4G; fig_4G(); quit()"

4K.eps:
	cd scripts && python -c "from figures import fig_4K; fig_4K(); quit()"



# delete simulated data and figures
clean: clean_data clean_figures

clean_data:
	cd data_simulated && rm -f *

clean_figures:
	cd figures && rm -rf *

.PHONY: clean clean_data clean_figures
