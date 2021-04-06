CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)


make_all: erf.png gam.png loggam.png complexgam.png


erf.png: erfout.txt erf_table.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "erf(x)";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"erf_table.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

gam.png: gammaout.txt gamma_table.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "Gamma(x)";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "homemade Gamma"\
	,"gamma_table.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

loggam.png: loggammaout.txt log_gamma_table.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "ln(Gamma(x))";\
	set title "Log of Gamma function";\
	plot \
	 "$<" using 1:2 with line title "lgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-ln-gamma"\
	,"$<" using 1:4 with line title "log of homemade Gamma"\
	,"log_gamma_table.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot


complexgam.png: xyz_out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set nokey;\
	set xlabel "real";\
	set ylabel "imag";\
	set zlabel "abs(Gamma(z))";\
	set title "Absolute value of complex Gamma function";\
	set xrange [-4.2:4.2];\
	set xtics -4,1,4 ; set mxtics 0.5;\
	set yrange [-4.2:4.2];\
	set ytics -4,1,4 ; set mytics 0.5;\
	set zrange [0:6];\
	set ztics 1 ; set mztics 0.5;\
	set size 8 square;\
	set view -40, 20;\
	set grid;\
	set samples grid 120x120;\
	plot 3d "xyz_out.txt" using 1:2:3 with surface fillc red;\
	' | tee test.gpi | pyxplot


erfout.txt: plotErf
	./$< > $@

gammaout.txt: plotGam
	./$< > $@

loggammaout.txt: plotLogGam
	./$< > $@


xyz_out.txt: complexGamma
	./$<


clean:
	$(RM) *.png plotErf plotGam plotLogGam complexGamma *out.txt
