CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt err_est.png Makefile
	@cat out.txt

out.txt: main
	./$< > $@

err_est.png: out.data.txt Makefile
	@echo '\
	set term png ;\
	set title "Scaling of error estimate as function of N" ;\
	set xlabel "Log(N)" ;\
	set ylabel "Log(error estimate)" ;\
	set key outside ;\
	set out "$@" ;\
	a=-0.5 ;\
	f(x)=a*x ;\
	fit f(x) "$<" using 1:2 via a ;\
	b=-0.5 ;\
	g(x)=b*x ;\
	fit g(x) "$<" using 1:3 via b ;\
	plot \
	 "$<" using 1:2 with linespoints title "plainmc" \
	, "$<" using 1:3 with linespoints title "quasimc" \
	, f(x) with lines title "plain fit:  %g*log(N)"%(a) \
	, g(x) with lines title "quasi fit:  %g*log(N)"%(b) \
	, -0.5*x with lines title "1/sqrt(N),    -0.5*log(N)" ;\
	'| pyxplot

main: main.c montecarlo.o

clean: 
	$(RM) main *.o out.* *.png