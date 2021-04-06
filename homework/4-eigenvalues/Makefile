CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt quantum.png Makefile

out.txt: main
	./main > $@
	cat $@

quantum.png: data.txt Makefile
	@echo '\
	set term png ;\
	set title "Plot of first 4 eigenfunctions for particle in a box" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set xrange [0:1] ;\
	set key outside ;\
	set out "$@" ;\
	a=0.3086067 ;\
	plot \
	 "data.txt" using 1:2 with points color black title "1. eigenstate" \
	, a*sin(x*(2*pi)/2) with lines color black title "sin(x/2)" \
	, "data.txt" using 1:3 with points color red title "2. eigenstate" \
	, -a*sin(2*x*(2*pi)/2) with lines color red title "sin(2*x/2)" \
	, "data.txt" using 1:4 with points color blue title "3. eigenstate" \
	, a*sin(3*x*(2*pi)/2) with lines color blue title "sin(3*x/2)" \
	, "data.txt" using 1:5 with points color purple title "4. eigenstate" \
	, -a*sin(4*x*(2*pi)/2) with lines color purple title "sin(4*x/2)" ;\
	'| pyxplot

data.txt: quantum
	./quantum

quantum: quantum.o jacobi.o

main: main.o jacobi.o

clean: 
	$(RM) main quantum  *.o *.txt *.png