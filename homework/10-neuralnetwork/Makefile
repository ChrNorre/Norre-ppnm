CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt interp.png Makefile
	

interp.png: out.txt Makefile
	@echo '\
	set term png ;\
	set title "Interpolation via neural network" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set key outside;\
	set out "$@" ;\
	plot \
	 "out.txt" index 1 using 1:2 with lines linewidth 2 color black title "Interpolation" \
	, "out.txt" index 1 using 1:3 with lines linewidth 2 color blue title "Integral from interp" \
	, "out.txt" index 1 using 1:4 with lines linewidth 2 color blue title "Derivative from interp" \
	, "out.txt" index 1 using 1:5 with lines linetype 2 linewidth 2 color red title "Derivative of point function" \
	, "out.txt" index 0 using 1:2 with points pointtype 3 linewidth 2 color red title "points" ;\
	'| pyxplot



out.txt: main Makefile
	./$< > $@
	
main: main.c ann.o minimum.o


clean: 
	$(RM) main *.o out.* *.png