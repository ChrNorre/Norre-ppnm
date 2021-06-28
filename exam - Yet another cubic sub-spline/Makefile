CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: cub_sub_example1.png cub_sub_example2.png Makefile

cub_sub_example2.png: out.data2.txt Makefile
	@echo '\
	data = "out.data2.txt" ;\
	set term png ;\
	set title "Sub-Spline" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set key top left;\
	set out "$@" ;\
	set multiplot ;\
	set width 8 ;\
	plot \
	 data index 1 using 1:2 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:3 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:4 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:2 with points pointtype 3 linewidth 2 color red title "points" ;\
	set origin 10,0 ;\
	set title "Derivative of splines" ;\
	set xlabel "x" ;\
	set ylabel "dy/dx" ;\
	set key bottom left;\
	plot \
	 data index 1 using 1:5 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:6 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:7 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:(0*$$2) with points pointtype 3 linewidth 2 color red title "points" ;\
	set origin 20,0 ;\
	set title "Integral of splines" ;\
	set xlabel "x" ;\
	set ylabel "Integral of y" ;\
	set key top left;\
	plot \
	 data index 1 using 1:8 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:9 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:10 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:(0*$$2) with points pointtype 3 linewidth 2 color red title "points" ;\
	'| pyxplot

cub_sub_example1.png: out.data1.txt Makefile
	@echo '\
	data = "out.data1.txt" ;\
	set term png ;\
	set title "Sub-Spline" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set key top left;\
	set out "$@" ;\
	set multiplot ;\
	set width 8 ;\
	plot \
	 data index 1 using 1:2 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:3 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:4 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:2 with points pointtype 3 linewidth 2 color red title "points" ;\
	set origin 10,0 ;\
	set title "Derivative of splines" ;\
	set xlabel "x" ;\
	set ylabel "dy/dx" ;\
	plot \
	 data index 1 using 1:5 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:6 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:7 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:(0*$$2) with points pointtype 3 linewidth 2 color red title "points" ;\
	set origin 20,0 ;\
	set title "Integral of splines" ;\
	set xlabel "x" ;\
	set ylabel "Integral of y" ;\
	plot \
	 data index 1 using 1:8 with lines linewidth 2 color black title "My sub-spline" \
	, data index 1 using 1:9 with lines linewidth 2 linetype 2 color purple title "gsl-cubic" \
	, data index 1 using 1:10 with lines linewidth 2 linetype 2 color red title "gsl-akima" \
	, data index 0 using 1:(0*$$2) with points pointtype 3 linewidth 2 color red title "points" ;\
	'| pyxplot


out.data1.txt out.data2.txt: main Makefile
	./main
	

main: main.c spline.o


clean: 
	$(RM) main *.o out.* *.png