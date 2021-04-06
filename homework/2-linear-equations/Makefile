CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out_tests.txt times.png Makefile


times.png: out_gsl.txt out_qr.txt
	echo '\
	set term png ;\
	set title "QR-decomposition times as function of matrix size" ;\
	set xlabel "matrix size N" ;\
	set ylabel "time, sec" ;\
	set key left top ;\
	set out "$@" ;\
	a=1 ;\
	f(x)=a*x**3 ;\
	fit f(x) "out_qr.txt" via a ;\
	b=1 ;\
	g(x)=b*x**3 ;\
	fit g(x) "out_gsl.txt" via b ;\
	plot \
	 "out_qr.txt" with linespoints title "my qr decomp" \
	, "out_gsl.txt" with linespoints title "gsl qr decomp" \
	, f(x) with lines title "qr fit:  const*NÂ³" \
	, g(x) with lines title "gsl fit:  const*NÂ³" ;\
	'| gnuplot


out_tests.txt: main Makefile
	./main > out_tests.txt

out_gsl.txt out_qr.txt: qr_time gsl_time Makefile
	>out_gsl.txt
	for n in $$(seq 300 25 700);\
	do \time --format "$$n %U" --append --output out_gsl.txt ./gsl_time $$n ;\
	done

	>out_qr.txt
	for n in $$(seq 300 25 700);\
	do \time --format "$$n %U" --append --output out_qr.txt ./qr_time $$n ;\
	done
#i would like to have the times in the same file, but \time has an automatic newline, and i dont know how to remove it

qr_time: qr_time.o qr.o

gsl_time: gsl_time.o

main: main.o qr.o tests.o

clean: 
	$(RM) main gsl_time qr_time *.o *.txt *.png