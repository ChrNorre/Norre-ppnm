CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt decay.png Makefile

decay.png: data.txt fit.txt
	echo '\
	set term png ;\
	set title "plot of data and exponential fit" ;\
	set xlabel "time [days]" ;\
	set ylabel "log(activity)" ;\
	set key left bottom ;\
	set out "$@" ;\
	plot \
	 "data.txt" using 1:2:3 with yerrorbars title "data" \
	, "fit.txt" using 1:3 with lines linetype 2 title "lower uncertainty in lambda" \
	, "fit.txt" using 1:4 with lines linetype 2 title "upper uncertainty in lambda" \
	, "fit.txt" using 1:2 with lines title "fit: ln(y)=ln(a)-lambda*t";\
	'| pyxplot

out.txt: main Makefile
	./main > out.txt
	cat out.txt

main: main.o qr.o least_squares_fit.o

clean: 
	$(RM) main  *.o out.txt decay.png