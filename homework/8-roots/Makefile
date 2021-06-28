CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt hydrogen.png Makefile
	@cat out.txt

hydrogen.png: Hydrogen.txt Makefile
	@echo '\
	set term png ;\
	set title "Radial hydrogen wavefunction" ;\
	set xlabel "r [Bohr radius]" ;\
	set ylabel "$$F(r)$$" ;\
	set key bottom left;\
	set out "$@" ;\
	plot \
	 "Hydrogen.txt" using 1:2 with lines linewidth 2 color black title "Rooted function" \
	, x*exp(-x) with lines linewidth 2 linetype 2 color red title "r*exp(-r)" ;\
	'| pyxplot


out.txt Hydrogen.txt: main Makefile
	./$< > $@

main: main.c qr.o rootfinder.o ode.o

clean: 
	$(RM) main *.o out.* *.png