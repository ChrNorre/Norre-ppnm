CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt SIR.png Makefile
	@cat out.txt

out.txt SIR_path.txt: main three_body
	./$< > $@
	./three_body >> $@

SIR.png: SIR_path.txt Makefile
	@echo '\
	set term png ;\
	set title "SIR, Tc=1, Tr=14" ;\
	set xlabel "Days" ;\
	set ylabel "Number of people" ;\
	set key bottom left;\
	set out "$@" ;\
	plot \
	 "SIR_path.txt" using 1:(($$2+$$3+$$4)/5800000) with lines color black title "Suseptible below" \
	, "SIR_path.txt" using 1:(($$3+$$4)/5800000) with lines color red title "Infected below" \
	, "SIR_path.txt" using 1:(($$4)/5800000) with lines color blue title "Removed below" ;\
	'| pyxplot

main: main.c ode.c

three_body: three_body.c ode.c

clean: 
	$(RM) main three_body *.o *.txt *.png