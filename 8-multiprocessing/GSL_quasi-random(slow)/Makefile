CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)



default: pthreads openMP piestimate fit.png Makefile
	./pthreads
	./openMP	


outpi.txt: piestimate Makefile
	./piestimate $$(seq 1 0.1 7) > outpi.txt


fit.png: outpi.txt Makefile
	echo '\
	f(x) = a1*x;\
	g(x) = -0.5*x;\
	fit f() '$<' via a1;\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set xlabel "log$$_{10}(N)$$";\
	set ylabel "log$$_{10}|\pi - pi(N)|$$";\
	set title "Monte Carlo for $$\pi$$";\
	plot \
	 "$<" with points title "Pi estimates"\
	,f(x) with line title "fitted function"\
	,g(x) with line title "log$$_{10}(1/\sqrt{N})$$"\
	' | tee test.gpi | pyxplot


pthreads: pthreads.c pi.c

openMP: openMP.c pi.c

piestimate: piestimate.c pi.c

clean: 
	$(RM) fit.png openMP outpi.txt piestimate pthreads test.gpi pi
