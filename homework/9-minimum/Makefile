CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt Makefile
	@cat out.txt

out.txt: main Makefile
	./$< > $@

main: main.c minimum.o

clean: 
	$(RM) main *.o *.txt *.png