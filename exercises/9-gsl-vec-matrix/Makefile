CFLAGS = -Wall -Werror -O1 -std=gnu11
LDLIBS = -lm
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: main hilbert
	./main
	./hilbert

clean:
	$(RM) lecture main hilbert 