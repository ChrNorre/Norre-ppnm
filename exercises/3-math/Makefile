CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: mathex
	./mathex > out.txt

mathex: mathex.o
	$(CC) $(LDLIBS) mathex.o -o mathex

mathex.o: mathex.c
	$(CC) $(CFLAGS) -c mathex.c -o mathex.o

clean: 
	rm -f mathex mathex.o out.txt
