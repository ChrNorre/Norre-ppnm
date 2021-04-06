CC = gcc
CFLAGS = -O -Wall -std=gnu11 -fwrapv
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon
	./epsilon > out.txt

epsilon: epsilon.o equal.o
	$(CC) $(LDFLAGS) epsilon.o equal.o -o epsilon $(LDLIBS)

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c -o epsilon.o

equal.o: equal.c
	$(CC) $(CFLAGS) -c equal.c -o equal.o

clean: 
	$(RM) epsilon epsilon.o out.txt


