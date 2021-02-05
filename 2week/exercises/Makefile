
CC = gcc
CFLAGS = -O -std=gnu11

default: out.txt
	cat out.txt

out.txt: Hello
	./Hello > out.txt

Hello: Hello.o
	$(CC) -o Hello Hello.o

Hello.o: Hello.c
	$(CC) $(CFLAGS) -c Hello.c

clean: 
	rm -f Hello.o Hello out.txt
