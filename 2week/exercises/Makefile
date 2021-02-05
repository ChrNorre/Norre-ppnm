

default: out.txt
	cat out.txt

out.txt: Hello
	./Hello > out.txt

Hello: Hello.o
	gcc -o Hello Hello.o

Hello.o: Hello.c
	gcc -c Hello.c

clean: 
	rm -f Hello.o Hello out.txt
