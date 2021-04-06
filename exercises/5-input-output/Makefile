CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

make_out: mainA mainB mainC input.txt
	./mainA $$(cat input.txt) > outA.txt
	./mainB < input.txt > outB.txt
	./mainC input.txt outC.txt
		


mainA: mainA.c
mainB: mainB.c
mainC: mainC.c

input.txt: Makefile
	cat /dev/null > input.txt
	for x in $$(seq 0 0.3 3); do echo "$$x" >> input.txt; done

clean: 
	$(RM) out* mainA mainB mainC in*

out:
	cat outA.txt
	cat outB.txt
	cat outC.txt

