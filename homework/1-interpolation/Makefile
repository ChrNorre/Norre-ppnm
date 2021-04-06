CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: linear.png quadratic.png cubic.png

linear.png: linear.gpi linear.txt
	pyxplot $<

quadratic.png: quadratic.gpi quadratic.txt
	pyxplot $<

cubic.png: cubic.gpi cubic.txt
	pyxplot $<


main: main.o linear.o quadratic.o cubic.o

linear.txt quadratic.txt cubic.txt: main
	./$<



clean: 
	$(RM) main main.o linear.o quadratic.o cubic.o *.txt *.png