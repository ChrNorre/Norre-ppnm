CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: cub_sub_example1.png cub_sub_example2.png Makefile


#It took me a couple hours to figure out that pyxplot gets VERY angry,
# if the scipt file is not .ppl format. I wouldnt think i matters since it gets concatenated
# and piped to pyxplot in what i assume is a string format
cub_sub_example2.png: plot.ppl out.data2.txt Makefile
	cat $< |sed "s/num/2/g" | pyxplot 

cub_sub_example1.png: plot.ppl out.data1.txt Makefile
	cat $< |sed "s/num/1/g" | pyxplot 


out.data1.txt out.data2.txt: main Makefile
	./main
	

main: main.c spline.o


clean: 
	$(RM) main *.o out.* *.png