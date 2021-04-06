CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm


default: exp_pyx.pdf exp_diff_pyx.pdf exp_other.pdf

exp_other.pdf: exp_other.gpi outdata.txt
	pyxplot $<

exp_pyx.pdf: exp_pyx.gpi outdata.txt
	pyxplot $<

exp_diff_pyx.pdf: exp_diff_pyx.gpi outdata.txt
	pyxplot $<

outdata.txt: main
	./$< > $@

