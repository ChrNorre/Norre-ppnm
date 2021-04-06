CFLAGS = -o -w -std=gnu99
#CFLAGS += 'gsl-config --cflags' 		#old way which does not work for me
#use $(shell ______) or $$(______)
CFLAGS += $(shell gsl-config --cflags)
#LDLIBS += $(gsl-config --libs)			#$ is a make variable call
LDLIBS += $$(gsl-config --libs)			


#we likely want output to be a .txt file
#all: output ; cat $@				#$< is first prerequisite, $@ is target (cant cat target)
all: output.txt ; cat $<				
#output: main : $@ > $< 			#used colon. wrong order of $@ and $<
output.txt: main ; ./$< > $@
#clean: ; $RM main output			#variable is called RM, but $ only takes the first character, thus use ()
clean: ; $(RM) main output
