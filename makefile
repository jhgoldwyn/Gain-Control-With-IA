all: ga gaCable

ga: ga.o memory.o 
	gcc -g ga.o memory.o -o ga -lgsl -lgslcblas -lm
ga.o: ga.c memory.h 
	gcc -c -g ga.c

gaCable: gaCable.o memory.o 
	gcc -g gaCable.o memory.o -o gaCable -lgsl -lgslcblas -lm
gaCable.o: gaCable.c memory.h 
	gcc -c -g gaCable.c