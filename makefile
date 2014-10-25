TARGETS:=main

all: ${TARGETS}

main : src/probleme.o src/maillage.o main.o
	g++ -o main src/probleme.o src/maillage.o main.o

%.o: %.c
	g++ -I ./eigen -g  -c $< -o $@ 
	


clean:
	find -name '*~' -exec rm {}\;
	find -name '*.o' -exec rm {} \;

clobber: clean
	rm -rf ${TARGETS}
