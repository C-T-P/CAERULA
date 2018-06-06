objects = Main.o BASIS.o CONTRACT.o I3NSERT.o
path-to-gsl = /usr/include/gsl

Spectrum: $(objects)
	g++ -std=c++11 -o Spectrum -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm Main.o BASIS.o CONTRACT.o I3NSERT.o

Main.o:	Main.C
	g++ -std=c++11 -c -Wall Main.C

BASIS.o: BASIS.C
	g++ -std=c++11 -c -Wall -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm BASIS.C

CONTRACT.o: CONTRACT.C
	g++ -std=c++11 -c -Wall CONTRACT.C

I3NSERT.o: I3NSERT.C
	g++ -std=c++11 -c -Wall I3NSERT.C

clean:
	rm -f $(objects)

force:
	touch Main.C BASIS.C CONTRACT.C I3NSERT.C
	make
