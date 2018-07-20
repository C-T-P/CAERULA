objects = Main.o BASIS.o CONTRACT.o I3NSERT.o colourtools.o
path-to-gsl = /usr/include/gsl

CXX = g++ -g -std=c++11

Spectrum: $(objects)
	$(CXX) -o Spectrum -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm $(objects)

Main.o:	Main.C
	$(CXX) -c -Wall Main.C

BASIS.o: BASIS.C
	$(CXX) -c -Wall -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm BASIS.C

CONTRACT.o: CONTRACT.C
	$(CXX) -c -Wall CONTRACT.C

I3NSERT.o: I3NSERT.C
	$(CXX) -c -Wall I3NSERT.C

tensortools.o: tensortools.C
	$(CXX) -c -Wall tensortools.C

clean:
	rm -f $(objects)

force:
	touch Main.C BASIS.C CONTRACT.C I3NSERT.C colourtools.C
	make
