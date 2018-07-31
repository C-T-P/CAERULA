o_objects = Main.o BASIS.o CONTRACT.o FCONTRACT.o I3NSERT.o colourtools.o
c_objects = Main.C BASIS.C CONTRACT.C FCONTRACT.C I3NSERT.C colourtools.C
path-to-gsl = /usr/include/gsl

CXX = g++ -g -std=c++11

Spectrum: $(o_objects)
	$(CXX) -o Spectrum -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm $(o_objects)

Main.o:	Main.C
	$(CXX) -c -Wall Main.C

BASIS.o: BASIS.C
	$(CXX) -c -Wall -L$(path-to-gsl) -I$(path-to-gsl) -lgsl -lgslcblas -lm BASIS.C

CONTRACT.o: CONTRACT.C
	$(CXX) -c -Wall CONTRACT.C

FCONTRACT.o: FCONTRACT.C
	$(CXX) -c -Wall FCONTRACT.C

I3NSERT.o: I3NSERT.C
	$(CXX) -c -Wall I3NSERT.C

tensortools.o: tensortools.C
	$(CXX) -c -Wall tensortools.C

clean:
	rm -f $(o_objects)

force:
	touch $(c_objects)
	make
