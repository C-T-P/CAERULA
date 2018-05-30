objects = Main.o BASIS.o CONTRACT.o I3NSERT.o

Spectrum: $(objects)
	g++ -std=c++11 -o Spectrum -L/usr/include/gsl -I/usr/include/gsl -lgsl -lgslcblas -lm Main.o BASIS.o CONTRACT.o I3NSERT.o

Main.o:	Main.C
	g++ -std=c++11 -c -Wall Main.C

BASIS.o: BASIS.C
	g++ -std=c++11 -c -Wall -L/usr/include/gsl -I/usr/include/gsl -lgsl -lgslcblas -lm BASIS.C

CONTRACT.o: CONTRACT.C
	g++ -std=c++11 -c -Wall CONTRACT.C

I3NSERT.o: I3NSERT.C
	g++ -std=c++11 -c -Wall I3NSERT.C

clean:
	rm -f $(objects)
