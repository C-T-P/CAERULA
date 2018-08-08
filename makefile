o_objects = Main.o BASIS.o CONTRACT.o FCONTRACT.o I3NSERT.o colourtools.o trace_basis.o f_basis.o
c_objects = Main.C BASIS.C CONTRACT.C FCONTRACT.C I3NSERT.C colourtools.C trace_basis.C f_basis.C

CXX = g++ -g -std=c++11

Spectrum: $(o_objects)
	$(CXX) -o Spectrum -lm $(o_objects)

Main.o:	Main.C
	$(CXX) -c -Wall Main.C

BASIS.o: BASIS.C
	$(CXX) -c -Wall -lm BASIS.C

CONTRACT.o: CONTRACT.C
	$(CXX) -c -Wall CONTRACT.C

FCONTRACT.o: FCONTRACT.C
	$(CXX) -c -Wall FCONTRACT.C

I3NSERT.o: I3NSERT.C
	$(CXX) -c -Wall I3NSERT.C

colourtools.o: colourtools.C
	$(CXX) -c -Wall colourtools.C

trace_basis.o: trace_basis.C
	$(CXX) -c -Wall trace_basis.C

f_basis.o: f_basis.C
	$(CXX) -c -Wall f_basis.C

clean:
	rm -f $(o_objects)

force:
	touch $(c_objects)
	make
