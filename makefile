o_objects = Main.o Insert.o Colourtools.o CBasis.o TraceBasis.o FBasis.o GenBasis.o MultipletBasis.o
c_objects = Main.C Insert.C Colourtools.C CBasis.C TraceBasis.C FBasis.C GenBasis.C MultipletBasis.C

CXX = g++ -Wall -std=c++11

Spectrum: $(o_objects)
	$(CXX) $(CXXFLAGS) -o Spectrum -lm $(o_objects)

Main.o:	Main.C
	$(CXX) -c -Wall Main.C

Insert.o: Insert.C
	$(CXX) -c -Wall Insert.C

Colourtools.o: Colourtools.C
	$(CXX) -c -Wall Colourtools.C

CBasis.o: CBasis.C
	$(CXX) -c -Wall CBasis.C

TraceBasis.o: TraceBasis.C
	$(CXX) -c -Wall TraceBasis.C

FBasis.o: FBasis.C
	$(CXX) -c -Wall FBasis.C

GenBasis.o: GenBasis.C
	$(CXX) -c -Wall GenBasis.C

MultipletBasis.o: MultipletBasis.C
	$(CXX) -c -Wall MultipletBasis.C

clean:
	rm -f $(o_objects)

force:
	touch $(c_objects)
	make
