o_objects = Main.o Spectrum.o Insert.o Colourtools.o CBasis.o TraceBasis.o FBasis.o GenBasis.o MultipletBasis.o
c_objects = Main.C Spectrum.C Insert.C Colourtools.C CBasis.C TraceBasis.C FBasis.C GenBasis.C MultipletBasis.C

CXX = g++ -Wall -std=c++11

Spectrum: $(o_objects)
	$(CXX) $(CXXFLAGS) -o Spectrum -lm $(o_objects)

Main.o:	Main.C
	$(CXX) -c Main.C

Spectrum.o: Spectrum.C
	$(CXX) -c Spectrum.C 

Insert.o: Insert.C
	$(CXX) -c Insert.C

Colourtools.o: Colourtools.C
	$(CXX) -c Colourtools.C

CBasis.o: CBasis.C
	$(CXX) -c CBasis.C

TraceBasis.o: TraceBasis.C
	$(CXX) -c TraceBasis.C

FBasis.o: FBasis.C
	$(CXX) -c FBasis.C

GenBasis.o: GenBasis.C
	$(CXX) -c GenBasis.C

MultipletBasis.o: MultipletBasis.C
	$(CXX) -c MultipletBasis.C

clean:
	rm -f $(o_objects)

force:
	touch $(c_objects)
	make
