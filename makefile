o_objects = Main.o insert.o colourtools.o c_basis.o trace_basis.o f_basis.o gen_basis.o multiplet_basis.o
c_objects = Main.C insert.C colourtools.C c_basis.C trace_basis.C f_basis.C gen_basis.C multiplet_basis.C

CXX = g++ -g -std=c++11

Spectrum: $(o_objects)
	$(CXX) $(CXXFLAGS) -o Spectrum -lm $(o_objects)

Main.o:	Main.C
	$(CXX) -c -Wall Main.C

insert.o: insert.C
	$(CXX) -c -Wall insert.C

colourtools.o: colourtools.C
	$(CXX) -c -Wall colourtools.C

c_basis.o: c_basis.C
	$(CXX) -c -Wall c_basis.C

trace_basis.o: trace_basis.C
	$(CXX) -c -Wall trace_basis.C

f_basis.o: f_basis.C
	$(CXX) -c -Wall f_basis.C

gen_basis.o: gen_basis.C
	$(CXX) -c -Wall gen_basis.C

multiplet_basis.o: multiplet_basis.C
	$(CXX) -c -Wall multiplet_basis.C

clean:
	rm -f $(o_objects)

force:
	touch $(c_objects)
	make
