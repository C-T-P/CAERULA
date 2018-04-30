Spectrum: Main.o BASIS.o CONTRACT.o
	g++ -std=c++11 -o Spectrum Main.o BASIS.o CONTRACT.o

Main.o:	Main.C
	g++ -std=c++11 -c -Wall Main.C

BASIS.o: BASIS.C
	g++ -std=c++11 -c -Wall BASIS.C

CONTRACT.o: CONTRACT.C
	g++ -std=c++11 -c -Wall CONTRACT.C
