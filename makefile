Spectrum: Main.o BASIS.o CONTRACT.o
	g++ -o Spectrum Main.o BASIS.o CONTRACT.o

Main.o:	Main.C
	g++ -c -Wall Main.C

BASIS.o: BASIS.C
	g++ -c -Wall BASIS.C

CONTRACT.o: CONTRACT.C
	g++ -c -Wall CONTRACT.C
