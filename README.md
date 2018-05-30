To run the program call ./Spectrum path-to-run-file/run-file.txt, where in the run file (see examples folder for prototypes) the process must be specified by its incoming and outgoing legs and a basis for the process has to be given.
The program then calculates and gives out the soft matrix and colour change matrices for all possible gluon insertions between two legs.

NOTE: If the GNU scientific library is installed on a different path than /usr/include/gsl, the path variable has to be updated to the actual path to GSL in the makefile.