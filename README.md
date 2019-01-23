# Spectrum
Spectrum is a C++ package to calculate matrices in colour space describing the gluon exchange between two partonic legs needed for the construction of the soft anomalous dimension matrix. 

## COMPILING
Just run 'make' in the Spectrum directory.

## USAGE
Spectrum can 

* evaluate ('-e') and simplify ('-s') strings of colour amplitudes with fully or partially contracted indices
* compute the soft matrix and all colour change matrices from files ('-f'), see /examples for reference input files
* compute the soft matrix and all colour change matrices for automatically constructed trace (default) and adjoint basis ('-adj') by specifying the number of quark pairs ('-nqp') and gluons ('-ng') in the process

See './Spectrum -help' for further informations.

## AUTHOR
Christian Tobias Preuss

## LICENSE
The Spectrum package is licensed under the GNU GPL version 3 - see the [COPYING.md](COPYING.md) for details.
