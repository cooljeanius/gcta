PROGRAM: GCTA

DESCRIPTION: Genome-wide Complex Trait Analysis

AUTHOR: Jian Yang, Hong Lee, Mike Goddard and Peter Visscher

CONTACT: jian.yang@uq.edu.au

YEAR: 2009-2011

LICENSE: Released under GNU General Public License, v2 (see COPYING.txt)

DOCUMENTATION: http://gump.qimr.edu.au/gcta/

LIBRARY: require the library EIGEN (version 3.0.3) support (see http://eigen.tuxfamily.org/). You will need to specify the path of library EIGEN in the Makefile (keyword EigenLib).


COMPILATION: You will need a standard C/C++ compiler such as GNU gcc and type

make -f Makefile

(at least, that was the original compilation process; this fork is an attempt
to autotools-ize the build process...)

OTHER NOTES: If you are viewing this on GitHub, these are not the original GCTA
sources. This fork is just a random person messing around with it in the process
of attempting to package it for a package manager.
