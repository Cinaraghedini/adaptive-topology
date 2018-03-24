###########################################################################
# Makefile for Sedumi 1.05, makes MEX-files out of C-files.
# For UNIX (e.g. Linux/sun/hp900/dec-alpha/ibm-rs).
# See "Makefile.mak" for DOS version.
# Use "make". Uninstall with "make clean".
#
# See the file `Install.unix' for details.
#
# IMPORTANT: you should have the variable MATLABROOT (below) point
#   to the directory where Matlab is located, AND define the variable
#   EXARCH corresponding to your computer architecture.
#   You may also need to specify a MEXOPTS file w.r.t your C-compiler.
#
# This file is part of SeDuMi 1.05
# Copyright (C) 2001 Jos F. Sturm
#   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
#   Supported by the Netherlands Organization for Scientific Research (NWO).
# Affiliation SeDuMi 1.03 and 1.04Beta (2000):
#   Dept. Quantitative Economics, Maastricht University, the Netherlands.
# Affiliations up to SeDuMi 1.02 (AUG1998):
#   CRL, McMaster University, Canada.
#   Supported by the Netherlands Organization for Scientific Research (NWO).
#
###########################################################################

# The USER should provide the root directory of the MATLAB installation:
MATLABROOT=/usr/local/matlab5

# Determine extension of mex-binaries. E.g., if the extension is ".mexsol",
# then the USER should set EXARCH=sol   (i.e., without .mex).
# Uncomment the EXARCH= line correspomnding to your system.

EXARCH=sol                        # Sun (sun4, sol2)
#EXARCH=lx                         # Linux (i86, PII, etc)
#EXARCH=axp                        # DEC Alpha
#EXARCH=hp7                        # hp700
#EXARCH=rs6                        # ibm_rs
#EXARCH=sg                         # sgi
#EXARCH=sg64                       # sgi64
#EXARCH=4                          # Outdated (Sun)

# use MATLABROOT variable to find location of MEX compiler.
MEX=$(MATLABROOT)/bin/mex

# Compiler flags for final version. Make empty for debugging version.
MEXFLAGS=-O -DNDEBUG
#MEXFLAGS=

# standard options are provided in this file, which is sourced prior
# to every mex command.  You may want to add your own compiler specific
# options to this file, or rewrite the file itself to suit your needs.
MEXOPTS=
#MEXOPTS=-f my_mexopts.sh
#MEXOPTS=-f $(MATLABROOT)/bin/mexopts.sh

# Extension of MEX-executables is mex+exarch (e.g. mexsol).
ST=mex$(EXARCH)

#
# C to MEX compilation rule
# Currently not needed: use explicit rules. (Hence outcommented)
#
#%.$(ST):%.c
#	$(MEX) $(MEXOPTS) $(MEXFLAGS) $^

# Having made the UNIX-specific macro definitions, proceed with the
# common makefile for SeDuMi.
include Makefile.sedumi
