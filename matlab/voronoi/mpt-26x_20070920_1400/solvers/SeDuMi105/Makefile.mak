###########################################################################
# Makefile for SeDuMi, makes DLL-files out of C-files.
# For DOS (i.e. NT/win95/win98/2000/XP/etc.). See "Makefile" for UNIX version.
# Use  "make -f Makefile.mak" (GNU MINGW, Borland, lcc, etc.) or
# "nmake -f Makefile.mak" (Microsoft)
# Uninstall by issuing "del *.dll" from the command prompt.
#
# NOTE: at the DOS prompt (or Control Panel/System/Advanced), type
#   set MATLABROOT=c:\matlab
#   set MATLAB=%MATLABROOT%
# before invoking nmake. The above assumes MATLAB is installed under
# c:\matlab. Another popular place is:
#   set MATLABROOT="C:\Windows\Application Data\Mathworks\Matlab"
#
# See the file `Install.dos' for details.
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

##### Instead of setting the MATLABROOT environment variable,
##### you may define it below: (by uncommenting)
#MATLABROOT = C:\MATLABR11

# If your system does not have PERL, then use "PERL=call".
#PERL=perl
PERL=call
MEX=$(PERL) $(MATLABROOT)\bin\mex.bat
# Compiler flags for final version. Make empty for debugging version.
MEXFLAGS=-O -DNDEBUG -I$(MATLABROOT)\extern\include
#MEXFLAGS=-I$(MATLABROOT)\extern\include
#MEXFLAGS=

# standard options are provided in this file, which is sourced prior
# to every mex command.  You may want to add your own compiler specific
# options to this file, or rewrite the file itself to suit your needs.
# WARNING: users of certain MATLAB versions should uncomment some of the
# definitions below. Usually, mex finds the default location by itself, though.
MEXOPTS=
#MEXOPTS=-f $(MATLABROOT)\bin\mexopts.bat
#MEXOPTS=-f $(MATLABROOT)\mexopts.bat

# Extension of DLL-executable libraries is dll
ST=dll

#
# C to MEX compilation rule
#
.c.dll:
	$(MEX) $(MEXOPTS) $(MEXFLAGS) $*.c

# To uninstall the mex-files, type "del *.dll"

# Having made the DOS-specific macro definitions, proceed with the
# common makefile for SeDuMi.
# WARNING: users of the lcc compiler Make (delivered with Matlab R11) must
# include "Makefile.sedumi" physically (using a text editor).
include Makefile.sedumi
