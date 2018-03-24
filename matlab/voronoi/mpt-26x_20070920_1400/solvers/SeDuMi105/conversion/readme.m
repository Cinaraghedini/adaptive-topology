% Conversion to SeDuMi.
%   Conversion routines for reading in problems from SDPPack, SDPA
%   or MPS formats. Some routines require SDPPACK and/or LIPSOL to
%   be installed as MATLAB Toolboxes, or need SDPPACK and SDPLIB
%   environment variables.
%
%   Copyright (C) 1998 Jos F. Sturm
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
% 
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
% Conversion to SeDuMi
%   feascpx     -  Random feasible problem with complex entries
%   feasreal    -  Random feasible problem with real entries
%   frompack    -  SDPPACK to SeDuMi conversion
%   fromsdpa    -  Read problem stored in SDPA format, such as SDPLIB set.
%   getproblem  -  Read in problem from MPS, SDPA or SDPPack.
%   prelp       -  Read problem stored in MPS file, such as NETLIB set.
%   runsp       -  Execute SP on an SDP problem in (At,b,c,K) format.
%
