% GLPKMEX MEX Interface for the GLPK Callable Library
%
% Copyright (C) 2004-2005, Nicolo' Giorgetti, All right reserved.
% E-mail: <giorgetti@dii.unisi.it>.
% 
% This file is part of GLPK.
% 
% GLPKMEX is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This part of code is distributed with the FUTHER condition that it is 
% possible to link it to the Matlab libraries and/or use it inside the Matlab 
% environment.
%
% GLPKMEX is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
% License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GLPKMEX; see the file COPYING. If not, write to the Free
% Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
% 02111-1307, USA.

disp('IMPORTANT !!!!! You have to update GLPK path for compiling the MEX interface !!');
pause;

filename='glpkmex.c';

%%%%%% IMPORTANT %%%%%%%%%%
% Update this path
GLPKpath='c:\cygwin\home\giorgetti\software\glpk-4.7';
%%%%%%%%%%%%%%%%%%%%%%%%%%

include=[GLPKpath '\include'];
lib=[GLPKpath '\src'];
library=[lib '\libglpk.a'];


cygwin='-f mexopts.bat ';

cmd=['-I' include ' ' filename ' ' library];
%eval(['mex ' cmd]);
eval(['mex ' cygwin cmd]);
