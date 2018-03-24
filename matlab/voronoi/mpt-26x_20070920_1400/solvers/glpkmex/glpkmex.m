% GLPKMEX.M MEX Interface for the GNU GLPK library
%     Copyright (C) 2001-2005  Nicolo' Giorgetti
%
%
%  This file is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
% This part of code is distributed with the FUTHER condition that it is 
% possible to link it to the Matlab libraries and/or use it inside the Matlab 
% environment.
%
%  This file is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
%  License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this file; see the file COPYING. If not, write to the Free
%  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
%  02111-1307, USA.
%
%    Nicolo' Giorgetti (giorgetti@dii.unisi.it)
%    DII, University of Siena, Italy
%
%
%
% This routine calls the glpk library to solve a LP/MIP problem. A typical
% LP problem has following structure:
%
%            [min|max] C'x
%             s.t.
%                  Ax ["="|"<="|">="] b
%                  {x <= UB}
%                  {x >= LB}
%
% The calling syntax is:
% [XMIN,FMIN,STATUS,EXTRA]=glpkmex(SENSE,C,...
%                                       A,B,CTYPE,LB,UB,...
%                                       VARTYPE,PARAM,LPSOLVER,SAVE)
%
% SENSE:     indicates whether the problem is a minimization
%            or maximization problem.
%            SENSE = 1 minimize
%            SENSE = -1 maximize.
%
% C:         A column array containing the objective function
%            coefficients.
%
% A:         A matrix containing the constraints coefficients. A
%            may be a sparse matrix (see 'help sparse' in matlab,
%            for details).
%
% B:         A column array containing the right-hand side value for
%            each constraint in the constraint matrix.
%
% CTYPE      A column array containing the sense of each constraint
%            in the constraint matrix.
%            CTYPE(i) = 'F'  Free (unbounded) variable
%            CTYPE(i) = 'U'  "<=" Variable with upper bound
%            CTYPE(i) = 'S'  "="  Fixed Variable
%            CTYPE(i) = 'L'  ">=" Variable with lower bound
%            CTYPE(i) = 'D'  Double-bounded variable
%            (This is case sensitive).
%
% LB         An array of at least length numcols containing the lower
%            bound on each of the variables.
%
% UB         An array of at least length numcols containing the upper
%            bound on each of the variables.
%
% VARTYPE    A column array containing the types of the variables.
%            VARTYPE(i) = 'C' continuous variable
%            VARTYPE(i) = 'I' Integer variable
%            (This is case sensitive).
%
% PARAM      A structure containing some parameters used to define
%            the behavior of solver. For more details type
%            HELP GLPKPARAMS.
%
% LPSOLVER   Selects which solver using to solve LP problems.
%             LPSOLVER=1  Revised Simplex Method
%             LPSOLVER=2  Interior Point Method
%            If the problem is a MIP problem this flag will be ignored.
%
% SAVE       Saves a copy of the problem if SAVE<>0.
%            The file name can not be specified and defaults to "outpb.lp".
%            The output file is CPLEX LP format.
%
%
% XMIN       The optimizer.
%
% FMIN       The optimum.
%
% STATUS     Status of the optimization.
%      
%               - Simplex Method -
%               Value   Code     
%               180     LPX_OPT     solution is optimal   
%               181     LPX_FEAS    solution is feasible
%               182     LPX_INFEAS  solution is infeasible
%               183     LPX_NOFEAS  problem has no feasible solution
%               184     LPX_UNBND   problem has no unbounded solution
%               185     LPX_UNDEF   solution status is undefined
%
%               - Interior Point Method -
%               Value   Code
%               150     LPX_T_UNDEF the interior point method is undefined
%               151     LPX_T_OPT   the interior point method is optimal
%               *  Note that additional status codes may appear in
%               the future versions of this routine *
%     
%               - Mixed Integer Method -
%               Value   Code
%               170     LPX_I_UNDEF  the status is undefined
%               171     LPX_I_OPT    the solution is integer optimal
%               172     LPX_I_FEAS   solution integer feasible but
%                                    its optimality has not been proven
%               173     LPX_I_NOFEAS no integer feasible solution
%               
% EXTRA      A data structure containing the following fields:
%              LAMBDA     Dual variables
%              REDCOSTS   Reduced Costs
%              TIME       Time (in seconds) used for solving LP/MIP problem in seconds.
%              MEM        Memory (in bytes) used for solving LP/MIP problem.
%
%
% In case of error the glpkmex returns one of the
% following codes (these codes are in STATUS). For more informations on
% the causes of these codes refer to the GLPK reference manual.
%
%  Value  Code
%  204    LPX_E_FAULT  unable to start the search
%  205    LPX_E_OBJLL  objective function lower limit reached
%  206    LPX_E_OBJUL  objective function upper limit reached
%  207    LPX_E_ITLIM  iterations limit exhausted
%  208    LPX_E_TMLIM  time limit exhausted
%  209    LPX_E_NOFEAS no feasible solution
%  210    LPX_E_INSTAB numerical instability
%  211    LPX_E_SING   problems with basis matrix
%  212    LPX_E_NOCONV no convergence (interior)
%  213    LPX_E_NOPFS  no primal feas. sol. (LP presolver)
%  214    LPX_E_NODFS  no dual feas. sol.   (LP presolver)
%