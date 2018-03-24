function info = catsdp(F)
  
% Author Johan Löfberg 
% $Id: lmiinfo.m,v 1.3 2005/02/04 10:10:27 johanl Exp $

info.sdp = [];
info.lin = [];
info.equ = [];
info.soc = [];
info.rlc = [];

Counter = size(F.clauses,2);
for i = 1:Counter
    Fi = F.clauses{i}.data;
    switch  F.clauses{i}.type;
        case {1,9}
            info.sdp = [info.sdp;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 2
            info.lin = [info.lin;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 3
            info.equ = [info.equ;size(Fi,1) size(Fi,2) F.LMIid(i)];    
        case 4
            info.soc = [info.soc;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 5
            info.rlc = [info.rlc;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case {7,8}
        otherwise
            error('Error in lmiinfo. Please report bug');
    end
end
