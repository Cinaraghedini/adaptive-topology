function YESNO = is(F,property)
%IS   Check property of constraint.
%   d = IS(x,property) returns 1 if 'property' holds
%
%   Properties possible to test are: 'elementwise', 'sdp', 
%   'socc', 'equality', 'lmi', 'linear', 'kyp', 'sos'

% Author Johan Löfberg
% $Id: is.m,v 1.17 2006/12/14 15:37:52 joloef Exp $

if isempty(F.clauses)
    YESNO = 0;
else
    %   for i = 1:length(F.clauses)
    %       Fi = F.clauses{i};
    YESNO=zeros(length(F.clauses),1);
    switch property
        case 'equality'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==3);
            end
        case {'element-wise','elementwise'}
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==2);
            end
        case {'socc','socp'}
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==4);
            end
        case 'sdp'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type==1;
            end
        case 'lmi'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = full((((Fi.type==1) | (Fi.type==9)) | ((Fi.type==2) & (prod(size(Fi.data))==1))) & (islinear(Fi.data)));
            end
        case 'linear'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = islinear(Fi.data);
            end
        case 'kyp'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==9);
            end
        case 'sos'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==11);
            end
        case 'eig'
            YESNO(i,1) = (Fi.type==10);
        case 'sigmonial'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                monomtable = yalmip('monomtable');
                monomtable = monomtable(getvariables(Fi.data),:);
                YESNO(i,1) = any(find(any(0>monomtable,2) | any(monomtable-fix(monomtable),2)));
            end
        case 'binary'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  8;
            end
        case 'integer'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  7;
            end
       case 'parametric'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  13;
            end            
        case 'uncertain'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  15;
            end
            
        case 'logic'
            allextvars = yalmip('extvariables');
            if isempty(allextvars)
                YESNO(i,1) = 0;
            else
                for i = 1:length(F.clauses)
                    Fi = F.clauses{i};
                    xi = getvariables(Fi.data);
                    lgc = find(ismembc(xi,allextvars));
                    if ~isempty(lgc)
                        for j = lgc
                            variable = xi(j);
                            extstruct = yalmip('extstruct',11);
                            if isequal(extstruct.fcn,'or') | isequal(extstruct.fcn,'and')
                                YESNO(i,1) = 1;
                                break
                            end
                        end
                    end
                end                
            end
        case 'lowrank'
            YESNO(i,1) = Fi.type == 14;
        case 'complex'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'complex');
            end
        case 'interval'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'interval');
            end            
        case 'real'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'real');
            end
        otherwise
            YESNO = error('Huh?');
    end
    %  end
end

YESNO = full(YESNO);