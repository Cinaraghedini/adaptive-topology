function [Q,c,f,x,info] = vecquaddecomp(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

% Author Johan Löfberg
% $Id: vecquaddecomp.m,v 1.2 2006/08/10 14:09:06 joloef Exp $

[n,m]=size(p);
info = 0;

Q = cell(n*m,1);
c = cell(n*m,1);
f = cell(n*m,1);

% Involved in polynomial expression
[mt,variabletype] = yalmip('monomtable');
x_lin = getvariables(p);
x_var = find(any(mt(x_lin,:),1));
if nargin==2
    x_var = union(x_var,depends(z));
end
x = recover(x_var);


basis = getbase(p);
vars = getvariables(p);
mm = mt;
for index = 1:n*m
    mt = mm;
    base = basis(index,:);
    %    all(variabletype(vars(find(base(2:end)))) == 0)
    if all(variabletype(vars(find(base(2:end)))) == 0)% is(p(index),'linear')
        n = length(x);
        Q{index} = spalloc(n,n,0);
        fc = basis(index,:);
        f{index} = fc(1);
        c{index} = zeros(length(x),1);
        for i = 1:length(vars)
            c{index}(find(vars(i)==x_var)) = fc(1+i);
        end
    else
        % degrees = sum(mt(x_lin,:),2);
        if all(variabletype(vars(find(base(2:end)))) <= 2)%all(degrees<=2)
            Qtemp = spalloc(length(x),length(x),2*nnz(base));
            ctemp = zeros(length(x),1);

            if nnz(base(1))==0
                f{index} = 0;
                base = base(2:end);
            else
                f{index} = base(1);
                base = base(2:end);
            end

            mt = mt(x_lin,x_var);

            [jj,ii,kk] = find(mt');ii = [ii(:) ;0];
            top = 1;
            for i = 1:length(x_lin)
                if ii(top) == ii(top+1)
                    j = [jj(top) jj(top+1)];
                    top = top + 2;
                else
                    j = [jj(top)];
                    top = top + 1;
                end              
                if length(j)==1    % One variable
                    if kk(top-1)==1
                        ctemp(j)=base(i);
                    else
                        Qtemp(j,j)=Qtemp(j,j) + base(i)/2;
                    end
                else
                    Qtemp(j(1),j(2))=Qtemp(j(1),j(2)) + base(i)/2;
                end
            end
            Q{index} = Qtemp+Qtemp';
            c{index} = ctemp;
        else
            info = 1;
            x = [];
        end
    end
end

