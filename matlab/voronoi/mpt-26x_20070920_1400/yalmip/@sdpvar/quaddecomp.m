function [Q,c,f,x,info] = quaddecomp(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

% Author Johan Löfberg 
% $Id: quaddecomp.m,v 1.12 2007/07/26 19:55:20 joloef Exp $   

[n,m]=size(p);
info = 0;

% Is it a scalar
if (n*m==1)
    % Involved in polynomial expression
    [mt,variabletype] = yalmip('monomtable');
    x_lin = getvariables(p);
    x_var = find(any(mt(x_lin,:),1));    
    if nargin==2
        x_var = union(x_var,depends(z));
    end
    x = recover(x_var);  
    if all(variabletype(x_lin) ==0)% is(p,'linear')
        n = length(x);
        Q = spalloc(n,n,0);
        fc = getbase(p);
        f = fc(1);
        if nargin==2
            vars = getvariables(p);
            c = zeros(length(x),1);
            for i = 1:length(vars)
                c(find(vars(i)==x_var)) = fc(1+i);
            end
        else
        c = fc(2:end);c=c(:);
        end
        return
    end    
%    degrees = sum(mt(x_lin,:),2);
    variabletype = variabletype(x_lin);
    if all(variabletype<=2)%all(degrees<=2)
        %Q = spalloc(length(x),length(x),2*nnz(getbase(p)));
        %c = zeros(length(x),1);
        base = getbase(p);
        if nnz(base(1))==0
            f = 0;
            base = base(2:end);
        else
            f = base(1);
            base = base(2:end);
        end
        
        mt = mt(x_lin,x_var);

        if 1

            quads   = find (variabletype == 2);% & base);
            bilins  = find (variabletype == 1);% & base);
            consts  = find (variabletype == 0);% & base);
            
            [varsC,aux1,aux2] = find(mt(consts,:)');            
            [varsQ,aux1,aux2] = find(mt(quads,:)');
            [varsB,aux1,aux2] = find(mt(bilins,:)');

            if isempty(varsQ)
                varsQ = [];
            end
            if isempty(varsB)
                varsB = [];
            end
            if isempty(varsC)
                varsC = [];
            end
            
            c = sparse(varsC,1,base(consts),length(x),1);
            %c(varsC) = base(consts);c = c(:);
            %c1 = c;
            ii = [varsQ ; varsB(1:2:end) ; varsB(2:2:end)];
            jj = [varsQ ; varsB(2:2:end) ; varsB(1:2:end)];
            kk = [base(quads)  base(bilins)/2  base(bilins)/2];
            Q = sparse(ii,jj,kk,length(x),length(x));
            %Q1 = sparse(varsQ,varsQ,base(quads),length(x),length(x));
            %Q2 = sparse(varsB(1:2:end),varsB(2:2:end),base(bilins)/2,length(x),length(x));
            %Q3 = sparse(varsB(2:2:end),varsB(1:2:end),base(bilins)/2,length(x),length(x));
            %Q = Q1 + (Q2+Q3);

         else
            Q = spalloc(length(x),length(x),2*nnz(getbase(p)));
            c = zeros(length(x),1);
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
                        c(j)=base(i);
                    else
                        Q(j,j)=Q(j,j) + base(i)/2;
                    end
                else
                    Q(j(1),j(2))=Q(j(1),j(2)) + base(i)/2;
                end
            end
            Q = Q+Q';
        end
        
    else
        if nargout==5
            info = 1;
            Q = [];
            c = [];
            f = [];
            x = [];
        else
            error('Function is not quadratic');
        end
    end  

else
    if nargout==5
        info = 1;
        Q = [];
        c = [];
        f = [];
        x = [];
    else
        error('Function is not scalar');
    end
end
