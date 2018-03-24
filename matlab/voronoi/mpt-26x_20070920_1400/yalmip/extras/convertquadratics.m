function [Fconv,no_changed,infeasible] = convertquadratics(F)
%CONVERTQUADRATICS Internal function to extract quadratic constraints

% Author Johan L�fberg
% $Id: convertquadratics.m,v 1.13 2007/07/28 10:04:54 joloef Exp $

% ******************************
% LINEAR?
% ******************************

infeasible = 0;
itslinear = islinear(F);
if itslinear
    Fconv = F;
    no_changed = 0;
    return
end

itssigmonial = issigmonial(F);
if itssigmonial
    Fconv = F;
    no_changed = 0;
    return
end

Fconv = lmi;
no_changed = 0;
for i = 1:1:length(F)
    if is(F(i),'element-wise') & ~is(F(i),'linear') & ~is(F(i),'sigmonial')
        % f-c'*x-x'*Q*x>0
        fi = sdpvar(F(i));fi = fi(:);
        for j = 1:length(fi)
            fij = fi(j);
            if isa(fij,'double')
                if fij < 0
                    infeasible = 1;                    
                    warning('The problem is trivially infeasible. You have a term POSITIVE NUMBER < 0.');
                    return
                end
            end
            [Q,c,f,x,info] = quaddecomp(fij);
            if info==0
                if nnz(Q)==0
                    % Oh, linear,...
                    Fconv = Fconv + set(fi(j));
                else
                    % Yes, quadratic, but convex?
                    % Change sign definitions
                    Q = -Q;
                    c = -c;
                    f = -f;
                    % Semi-definite case when only part of x in Q
                    % Occurs, e.g, in constraints like y'*Q*y < t
                    used = find(any(Q));Qred=Q(:,used);Qred = Qred(used,:);xred = x(used);
                    [R,p]=chol(Qred);
                    if p
                        % Safety check to account for low rank problems
                        if all(eig(Qred)>=0)
                            [u,s,v]=svd(full(Qred));
                            r=find(diag(s)>1e-12);
                            R=(u(:,r)*sqrt(s(r,r)))';
                            p=0;
                        end
                    end
                    if p==0
                        % Write as second order cone
                        d = -c'*x-f;
                        if isa(d,'double') & d==1
                            Fconv=Fconv + lmi(cone([2*R*xred],1+d));
                        else
                            Fconv=Fconv + lmi(cone([2*R*xred;1-d],1+d));
                        end
                        no_changed = no_changed + 1;
                    else
                        Fconv = Fconv + set(fi(j));
                    end
                end
            else
                Fconv = Fconv + set(fi(j));
            end
        end
    else
        Fconv = Fconv + F(i);
    end
end
