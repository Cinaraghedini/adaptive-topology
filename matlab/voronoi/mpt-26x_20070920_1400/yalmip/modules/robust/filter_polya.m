function [F_xw,F_polya] = polya_constraints(F_xw,w,N)

F_polya = [];
Fvars = getvariables(F_xw);
wvars = getvariables(w);
[mt,vt] = yalmip('monomtable');
if ~(N==ceil(N)) & (N>=0)
    error('The power in robust.polya must be a non-negative integer');
end
F_new = [];
if any(sum(mt(Fvars,wvars),2)>1)
    removeF = zeros(length(F_xw),1);
    for i = 1:length(F_xw)
        Fi = sdpvar(F_xw(i));
        if length(Fi)>1 & is(Fi,'symmetric')
            % FIXME: SUUUUUPER SLOW
            P = polyapolynomial(sdpvar(Fi),w,N);
            C = [];
            V = [];
            for ii=1:length(P)
                t2 = [];
                for jj=1:length(P)
                    if isa(P(ii,jj),'double')
                        cc = cc*0;
                    else
                        [cc,vv] = (coefficients(P(ii,jj),w));
                    end
                    try
                        C = [C cc];
                        V = [V vv];
                    catch
                        error('Polya filter not yet implemented for all SDP cone cases. Please report bug')
                    end
                end
            end
            if ~isa(diff(V'),'double')
                error('Polya filter not yet implemented for all SDP cone cases. Please report bug')                
            end
            for k = 1:size(C,1)
                F_new = F_new + set(reshape(C(k,:),size(P,1),size(P,1)) > 0);
            end
            removeF(i) = 1;
        else
            Fi = Fi(:);
            removeFi = zeros(length(Fi),1);
            for k = 1:length(Fi)
                if any(sum(mt(getvariables(Fi(k)),wvars),2)>1)
                    p_polya = polyapolynomial(Fi(k),w,N);
                    F_polya = F_polya + set(coefficients(p_polya,w) > 0);
                    % this element has been taken care of
                    removeFi(k) = 1;
                end
            end
            if all(removeFi)
                % all elements removed, so we can remove the whole
                % constraint
                removeF(i) = 1;
            else
                % Keep some of the elements
                F_xw(i) = set(Fi(find(~removeFi)) > 0);
            end
        end
    end
    F_xw(find(removeF)) = [];
end
F_polya = F_polya + F_new;


function  pi_monoms = homogenize_(pi_monoms,w,Nmax,Nj)
pi_monoms = pi_monoms*sum(w)^(Nmax - Nj);

function P = polyapolynomial(p,w,N)
for i = 1:size(p,1)
    for j = 1:size(p,2)
        [pi_coeffs{i,j},pi_monoms{i,j}] = coefficients(p(i,j),w);
    end
end
Nmax = -inf;
for i = 1:size(p,1)
    for j = 1:size(p,2)
        for k = 1:length(pi_monoms{i,j})
            deg_pi_monom{i,j}(k) = degree(pi_monoms{i,j}(k));
            Nmax = max(Nmax,deg_pi_monom{i,j}(k));
        end
    end
end
for i = 1:size(p,1)
    for j = 1:size(p,2)
        if isa(pi_monoms{i,j},'sdpvar')
            for k = 1:length(pi_monoms{i,j})
                pi_monoms{i,j}(k) = homogenize_(pi_monoms{i,j}(k),w,Nmax,deg_pi_monom{i,j}(k));
            end
        end
    end
end
P = [];
sumNmax = sum(w)^(N + Nmax);
sumN    = sum(w)^(N);
for i = 1:size(p,1)
    temp = [];
    for j = 1:size(p,2)
        if isa(pi_monoms{i,j},'sdpvar')
            pij = (pi_coeffs{i,j}'*pi_monoms{i,j})*sumN;
        else
            pij = p(i,j)*sumNmax;
        end
        temp = [temp pij];
    end
    P = [P;temp];
end