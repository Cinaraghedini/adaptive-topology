function [F,feasible] = filter_boxmodel(F_xw,lower,upper,x,w)

% Creates robustified version of the uncertain set of linear inequalities
% s.t A(w)*x <= b(w) for all L <= w <= U

feasible = 1;

if length(F_xw) == 0
    F = [];
    return
end

X = sdpvar(F_xw);
b = [];
A = [];
% Some pre-calc
xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
 all_f = [];
 all_c_w = [];
% all_c_x = [];
% all_Q_xw = [];
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};
    if nonquadratic
        error('Constraints can be at most quadratic, with the linear term uncertain');
    end
    Q_ww = Q(wind,wind);
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);

     all_f = [all_f;f];
     all_c_w = [all_c_w;c_w'];
%     all_c_x = [all_c_x;-c_x'];
%     all_Q_xw = [all_Q_xw;-2*Q_xw];
%    b = [b;f + c_w'*w];
    A = [A;-(c_x + 2*Q_xw*w)'];
end
b = all_f + all_c_w*w;
% A = all_c_x + w'*all_Q_xw
% Linear uncertain constraint is (Bbetai*x + cdi) >= 0 for all w, or
% (bi' + (Bi*w)')*x + (ci'*w + di).
cd    = b;
Bbeta = -A;

F = set([]);
top = 1;

% To speed up the construction, compute the ci vectors for all constraints
% in one call ci_basis = [c1 c2 ...]
ci_basis = basis(cd',w);

for i = 1:length(b)
    cdi = cd(i);
    Bbetai = Bbeta(i,:);

    if (nnz(ci_basis(:,i))==0) & isa(Bbetai,'double')
        % This constraint row is constant
        F = F + set(Bbetai*x + cdi >= 0);
    else
        ci = ci_basis(:,i);

        di = basis(cdi,0);
        if isa(Bbetai,'double')
            Bi = zeros(1,length(w));
        else
            Bi = basis(Bbetai,w)';
        end
        bi = basis(Bbeta(i,:),0)';
        % Scale to -1,1 uncertainty
        T = diag((upper-lower))/2;
        e = (upper+lower)/2;
        if nnz(Bi) == 0
            if nnz(bi)==0
                % Basically constant + w > 0
                if  (di+e'*ci) - norm(T*ci,1) < 0
                    error('Problem is trivially infeasible');
                    feasible = 0;
                    return
                end
            else
                F = F + set(bi'*x + (di+e'*ci) - norm(T*ci,1) > 0);
            end
        else
            non_zeroBirow = find(sum(abs(Bi'),2));
            zeroBirow = find(sum(abs(Bi'),2) == 0);
            if length(non_zeroBirow)>1
                t = sdpvar(length(non_zeroBirow),1);
                F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - sum(t) >= 0) + set(-t < T(non_zeroBirow,:)*(ci+Bi'*x) < t);
            else
                F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
                F = F + set((bi'+e'*Bi')*x + (di+e'*ci) + T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
            end
        end
    end
end


function b = basis(p,w)

if isequal(w,0)
    b = getbasematrix(p,0);
else
    n = length(w);
    if  isequal(getbase(w),[zeros(n,1) eye(n)])
        b = [];
        lmi_variables = getvariables(w);
        for i = 1:length(w)
            b = [b ; getbasematrix(p,lmi_variables(i))];
        end
    else
        b = [];
        for i = 1:length(w)
            b = [b ; getbasematrix(p,getvariables(w(i)))];
        end
    end
end
b = full(b);












