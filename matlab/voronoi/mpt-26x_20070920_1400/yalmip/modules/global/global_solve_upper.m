function output = global_solve_upper(p,p_original,x,options,uppersolver)

% The bounds and relaxed solutions have been computed w.r.t to the relaxed
% bilinear model. We only need the original bounds and variables.
p.lb = p.lb(1:length(p_original.c));
p.ub = p.ub(1:length(p_original.c));
x = x(1:length(p_original.c));

p_upper = p_original;

% ...expand the current node just slightly
p_upper.lb = p.lb;
p_upper.ub = p.ub;
fixed = find(abs([p.lb-p.ub]) < 1e-5);
p_upper.lb(fixed) = (p.lb(fixed) + p.ub(fixed) )/2;
p_upper.ub(fixed) = (p.lb(fixed) + p.ub(fixed) )/2;

% Pick an initial point (this can be a bit tricky...)
% Use relaxed point, shifted towards center of box
if all(x<=p.ub) & all(x>=p.lb)
    p_upper.x0 = 0.1*x + 0.9*(p.lb+p.ub)/2;
else
    p_upper.x0 = (p.lb+p.ub)/2;
end
% Shift towards interior for variables with unbounded lower or upper
lbinfbounds = find(isinf(p.lb));
ubinfbounds = find(isinf(p.ub));
p_upper.x0(ubinfbounds) = x(ubinfbounds)+0.01;
p_upper.x0(lbinfbounds) = x(lbinfbounds)-0.01;
ublbinfbounds = find(isinf(p.lb) & isinf(p.ub));
p_upper.x0(ublbinfbounds) = x(ublbinfbounds);

change_these_lb = setdiff(1:length(p.lb),fixed);
change_these_lb = setdiff(change_these_lb,lbinfbounds);
change_these_ub = setdiff(1:length(p.lb),fixed);
change_these_ub = setdiff(change_these_ub,lbinfbounds);

p_upper.lb(change_these_lb) = 0.99*p.lb(change_these_lb)+p_original.lb(change_these_lb)*0.01;
p_upper.ub(change_these_ub) = 0.99*p.ub(change_these_ub)+p_original.ub(change_these_ub)*0.01;
p_upper.lb(isinf(p_original.lb)) = p_upper.lb(isinf(p_original.lb)) - 0.001;
p_upper.ub(isinf(p_original.ub)) = p_upper.ub(isinf(p_original.ub)) + 0.001;

p_upper.options.saveduals = 0;

ub = p_upper.ub ;
lb = p_upper.lb ;

% Remove redundant equality constraints (important for fmincon)
if p_upper.K.f > 0
    Aeq = -p_upper.F_struc(1:1:p_upper.K.f,2:end);
    beq = p_upper.F_struc(1:1:p_upper.K.f,1);
    redundant = find(((Aeq>0).*Aeq*(p_upper.ub-p_upper.lb) - (beq-Aeq*p_upper.lb) <1e-6));
    p_upper.F_struc(redundant,:) = [];
    p_upper.K.f = p_upper.K.f - length(redundant);
end

% Solve upper bounding problem
p_upper.options.usex0 = 1;

output = feval(uppersolver,p_upper);
% Project into the box (numerical issue)
output.Primal(output.Primal<p_upper.lb) = p_upper.lb(output.Primal<p_upper.lb);
output.Primal(output.Primal>p_upper.ub) = p_upper.ub(output.Primal>p_upper.ub);


