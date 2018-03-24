function regress_weird2

clear sin % To avoid aby problems from code above
sdpvar x
strange = @(x) sdpfun(x,'@(x) sin(10*x)+abs(sin(x))+x');
sol = solvesdp(set(-pi < x < pi),strange(x),sdpsettings('solver','bmibnb'));
clear sin

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),-2.68320796393170, 1e-2);