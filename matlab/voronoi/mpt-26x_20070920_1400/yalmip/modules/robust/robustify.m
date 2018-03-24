function [F,h,failure] = robustify(F,h,ops,w)
%ROBUSTIFY  Derives robust counterpart.
%
% [Frobust,objrobust,failure] = ROBUSTIFY(F,h,options) is used to derive
% the robust counterpart of an uncertain YALMIP model.
%
%   min        h(x,w)
%   subject to
%           F(x,w) >(=) 0  for all w in W
%
% The constraints and objective have to satisfy a number of conditions for
% the robustification to be tractable. Please refer to the YALMIP Wiki for
% the current assumptions (this is constantly evolving)
% 
% Some options for the robustification strategies can be altered via the
% solver tag 'robust' in sdpsettings
%  'robust.lplp'  : Controls how linear constraints with affine
%                   parameterization in an uncertainty with polytopic
%                   description is handled. Can be either 'duality' or
%                   'enumeration'
%  'polya'        : Controls the relaxation order of polynomials. If set to
%                   NAN, the polynomials will be eliminated by forcing the 
%                   coefficients to zero
%
% See also UNCERTAIN

% Author Johan Löfberg
% $Id: robustify.m,v 1.39 2007/08/29 13:17:04 joloef Exp $
if nargin < 3
    ops = [];
end

if nargin < 4
    w = [];
end

if isempty(w)
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = recover(getvariables(sdpvar(F(find(unc_declarations)))));
        F = F(find(~unc_declarations));
    else
        error('There is no uncertainty definition in the model.')
    end
end

if isempty(ops)
    ops = sdpsettings;
end

% Figure out which variables are uncertain, certain, and lifted variables
% in the uncertainty description (this code is slow and buggy as ....)
[x,w,x_variables,w_variables,aux_variables,F,failure] = robust_classify_variables(F,h,ops,w);
if failure
    return
end

% Integer variables are OK in x, but not in the uncertainty (robustification
% is based on strong duality in w-space)
integervars = [yalmip('binvariables') yalmip('intvariables')];
ind = find(is(F,'integer') | is(F,'binary'));
if ~isempty(ind)
    integervars = [integervars getvariables(F(ind))];
    if any(ismember(w_variables,integervars))
        failure = 1;
        return
    end
end

% Find  uncertainty description, uncertain and certain constraints
F_w = set([]);
F_x = set([]);
F_xw = set([]);
for i = 1:length(F)
    if all(ismember(depends(F(i)),w_variables))
        % Uncertainty definition
        F_w = F_w + F(i);
    elseif all(ismember(depends(F(i)),x_variables))
        % Certain constraint
        F_x = F_x +  F(i);
    else
        % Uncertain constraint
        F_xw = F_xw + F(i);
    end
end

% Limitation in the modelling language...
if ~isempty(intersect(intersect(depends(F_xw),depends(F_w)),aux_variables))
    disp('You are most likely using a nonlinear operator to describe the');
    disp('uncertainty set (such as norm(w,1) <=1). This is currently not');
    disp('supported. Please model the constraint manually.');
    error('Uncertain model does not satisfy assumptions (nonlinear operator on uncertainty in uncertain constraint)');
end

if length(F_w)==0
    error('There is no uncertainty description in the model.');
end

% Some pre-calc
xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
% Analyze the objective and try to rewrite any uncertainty into the format
% assumed by YALMIP (
if ~isempty(h)
    
    [Q,c,f,dummy,nonquadratic] = vecquaddecomp(h,xw);
    Q = Q{1};
    c = c{1};
    f = f{1};

    if nonquadratic
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    
    Q_ww = Q(wind,wind);
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);
 
    if nnz(Q_ww) > 0
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    % Separate certain and uncertain terms, place uncertain terms in the
    % constraints instead
    if is(h,'linear')
        if isempty(intersect(getvariables(w),getvariables(h)))
            h_fixed = h;
        else
            sdpvar t
            F_xw = F_xw + set(h < t);
            h_fixed = t;
            x = [x;t];          
        end
    else
        h_fixed     = x'*Q_xx*x + c_x'*x + f;
        h_uncertain = 2*w'*Q_xw'*x + c_w'*w;
        if ~isa(h_uncertain,'double')
            sdpvar t
            F_xw = F_xw + set(h_uncertain < t);
            h_fixed = h_fixed + t;
            x = [x;t];          
        end
    end
else
    h_fixed = [];
end

% Convert quadratic constraints in uncertainty model to SOCPs. This will
% enable us to use duality based removal of uncertainties in linear
% inequalities
F_w = convertquadratics(F_w);

% Export uncertainty model to numerical format
ops.solver = '';
[aux1,aux2,aux3,Zmodel] = export(F_w,[],ops,[],[],1);

if ~isempty(Zmodel)
    if length(Zmodel.c) ~= length(w)
        error('Some uncertain variables are unconstrained.')
    end
else
    error('Failed when exporting a model of the uncertainty.')    
end

% The uncertainty model is in the complete w-space. However, it might
% happen that the uncertainty is separable in some sense. Find groups of
% uncertain variables that can be treated separately.
% FIXME: This structure is currently not exploited
% uncertaintyGroups = find_groups(Zmodel);

% Code will be added to detect uncertainty cases in a more general and
% modular way. Additional code will also be added to find hidden simple
% structures, such as norm(w,1)<1, which currently is treated as a general
% polytopic uncertainty, since the expansion hides the simplicity
% 'Box', 'Simplex', 'Conic', 'Polytopic', '2-norm', '1-norm', 'inf-norm'
% [uncertaintyTypes] = classify_uncertainty(Zmodel,uncertaintyGroups)


% Temporary conversion of bounded variables that have been modelled using
% the norm operator (remove the epigraph variable to ensure explicit
% maximization is used). This will be generalized in the next version
Zmodel = convertuncertainty(Zmodel);

% OK, we are done with the initial analysis of the involved variables, and
% check of the objective function. 

% At this point, we have to decide on the algorithm we should use for
% robustifying the constraints. There are a couple of alternatives,
% depending on uncertainty and constraints
% 1. Polya:       Polynomial uncertainty dependence, simplex uncertainty,
%                 can only be applied on LP constraints
% 2. Elimination: Last resort, tries to cancel all nonlinear uncertainties
%                 by setting coefficients to zero
% 3. Explicit:    Linear uncertainty dependence, box-model uncertainty, can
%                 only be applied on LP constraints
% 4a. Enumeration:Linear uncertainty dependence, polytopic uncertainty,
%                 arbitrary type of constraints (convex)
% 4b. Duality:    Linear uncertainty dependence, conic uncertainty, can
%                 only be applied on LP constraints

% Initialize a set of constraints that correspond to robustified
% constraints, and constraints related to Polyas theorem, and elimination
% of nonlinear terms and uncertain equalities
F_robust = set([]);

% Since nonlinearities have to be eliminated before standard approaches
% based on duality and enumeration, we start with these

% We begin by checking to see if the user wants to apply Polyas theorem.
% If that is the case, search for simplex structures, and apply Polyas.
% FIXME: Very slow and experimental
if ~isnan(ops.robust.polya)
    [aux,simplex_members] = find_simplex_models(Zmodel);
    if length(simplex_members)>0
         [F_xw, F_polya] = filter_polya(F_xw,w,ops.robust.polya);
         F_robust = F_robust + F_polya;
    end
end

% There might still be nonlinearities left in the model. These have to be
% removed. We remove all terms with w-degree larger than 1
[F_xw,F_elimination] = filter_eliminatation(F_xw,w,1);
F_robust = F_robust + F_elimination;

% Equality constraints can not be part of an uncertain problem. Any
% dependence w.r.t w in inequalities has to be removed
F_eq = F_xw(find(is(F_xw,'equality')));
F_xw = F_xw - F_eq;
[F_eq_left,F_eliminate_equality] = filter_eliminatation(F_eq,w,0);
F_robust = F_robust + F_eliminate_equality + F_eq_left;

% The problem should now be linear in the uncertainty, with no uncertain
% equality constraints. Hence, now we apply explicit maximization,
% enumeration or duality-based robustification 

% Check if the uncertainty model is a simple box model. 
[aux,lower,upper]     = find_simple_variable_bounds(Zmodel);
if all(~isinf(lower)) & all(~isinf(upper))
    F_lp = F_xw(find(is(F_xw,'elementwise')));
    F_xw = F_xw - F_lp;
    F_robust = F_robust + filter_boxmodel(F_lp,lower,upper,x,w);
end

% Pick out the uncertain linear equalities and robustify using duality.
% Duality is only option when uncertainty is conic. We use this approach
% otherwise only if user explicitly asks for it
conic = ~isequal(Zmodel.K.s,0) | ~isequal(Zmodel.K.q,0);
if conic | isequal(ops.robust.lplp,'duality')
    F_lp = F_xw(find(is(F_xw,'elementwise')));
    F_xw = F_xw - F_lp;
    F_robust = F_robust + filter_duality(F_lp,Zmodel,x,w);
end

% Robustify remaining uncertain LP/SOCP/SDP constraints and robustify by
% enumeration. 
F_conic = F_xw(find(is(F_xw,'sdp') | is(F_xw,'socc') | is(F_xw,'elementwise')));
F_xw = F_xw - F_conic;
F_robust = F_robust + filter_enumeration(F_conic,Zmodel,x,w,ops);

% If there is anything left now, it means that we do not support it (such
% as conic uncertainty in conic constraint)
if length(F_xw) > 0
    if any(~islinear(F_xw))
        error('There are some uncertain constraints that not are supported by YALMIP')
    else
        F_robust = F_robust + F_xw;
    end
end

% Return the robustfied model
F = F_robust+F_x;
h = h_fixed;

% The model has been expanded, so we have to remember this (trying to
% expand an expanded model leads to nonconvexity error)
F = expanded(F,1); % This is actually done already in expandmodel
h = expanded(h,1); % But this one has to be done manually


function groups = find_groups(model)

X = zeros(size(model.F_struc,2)-1);
top  = 1;
if model.K.f + model.K.l > 0
    A = model.F_struc(top:model.K.f+model.K.l,2:end);
    for i = 1:size(A,1)
        X(find(A(i,:)),find(A(i,:))) = 1;
    %    X(find(A(i,:)),i) = 1;
    end
    top = top + model.K.f + model.K.l;
end

if any(model.K.q) > 0
    for j = 1:length(model.K.q)
        A = model.F_struc(top:top+model.K.q(j)-1,2:end);top = top + model.K.q(j);
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

if any(model.K.s) > 0
    for j = 1:length(model.K.s)
        A = model.F_struc(top:top+model.K.s(j)^2-1,2:end);top = top + model.K.s(j)^2;
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

[a,b,c,d] = dmperm(X);
for i = 1:length(d)-1
    groups{i} = sort(a(d(i):d(i+1)-1));
end

function Zmodel = convertuncertainty(Zmodel);
% Temoporary hack, will be generalized once the framework for multiple
% uncertainty models is supported
% We are looking for k>t, -tw<t
if size(Zmodel,1) == 1+(size(Zmodel,2)-1)*2 & Zmodel.K.f==0 & Zmodel.K.l==size(Zmodel.F_struc,1)
    n = size(Zmodel.F_struc,2)-1;
    if isequal(Zmodel.F_struc(:,2:end),sparse([zeros(1,n-1) -1;[eye(n-1);-eye(n-1)] ones(2*(n-1),1)]))
        Zmodel.F_struc = [ones(2*n,1)*Zmodel.F_struc(1,1) [eye(n);-eye(n)]];
        Zmodel.K.l = 2*n;
    end
end