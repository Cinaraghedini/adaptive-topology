function [nsv, alpha, b0] = svc(X,Y,ker,C,solver)
%  SVC Support Vector Classification
%
%  Usage: [nsv alpha bias] = svc(X,Y,ker,C)
%
%  Parameters: X      - Training inputs
%              Y      - Training targets
%              ker    - kernel function
%              C      - upper bound (non-separable case)
%              nsv    - number of support vectors
%              alpha  - Lagrange Multipliers
%              b0     - bias term
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk), modified by Giancarlo
%  Ferrari-Trecate for using the solver for QPs in MPT

if (nargin <2 | nargin>5) % check correct number of arguments
    help svc
else

    %fprintf('Support Vector Classification\n')
    %fprintf('_____________________________\n')
    n = size(X,1);
    if (nargin<4) C=Inf;, end
    if (nargin<3) ker='linear';, end

    % tolerance for Support Vector Detection
    epsilon = svtol(C);

    % Construct the Kernel matrix
    %fprintf('Constructing ...\n');
    H = zeros(n,n);
    for i=1:n
        for j=1:n
            H(i,j) = Y(i)*Y(j)*svkernel(ker,X(i,:),X(j,:));
        end
    end
    c = -ones(n,1);

    % Add small amount of zero order regularisation to
    % avoid problems when Hessian is badly conditioned.
    H = H+1e-6*eye(size(H));

    % Set up the parameters for the Optimisation problem

    vlb = zeros(n,1);      % Set the bounds: alphas >= 0
    vub = C*ones(n,1);     %                 alphas <= C
    x0 = zeros(n,1);       % The starting point is [0 0 0   0]
    neqcstr = nobias(ker); % Set the number of equality constraints (1 or 0)
    if neqcstr
        A = Y';, b = 0;     % Set the constraint Ax = b
    else
        A = [];, b = [];
    end

    % Solve the Optimisation Problem

    %fprintf('Optimising ...\n');
    st = cputime;

    %  this is the original line in the file by Gunn for solving svc
    %  i keep it just for comparison when trying new QP solvers
    %  [alpha lambda how] = qp(H, c, A, b, vlb, vub, x0, neqcstr);


    nconstr=size(b,1);

    %% New part added by G. Ferrari-Trecate
    %% With rescue=1, mpt tries different QP solvers, if one reports
    %% infeasibility or unboundedness
    rescue=1;
    if neqcstr==0
        [alpha,lambda,how,exitflag,fopt]=mpt_solveQP(H,c,[A;-eye(n);eye(n)],[b;vlb;vub],[],[],x0,solver,[],rescue);
        switch(lower(how))
            case{'unbounded', 'infeasible'}
                fprintf('\n =========================================================== ');
                fprintf('\n ERROR IN COMPUTING ONE HYPERPLANE SEPARATING THE REGIONS:');
                fprintf('\n mpt_solveQP reported infeasibility or unboundedness');
                fprintf('\n in solving an SVC problem');
                fprintf('\n Try to change the QP solver (modifying idpar.QPsolver) or');
                fprintf('\n try with other classification algorithms (e.g. MRLP or PSVC)');
                fprintf('\n hit_regression.m ends here. ');
                fprintf('\n ============================================================= \n');
                error('hit_regression:end','hit_regression ends here.')
        end
    elseif neqcstr==1 & nconstr==1
        [alpha,lambda,how,exitflag,fopt]=mpt_solveQP(H,c,[-eye(n);eye(n)],[vlb;vub],[A],[b],x0,solver,[],rescue);
        switch(lower(how))
            case{'unbounded', 'infeasible'}
                fprintf('\n =========================================================== ');
                fprintf('\n ERROR IN COMPUTING ONE HYPERPLANE SEPARATING THE REGIONS:');
                fprintf('\n mpt_solveQP reported infeasibility or unboundedness');
                fprintf('\n in solving an SVC problem');
                fprintf('\n Try to change the QP solver (modifying idpar.QPsolver) or');
                fprintf('\n try with other classification algorithms (e.g. MRLP or PSVC)');
                fprintf('\n ============================================================= \n');
                error('hit_regression:end','hit_regression ends here.')
        end
    else
        [alpha,lambda,how,exitflag,fopt]=mpt_solveQP(H,c,[A(2:nconstr,:);-eye(n);eye(n)],[b(2:nconstr);vlb;vub],A(1,:),b(1),x0,solver,[],rescue);
        switch(lower(how))
            case{'unbounded', 'infeasible'}
                fprintf('\n =========================================================== ');
                fprintf('\n ERROR IN COMPUTING ONE HYPERPLANE SEPARATING THE REGIONS:');
                fprintf('\n mpt_solveQP reported infeasibility or unboundedness');
                fprintf('\n in solving an SVC problem');
                fprintf('\n Try to change the QP solver (modifying idpar.QPsolver) or');
                fprintf('\n try with other classification algorithms (e.g. MRLP or PSVC)');
                fprintf('\n ============================================================= \n');
                error('hit_regression:end','hit_regression ends here.')
        end
    end

    %% End of New part

    %fprintf('Execution time: %4.1f seconds\n',cputime - st);
    %fprintf('Status : %s\n',how);
    %w2 = alpha'*H*alpha;
    %fprintf('|w0|^2    : %f\n',w2);
    %fprintf('Margin    : %f\n',2/sqrt(w2));
    %fprintf('Sum alpha : %f\n',sum(alpha));


    % Compute the number of Support Vectors
    svi = find( alpha > epsilon);
    nsv = length(svi);
    %fprintf('Support Vectors : %d (%3.1f%%)\n',nsv,100*nsv/n);

    % Implicit bias, b0
    b0 = 0;

    % Explicit bias, b0
    if nobias(ker) ~= 0
        % find b0 from average of support vectors on margin
        % SVs on margin have alphas: 0 < alpha < C
        svii = find( alpha > epsilon & alpha < (C - epsilon));
        if length(svii) > 0
            b0 =  (1/length(svii))*sum(Y(svii) - H(svii,svi)*alpha(svi).*Y(svii));
        else
            fprintf('No support vectors on margin - cannot compute bias.\n');
        end
    end

end

