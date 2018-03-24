function [f,df] = fmincon_fun(x,model)

xevaled = zeros(1,length(model.c));
xevaled(model.linearindicies) = x;

% Apply the precomputed evaluation scheme (if necessary)
if ~model.SimpleLinearObjective
    xevaled = apply_recursive_evaluation(model,xevaled);
end

xevaled = xevaled(:);
if model.SimpleLinearObjective
    f = model.f + model.c'*xevaled;
else
    f = model.f + (model.c'+xevaled'*model.Q)*xevaled;
end

if model.SimpleLinearObjective
    df = model.c(model.linearindicies);
elseif model.SimpleQuadraticObjective
    df = model.c(model.linearindicies) + 2*model.Q(model.linearindicies,model.linearindicies)*x;
elseif model.SimpleNonlinearObjective
    requested = model.c | any(model.Q,2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    df = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    X = repmat(xevaled,size(mtNonlinear,1),1); 
    r = find(mtNonlinear );
    X = X(r);
    Xones = ones(size(mtNonlinear,1),size(mtNonlinear,2));
    for i = 1:length(linearindicies)
        if requested(i)
           % xevaled = zeros(1,n);
           % xevaled(linearindicies) = x;
            mt = mtNonlinear;
            oldpower = mtNonlinear(:,linearindicies(i));
            mt(:,linearindicies(i)) = mt(:,linearindicies(i))-1;
%            r = find(mtNonlinear );
%            Z = X(r).^mt(r);
            Z = X.^mt(r);
            XX = Xones;%ones(size(mt,1),size(X,2));
            XX(r) = Z;
            xevaledNonLinear = prod(XX,2);
  %         xevaledNonLinear = prod(X.^mt,2);            
            xevaledNonLinear = xevaledNonLinear(:)'.*oldpower';xevaledNonLinear(isnan(xevaledNonLinear))=0;
            dx = zeros(1,n);
            dx(linearindicies(i)) = 1;
            dx(model.nonlinearindicies) = xevaledNonLinear;
            df = [df;model.c'*dx'];
        else
            df = [df;zeros(1,length(n))];
        end
    end
    df = real(df + 2*model.Q(model.linearindicies,model.linearindicies)*x);
elseif nargout > 1
    requested = model.c | any(model.Q,2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    dx = apply_recursive_differentiation(model,xevaled,requested);
    df = model.c'*dx+2*xevaled'*model.Q*dx;
else
    df = [];
end






