function [f,df] = ipopt_fun(x,modelIn)

persistent model

if nargin == 2
    model = modelIn;
    return
end

xevaled = zeros(1,length(model.c));
xevaled(model.linearindicies) = x;

% Apply the precomputed evaluation scheme (if necessary)
if ~model.SimpleLinearObjective
    xevaled = zeros(1,length(model.c));
    xevaled(model.linearindicies) = x;
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
    df = [];
    c = model.c;
    for i = 1:length(model.linearindicies)
        xevaled = zeros(1,length(c));
        xevaled(model.linearindicies) = x;
        mt = model.monomtable;
        oldpower = mt(:,model.linearindicies(i));
        mt(:,model.linearindicies(i)) = mt(:,model.linearindicies(i))-1;
        xevaled = prod(repmat(xevaled,size(mt,1),1).^mt,2);
        xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
        df = [df;c'*xevaled'];
    end
    df = df + 2*model.Q(model.linearindicies,model.linearindicies)*x;
else
    df = [];
end