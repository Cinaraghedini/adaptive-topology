function g = ipopt_con(x,params)

%g = datasaver(4,x);

xevaled = zeros(1,length(params.interfacedata.c));
xevaled(params.linearindicies) = x;

% Experimental support for arbitrary functions
if ~isempty(params.interfacedata.evalMap)
    for i = 1:length(params.interfacedata.evalMap)
        xevaled(params.interfacedata.evalVariables(i)) = feval( params.interfacedata.evalMap{i}.fcn,xevaled(params.interfacedata.evalMap{i}.variableIndex));
    end
end

xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^params.monomtable(params.nonlinearindicies,:),2);

xevaled = xevaled(:);
g = params.interfacedata.F_struc*[1;xevaled];
%;f = params.interfacedata.c'*xevaled+xevaled'*params.interfacedata.Q*xevaled;
