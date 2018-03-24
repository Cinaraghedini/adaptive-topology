function J = ipopt_congrad(x,params)

J = [];

xevaled = zeros(1,length(params.interfacedata.c));

for i = 1:length(x)    
    xevaled(params.linearindicies) = x;
    xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^params.monomtable(params.nonlinearindicies,:),2);
    

xevaled = xevaled(:);
g = params.interfacedata.F_struc*[0;xevaled];
%;f = params.interfacedata.c'*xevaled+xevaled'*params.interfacedata.Q*xevaled;
