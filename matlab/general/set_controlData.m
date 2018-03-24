function data = set_controlData(position,param,options)

[A] = initialize_matrixA(position,param,options);

ac(1:size(A,1),1) = algebraic_connectivity_New(A,param.normalized);

eigV=compute_eigVector(A,param.normalized);

data=[position(:,1);position(:,2);ac;eigV];