function [g,geq,dg,dgeq] = fmincon_congp(x,dummy1,dummy2,dummy4)

persistent prob
persistent B

if nargin == 4
    prob = dummy4;prob.b = full(prob.b);
    %B = spalloc(max(prob.map),size(prob.A,1),0); 
    [aa,bb,cc] = unique(prob.map);
    if aa(1)==1
        bb=[0 bb(1:end)'];
    end   
    indsi = [];
    indsj = [];
    vals = [];
    for i = 1:max(prob.map)  
        ind = [(bb(i)+1):(bb(i+1))];
        indsj = [indsj ind];
      %  indsi = [indsi repmat(i,1,length(ind))];
      %  vals = [vals prob.b(ind)'];        
      %  ind = find(prob.map==i);
      %  B(i,ind) = prob.b(ind)';
    end 
    indsi = prob.map(prob.map>0);
    vals = prob.b(prob.map>0);
    B = sparse(indsi,indsj,vals,max(prob.map),size(prob.A,1));
    
    return
end

z = prob.A*x;
expz = exp(z);
g = B*expz;

% ind = find(g<1e-2);
% if ~isempty(ind)
%     g(ind) = exp(log(1e-2)+(g(ind)-1e-2)/1e-2);
% end
g = log(g);

if  length(prob.h) > 0
    geq = log(prob.h) + prob.G*x;
    dgeq = prob.G';
else
    geq = [];
    dgeq = [];
end
% Should be correct, but it fails for some problems (test_gp_5)
dg = [];
%z = prob.A*x;
%expz = exp(z);
expz(isinf(expz)) = 1e5;
T=diag(sparse(1./(B*expz)));
U=diag(sparse(expz));
dg = (T*B*(U*prob.A))';