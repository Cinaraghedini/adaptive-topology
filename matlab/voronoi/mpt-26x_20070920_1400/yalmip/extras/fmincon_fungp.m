function [f,df] = fmincon_fungp(x,dummy1,dummy2,dummy4)

persistent prob
persistent A
persistent b

if nargin == 4
   prob = dummy4;
   ind = find(prob.map==0);
   A = prob.A(ind,:);
   b = prob.b(ind);
   return
end

z = A*x;
w = exp(z);
f = full(log(b'*w));
df = ((1/(sum(b.*w)))*(b.*w)'*A);