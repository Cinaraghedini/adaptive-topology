% this subroutine was written by Michael Kleder and was found on The
% Mathwork's MATLAB Central file exchange.
 
function B = allspath(A)
B=full(double(A));
B(B==0)=Inf;
C=ones(size(B));
iter=0;
while any(C(:))
    C=B;
    B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
    repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
    C=B-C;
end
B(logical(eye(length(B))))=0;
return