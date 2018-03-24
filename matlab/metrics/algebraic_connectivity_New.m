% filename: algebraic_connectivityNew.m
% Purpose:  computes the algebraic connectivity of a network
% Input: 
% - A - adjacency matrix 
% - normalized - 1 for normalized algebraic connectivity, otherwise 0
% Output: 
% - ac - algebraic connectivity

function ac  =  algebraic_connectivity_New(A,normalized)

D = degree(A); % degree A
ac = 0; 
if isempty(find(diag(D)==0))
    L  =  D - A;
    if normalized
        L = ((D)^(-1/2)) * L * ((D)^(-1/2));
    end
    [~,egValue]=eig(L);
    eGvalues = sort(diag(egValue));
    ac=eGvalues(2,1);    
end