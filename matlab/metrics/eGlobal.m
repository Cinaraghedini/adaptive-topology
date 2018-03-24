% filename: eGlobal.m
% Purpose:  computes the global efficiency of a graph
% Input: 
% - mat - adjacency matrix
% - mwI - distance matrix
% - ideal - 1 eGlobal is normalized, otherwise 0 
% Output: 
% - e -  network global efficiency
% Reference :
% Phys Rev Lett. 2001 Nov 5;87(19):198701. Epub 2001 Oct 17.
% Efficient behavior of small-world networks.
% Latora V1, Marchiori M.

function e = eGlobal(mat,mwI,ideal)

if isempty(find(mat~=0)) % there is no links in the network 
    e = 0;
    return
end

e=efficiency(mat); % computes the global efficiency

if (ideal  && e > 0)   % if is a normalized efficiency is required  
   
    e=e/efficiency(mwI) ; %efficiency is computed for the distance matrix, considering a fully connected network  
    
end
