% filename: efficiency.m
% Purpose:  computes the efficiency of a network 
% Input: 
% - md - adjacency matrix (weighted or unweighted)
% Output: 
% - e -  network efficiency
% Reference :
% Phys Rev Lett. 2001 Nov 5;87(19):198701. Epub 2001 Oct 17.
% Efficient behavior of small-world networks.
% Latora V1, Marchiori M.

function e= efficiency(md)

 mdist = md.^(-1);
 
 mdist(~isfinite(mdist))=0;
   
e = (sum(sum(mdist)))/(length(mdist)*(length(mdist)-1));