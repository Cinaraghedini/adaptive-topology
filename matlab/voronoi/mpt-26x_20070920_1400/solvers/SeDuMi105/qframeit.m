% x = qframeit(lab,frmq,K)
function x = qframeit(lab,frmq,K)
 lorN = length(K.q);
 if length(lab) > 2*lorN
   lab = lab(K.l+1:K.l+2*lorN);       % Take out Lorentz spectral values
 end
 x = [(lab(1:lorN) + lab(lorN+1:end))/sqrt(2);...
       qblkmul(lab(lorN+1:end) - lab(1:lorN),frmq,K.qblkstart)];
