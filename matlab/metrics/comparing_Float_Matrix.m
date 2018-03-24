function equal = comparing_Float_Matrix(xy1,xy2)

% compares if two matrix are equal given specific errors.
% xy1 and xy2 matrix  

absTol = 1e-5;   %absolute error
relTol = 0.005;   %realtive error 

absError = xy1(:)-xy2(:);
relError = absError./xy1(:);
relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
equal = all( (abs(absError) < absTol) & (abs(relError) < relTol) );