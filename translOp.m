function [Tc1c2F] = translOp(F, sampling, c1, c2, kappa)

%% TRANSLOP: Applies generalized translation operator to a far field matrix F.
%
% INPUT:    F           Far field matrix, nxhat*nd-array.
%           sampling    Structure containing information about the discretization.    
%           c1          Determines translation in row direction, vector of length 2.
%           c2          Determines translation in column direction, vector of length 2.
%           kappa       Wave number, >0.
%
% OUTPUT:   Tc1c2F  Translated far field matrix, nxhat*nd-array.  
%
% SYNTAX:   translOp(F, sampling, c1, c2, kappa)
%
% **************************************************************************************

xhat = sampling.xhat;
d = sampling.d.';

Tc1 = diag(exp(1i*kappa*[cos(xhat), sin(xhat)]*c1));
Tc2 = diag(exp(1i*kappa*[cos(d), sin(d)]*c2));
Tc1c2F = Tc1*F*Tc2';

end
