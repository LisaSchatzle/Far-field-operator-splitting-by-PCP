function F = addnoise(F, noiselevel)

%% ADDNOISE: Adds p%=noiselevel complex valued uniformly distributed additive error to a matrix.
%
% INPUT:    F           Unperturbed far field matrix, nxhat*nd-array.
%           noiselevel  Relative noise level, between 0 and 1.
% OUTPUT:   Fnoisy  Perturbed far field matrix, nxhat*nd-array.
%
% SYNTAX: addnoise(F, noiselevel)
%
% **********************************************************************************************

if noiselevel > 0
    errMat = rand(size(F))-.5+1i*(rand(size(F))-.5);
    errMat = errMat ./ norm(errMat,'fro');
    F = F + errMat * noiselevel * norm(F,'fro');
end
