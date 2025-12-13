function [M] = elementwise_soft_thresholding(M, epsilon)

%% Applies elementwise soft thresholding to matrix M
%
% INPUT:    M - Matrix with complex entries
%           epsilon - Threshold parameter
%
% OUTPUT:   M - M after elementwise thresholding has been applied
%
% *************************************************************************

absM = abs(M);
phaseM = M ./ absM;

Index = (absM < epsilon/2);

absM = absM - epsilon/2;
absM(Index) = 0;

M = absM .* phaseM;

M(Index) = 0;

end

