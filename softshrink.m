function F = softshrink(F, softshrink_mu, softshrink_type, nr_of_sources)

%% SOFTSHRINK: Applies the soft-shrinkage operator to the Matrix F and returns the result.
%
% INPUT:    F               Cell-Array containing matrices with real/complex entries.
%           softshrink_mu   Threshold parameter in the soft shrinkage operatore.
%           softshrink_type Type of soft-shrinkage to be used ('real'/'complex').
%           nr_of_sources   Size of cell-arrays.
%
% OUTPUT:  F   Cell-Array containing matrices with real/complex entries.
%
% SYNTAX: softshrink(F, softshrink_mu, softshrink_type, nr_of_sources)
%
% ******************************************************************************************

for iters = 1:nr_of_sources

f = F{iters};

switch softshrink_type
    
    case 'real'
        
        f = real(f);
        
        ind1 = f<=-softshrink_mu/2;
        ind2 = abs(f)<softshrink_mu/2;
        ind3 = f>=softshrink_mu/2;
        
        f(ind1) = f(ind1) + softshrink_mu(ind1)/2;
        f(ind2) = 0;
        f(ind3) = f(ind3) - softshrink_mu(ind3)/2;
        
    case 'complex'
        
        absf = abs(f);
        phasef = f ./ absf;
        
        ind2 = absf<softshrink_mu/2;
        
        absf = absf - softshrink_mu/2;
        
        f = absf .* phasef;
        f(ind2) = 0;
end

F{iters} = f;

end
