function [F, nr_of_iterations, Snorm] = RPCP(G, X, sampling, kappa, P_OmegaC, z, lambda, mu, kmax, tol)
%% Implementation of Robust Principal Component Pursuit via accelerated
%% proximal gradient method
%
% INPUT:    G - observed possibly noisy far field operator
%           X - initial guess, (nr_of_sparse_scatterers^2+1)-cell array 
%           sampling - structure containing observation and illumination
%                      directions
%           kappa - wave number
%           P_OmegaC - indicator matrix for observable index set if data
%                      completion shall be included in considertaion;
%                      if not choose ones(nxhat,nd)
%           z - vector containing positions of all sparse components,
%               (nr_of_sparse_scatterersx2)-dimensional
%           lambda - coupling parameter 
%           mu - ragularization paramter
%           X - initial guess, (1xnr_of_sparse_scatterers^2)-cell array 
%           kmax - maximal number of RPCP iterations 
%           tol - tolerance for stopping criterion
%
% OUTPUT:   F - (nr_of_sparse_scatterers^2+1)-cell array containing the RPCP
%               approximation, F{1} low-rank component,
%               F{k*nr_of_sparse_scatterers+2} k-th sparse component
%           nr_of_iterations - number of performed RPCP iterations
%           Res_norm - vector containing RPCP residuen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nxhat = sampling.nxhat;
nd = sampling.nd;

nr_of_sparse_scatterers = size(z,2);

% generates indices-vector to get all possible combinations (j,l) for 1\leq
% j,l\leq nr_of_sparse_scatterers
% indices for components approximating sparse scatterers are at rows indices(iterkk*nr_of_sparse_scatterers+1)
% indices specifies ordering of cell-arrays
[Z,Y] = meshgrid(1:nr_of_sparse_scatterers, 1:nr_of_sparse_scatterers);
indices = [reshape(Z,nr_of_sparse_scatterers^2,1),reshape(Y,nr_of_sparse_scatterers^2,1)]; % (J^2+1)*2-array that indicates the finite dimensional subspace where the corresponding block should lie in.
clear Z Y

H = X;
Hold = H;

F = X;
Fold = F;

t = 1;
told = 1;

if mu == 0
    mu_0 = 0.8*norm(G);
    mu_k = mu_0;
    mu_bar = 1e-5 * mu_k;
else
    mu_k = mu;
end

omega = 1/(nr_of_sparse_scatterers^2+1);

Res = G - P_OmegaC .* mysum(H,nr_of_sparse_scatterers^2+1);

Snorm = zeros(kmax,1);

nr_of_iterations = kmax;

for iterk = 1:kmax

    Fold = F;
    
    % Compute new iterate:

    [U,S,V] = svd(H{1} + omega*Res);
    F{1} = U*elementwise_soft_thresholding(S, omega*mu_k)*V';
    clear U V S

    for iteri = 1 : nr_of_sparse_scatterers^2
       incidentDir = indices(iteri,1);
       detectorPos = indices(iteri,2);
       Help1 = fft2(translOp(H{iteri+1} + omega*Res,sampling,z(:,incidentDir),z(:,detectorPos), kappa))/sqrt(nxhat*nd);
       Help2 = elementwise_soft_thresholding(Help1, omega*lambda*mu_k);
       F{iteri+1} = translOp(ifft2(Help2)*sqrt(nxhat*nd),sampling,-z(:,incidentDir),-z(:,detectorPos), kappa);
       clear Help1 Help2
    end

    % Stopping criterion:
    
    for iteri = 1 : (nr_of_sparse_scatterers^2+1)
        S = 1/omega * (H{iteri}-F{iteri}) + mysum(F,nr_of_sparse_scatterers^2+1) - mysum(H,nr_of_sparse_scatterers^2+1);
        Snorm(iterk) = Snorm(iterk) + norm(S,'fro')^2;
    end
    clear S
    Snorm(iterk) = sqrt(Snorm(iterk));
    if (Snorm(iterk) < tol && iterk > 1)
        Snorm = Snorm(1:iterk);
        nr_of_iterations = iterk;
        break
    end

    % weighting parameter for computing new auxiliary iterate:
    
    told = t;
    t = (1 + sqrt(4 * told^2 + 1))/2;
    
    Hold = H;

    % Compute new auxiliary iterate:
    
    H = myplus(myprod(1+(told-1)/t,F,nr_of_sparse_scatterers^2+1), myprod((1-told)/t,Fold,nr_of_sparse_scatterers^2+1), nr_of_sparse_scatterers^2+1);
    
    % Compute new residuum:

    Res = G - P_OmegaC .* mysum(H,nr_of_sparse_scatterers^2+1);
  
end

end

function Bsum = mysum(B, nr)

Bsum = B{1};
for iterk = 2:nr
    
    Bsum = Bsum + B{iterk};
    
end
end

function Z = myplus(X, Y, nr)

Z = cell(nr, 1);
for iterk = 1:nr
    
    Z{iterk} = X{iterk} + Y{iterk};
    
end
end

function Z = myprod(alpha, X, nr)

Z = cell(nr, 1);
for iterk = 1:nr
    
    Z{iterk} = alpha*X{iterk};
    
end
end
