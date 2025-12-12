%% Example 4.2:
%% Tests for varying noise level

clear all
close all
clc

%% Initializations:

k = 1; % wave number

nobs = 170; % number of illumination and observation directions

sampling.nxhat = nobs;  % number of observation directions, i.e. number of rows 
sampling.nd = nobs; % number of illumination directions, i.e. number of columns
sampling.xhat = (0:sampling.nxhat-1)'/sampling.nxhat*2*pi; % observation directions
sampling.d = (0:sampling.nd-1)/sampling.nd*2*pi; % illumination directions

nr_of_scatterers = 2; 

q = [1, -.5]; % values of contrast

% Radii of scatterers:
Rnut = 5;
Rkite = 5;
R = [Rkite, Rnut]; 

% Positions of scatterers:
zkite  = [-28;-30];
znut = [26;-3];
z = [zkite, znut];

% Parameters defining scatterers:
paramskite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4];
paramsnut = [1 .5 znut.' 45 4.5];
params = [paramskite; paramsnut];

noiselevel = 0.01 * [1:10]; % relative noiselevel
nr_rep = length(noiselevel);

nr_tests = 15;

%% Simulate far field data:

for iteri = 1:nr_rep

    [Fkite, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, q(1), k, nobs, 0);
    [Fall, ~, ~] = evaluateFarfieldNystrom({'kite','nut'}, params, q, k, nobs, 0);
    
    Fkite = Fkite*2*pi/nobs;
    Fall = Fall*2*pi/nobs;

end

%% Initializations for numerical reconstructions:

relerrLRC_b = zeros(nr_rep,1);
relerrkite_b = zeros(nr_rep,1);

relerrLRC_w = zeros(nr_rep,1);
relerrkite_w = zeros(nr_rep,1);

delta_abs = zeros(1,nr_rep); % absolute noiselevel

%% RPCP Reconstructions

for iteri = 1:nr_rep
   
   rng("default")
    
   relerrLRC = zeros(nr_tests,1);
   relerrkite = zeros(nr_tests,1);

   delta_abs(iteri) = noiselevel(iteri)*norm(Fall,'fro');

   for iterk = 1:nr_tests

       % Generate noisy data:
       G = addNoise(Fall, noiselevel(iteri));

       % Reconstruction with RPCP formulation:
       lambda = 0.125;
       mu = nobs*(2e-4 + iteri*1e-4)/lambda;
       kmax = 200;
       tol = 1e-4;

       X = {zeros(sampling.nxhat, sampling.nd), zeros(sampling.nxhat, sampling.nd)}; % initial guess

       [A,~,~] = RPCP(G, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda, mu, kmax, tol);
       
       ALRC = A{1};
       Akite = A{2};

       relerrkite(iterk) = norm(Fkite-Akite,'fro')/norm(Fkite,'fro');
       relerrLRC(iterk) = norm(Fall-Fkite-ALRC,'fro')/norm(Fall-Fkite,'fro');

       clear ALRC Akite
   end

   relerrLRC_b(iteri) = min(relerrLRC);
   relerrkite_b(iteri) = min(relerrkite);
    
   relerrLRC_w(iteri) = max(relerrLRC);
   relerrkite_w(iteri) = max(relerrkite);
end


%% Plot worst relative reconstruction errors:

f = figure();
f.Position = [100 200 950 450];

semilogy(noiselevel, relerrkite_w, '--*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(noiselevel, relerrLRC_w,'--o','Color', 'blue','LineWidth', 1.5)
hold on
semilogy(noiselevel, 2*noiselevel,  '-','Color', [0.6, 0.6, 0.6],'LineWidth', 1.5)
hold on 

xlabel('$\delta_{\mathrm{rel}}$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

legend({'$\varepsilon_{\mathrm{rel}}^{1}$', '$\varepsilon_{\mathrm{rel}}^{L}$','$\mathcal{O}(\delta_{\mathrm{rel}})$'}, ...
    'Interpreter','latex', 'Location','southeast', 'NumColumns', 3)

xlim([noiselevel(1) noiselevel(end)])
ylim([10^(-3) 1])

set(gca,'XTick',[0.02 0.04 0.06 0.08 0.1]);
set(gca,'YTick',[1e-2 1e-1 1e0]);

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_noiselevel_worst.eps -depsc