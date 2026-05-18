%% Example 4.4:
%% Tests for varying intensity of the contrast

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

q = [1, .5]; % values of contrast
sigma = [.2:.2:3]; % intensities of contrast

% Radii of scatterers:
Rnut = 5;
Rkite = 5;
R = [Rkite, Rnut]; 

% Positions of scatterers:
znut = [26;-3]; % position of nut-shaped scatterer fixed
v = [-2;-1]/norm([-2;-1]);
dist = 60; % distances of kite-shaped scatterer to nut-shaped one
zkite = znut + dist * v;

paramskite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4]; % parameters defining kite-shaped scatterer
paramsnut = [1 .5 znut.' 45 4.5]; % parameters defining nut-shaped scatterer
params = [paramskite; paramsnut];

nr_rep = length(sigma);

%% Simulate far field data:

Gexact = cell(1,nr_rep);
Fnut = cell(1,nr_rep);
Fkite = cell(1,nr_rep);
Fall = cell(1,nr_rep);
model_err = zeros(nr_rep,1);

for iteri = 1:nr_rep

    [Fkite{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, sigma(iteri)*q(1), k, nobs, 0);
    [Fnut{iteri}, ~, ~] = evaluateFarfieldNystrom({'nut'}, paramsnut, sigma(iteri)*q(2), k, nobs, 0);
    [Fall{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite','nut'}, params, sigma(iteri)*q, k, nobs, 0);
    
    Fkite{iteri} = Fkite{iteri}*2*pi/nobs;
    Fnut{iteri} = Fnut{iteri}*2*pi/nobs;
    Gexact{iteri} = Fall{iteri}*2*pi/nobs;

    model_err(iteri) = norm(Fkite{iteri} + Fnut{iteri} - Gexact{iteri},'fro') / norm(Gexact{iteri},'fro');

end

%% Initializations for numerical reconstructions:

Akite = cell(nr_rep,1);
ALRC = cell(nr_rep,1);

relerrkite = zeros(nr_rep,1);
relerrLRC = zeros(nr_rep,1);

%% RPCP reconstruction:

kmax = 400;

lambda = 0.125;

mu = [1e-5 2e-5 2.5e-5 3e-5 3.5e-5 4.5e-5 4e-5 4e-5 4e-5 2e-5 2e-5 2e-5 2e-5 2e-5 2e-5] * nobs/lambda;

tol = 1e-3;

X = {zeros(sampling.nd, sampling.nd), zeros(sampling.nd, sampling.nd)};

for iteri = 1:nr_rep

   [A,~,Snorm{iteri}] = RPCP(Gexact{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda, mu(iteri), kmax, tol);
   
   ALRC{iteri} = A{1};
   Akite{iteri} = A{2};

   clear A

   relerrkite(iteri) = norm(Fkite{iteri}-Akite{iteri},'fro')/norm(Fkite{iteri},'fro');
   relerrLRC(iteri) = norm(Gexact{iteri}-Fkite{iteri}-ALRC{iteri},'fro')/norm(Gexact{iteri}-Fkite{iteri},'fro');

end

%% Plot relative reconstruction errors:

f = figure();
f.Position = [100 200 950 450];

semilogy(sigma, relerrkite, '--*','Color', 'red','LineWidth', 1.5)
hold on 
semilogy(sigma, relerrLRC, '--o','Color', 'blue','LineWidth', 1.5)
hold on
semilogy(sigma, model_err, '-','Color', [0.6, 0.6, 0.6],'LineWidth', 1.5)

xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

legend({'$\varepsilon_{\mathrm{rel}}^{1}$', '$\varepsilon_{\mathrm{rel}}^{L}$'}, ...
    'Interpreter','latex', 'Location','southeast', 'NumColumns', 3)

xlim([sigma(1) sigma(end)])
ylim([10^(-3) 1])

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_q.eps -depsc
