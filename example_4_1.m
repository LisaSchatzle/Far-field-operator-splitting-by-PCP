%% Example 4.1:
%% Tests for varying number of illumination and oservation directions
%% Tests on choice of coupling parameter lambda

clear all
close all
clc

%% Initializations:

k = 1;  % wave number

nobs = 130:20:350; % number of illumination and observation directions

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
    
nr_rep = length(nobs); % number of test runs

%% Simulate far field data:

Fkite = cell(1,nr_rep);
Fall = cell(1,nr_rep);

for iteri = 1:nr_rep

    [Fkite{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, q(1), k, nobs(iteri), 0);
    [Fall{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite', 'nut'}, params, q, k, nobs(iteri), 0);
    
    Fkite{iteri} = Fkite{iteri}*2*pi/nobs(iteri);
    Fall{iteri} = Fall{iteri}*2*pi/nobs(iteri);

end

%% Initializations for numerical reconstructions:

Akite1 = cell(nr_rep,1);
ALRC1 = cell(nr_rep,1);

relerrkite1 = zeros(nr_rep,1);
relerrLRC1 = zeros(nr_rep,1);

Akite2 = cell(nr_rep,1);
ALRC2 = cell(nr_rep,1);

relerrkite2 = zeros(nr_rep,1);
relerrLRC2 = zeros(nr_rep,1);

%% RPCP for standard choice of lambda:

kmax = 200;
lambda = 1./(nobs.^(1/2));

for iteri = 1:nr_rep

    sampling.nxhat = nobs(iteri);  % number of observation directions, i.e. number of rows 
    sampling.nd = nobs(iteri); % number of illumination directions, i.e. number of columns
    sampling.xhat = (0:sampling.nxhat-1)'/sampling.nxhat*2*pi; % detector positions
    sampling.d = (0:sampling.nd-1)/sampling.nd*2*pi; % illumination directions

    mu = nobs(iteri)*3e-4/lambda(iteri);

    tol = 1e-4;

    X = {zeros(sampling.nxhat, sampling.nd), zeros(sampling.nxhat, sampling.nd)}; % initial guess

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda(iteri), mu, kmax, tol);
    
    ALRC1{iteri} = A{1};
    Akite1{iteri} = A{2};

    clear A sampling

    relerrkite1(iteri) = norm(Akite1{iteri}-Fkite{iteri},'fro')/norm(Fkite{iteri},'fro');
    relerrLRC1(iteri) = norm(Fall{iteri}-Fkite{iteri}-ALRC1{iteri},'fro')/norm(Fall{iteri}-Fkite{iteri},'fro');

end

%% RPCP with optimally chosen lambda:

kmax = 200;
lambda = nobs.*[10 8.8 7.4 6.7 6 5.6 5.1 4.8 4.4 4.2 4 3.8]*1e-4;

for iteri = 1:nr_rep

    sampling.nxhat = nobs(iteri);  % number of observation directions, i.e. number of rows 
    sampling.nd = nobs(iteri); % number of illumination directions, i.e. number of columns
    sampling.xhat = (0:sampling.nxhat-1)'/sampling.nxhat*2*pi; % detector positions
    sampling.d = (0:sampling.nd-1)/sampling.nd*2*pi; % illumination directions

    mu = nobs(iteri)*3e-4/lambda(iteri);

    tol = 1e-4;

    X = {zeros(sampling.nd, sampling.nd), zeros(sampling.nd, sampling.nd)}; % initial guess

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda(iteri), mu, kmax, tol);
    
    ALRC2{iteri} = A{1};
    Akite2{iteri} = A{2};

    clear A sampling

    relerrkite2(iteri) = norm(Akite2{iteri}-Fkite{iteri},'fro')/norm(Fkite{iteri},'fro');
    relerrLRC2(iteri) = norm(Fall{iteri}-Fkite{iteri}-ALRC2{iteri},'fro')/norm(Fall{iteri}-Fkite{iteri},'fro');
end

%% Plot relative reconstruction errors for both choices of lambda:

f = figure();
f.Position = [100 200 950 450];

semilogy(nobs, relerrkite1, ':*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(nobs, relerrLRC1, ':o','Color', 'blue','LineWidth', 1.5)
hold on 
semilogy(nobs, relerrkite2, '--*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(nobs, relerrLRC2, '--o','Color', 'blue','LineWidth', 1.5)
hold on 
semilogy(nobs, nobs/nobs(1)*3e-2, '-','Color', [0.6, 0.6, 0.6],'LineWidth', 1.5)

xlabel('$M$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

lgd = legend({'$\varepsilon_{\mathrm{rel}}^1,\lambda =M^{-1/2}$' , '$\varepsilon_{\mathrm{rel}}^{L},\lambda =M^{-1/2}$','$\varepsilon_{\mathrm{rel}}^1$, optimal $\lambda$', '$\varepsilon_{\mathrm{rel}}^{L}$, optimal $\lambda$','$\mathcal{O}(M)$'},'Interpreter','latex', 'NumColumns', 3);
lgd.Location = 'southeast';

xlim([nobs(1) nobs(end)])
ylim([10^(-4) 1])

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_N.eps -depsc

%% lambda test for nobs=170:

lambda = (0.07:0.005:0.15);
nr_rep = length(lambda);

Akite3 = cell(nr_rep,1);
ALRC3 = cell(nr_rep,1);

relerrkite3 = zeros(nr_rep,1);
relerrLRC3 = zeros(nr_rep,1);

iteri = 3;

kmax = 200;

sampling.nxhat = nobs(iteri);  % number of observation directions, i.e. number of rows 
sampling.nd = nobs(iteri); % number of illumination directions, i.e. number of columns
sampling.xhat = (0:sampling.nxhat-1)'/sampling.nxhat*2*pi; % detector positions
sampling.d = (0:sampling.nd-1)/sampling.nd*2*pi; % illumination directions

for iterl = 1:nr_rep

    mu = nobs(iteri)*3e-4/lambda(iteri);

    tol = 1e-4;

    X = {zeros(sampling.nd, sampling.nd), zeros(sampling.nd, sampling.nd)};

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda(iterl), mu, kmax, tol);
    
    ALRC3{iterl} = A{1};
    Akite3{iterl} = A{2};

    clear A 

    relerrkite3(iterl) = norm(Akite3{iterl}-Fkite{iteri},'fro')/norm(Fkite{iteri},'fro');
    relerrLRC3(iterl) = norm(Fall{iteri}-Fkite{iteri}-ALRC3{iterl},'fro')/norm(Fall{iteri}-Fkite{iteri},'fro');

end

%% Plot of relative reconstruction errors for nobs=170 depending on lambda:

f = figure();
f.Position = [100 200 950 450];

plot([0.0767 0.0767], [0 0.25], ':','Color', [0.6, 0.6, 0.6],'LineWidth', 1.5)
hold on
plot([0.125 0.125], [0 0.25], '--','Color', [0.6, 0.6, 0.6],'LineWidth', 1.5)
hold on
plot(lambda, relerrkite3, '-*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(lambda, relerrLRC3, '-o','Color', 'blue','LineWidth', 1.5)
hold on 

xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

lgd = legend({'','','$\varepsilon_{\mathrm{rel}}^1$' , '$\varepsilon_{\mathrm{rel}}^{L}$'},'Interpreter','latex', 'NumColumns', 1);
lgd.Location = 'northeast';

xlim([lambda(1) lambda(end)])
ylim([0 0.25])

set(gca,'XTick',[4.5 7.4]*1e-4);
set(gca,'XTickLabel',{'0.0767','0.1258'})

grid on

ax = gca;
ax.FontSize = 30;

print lambda_scan_M170.eps -depsc