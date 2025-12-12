%% Example 4.3:
%% Tests for varying distance between scatterer's components

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
znut = [26;-3]; % position of nut-shaped scatterer fixed
v = [-2;-1]/norm([-2;-1]);
dist = 30:5:80; % distances of kite-shaped scatterer to nut-shaped one

paramsnut = [1 .5 znut.' 45 4.5]; % parameters defining nut-shaped scatterer

nr_rep = length(dist);

%% Simulate far field data:

Fkite = cell(1,nr_rep);
Fall = cell(1,nr_rep);


for iteri = 1:nr_rep
    
    % Position of kite-shaped scatterer:
    zkite  = znut+dist(iteri)*v;
    
    % Parameters defining scatterers:
    paramskite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4];
    params = [paramskite; paramsnut];
    
    [Fkite{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, q(1), k, nobs, 0);
    [Fall{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite','nut'}, params, q, k, nobs, 0);

    clear params paramskite zkite z
    
    Fkite{iteri} = Fkite{iteri}*2*pi/nobs;
    Fall{iteri} = Fall{iteri}*2*pi/nobs;

end

%% Initializations for numerical reconstructions:

Akite = cell(nr_rep,1);
ALRC = cell(nr_rep,1);

relerrkite = zeros(nr_rep,1);
relerrLRC = zeros(nr_rep,1);

%% RPCP reconstruction:

kmax = 200;
lambda = 0.125;
tol = 1e-4;
X = {zeros(sampling.nxhat, sampling.nd), zeros(sampling.nxhat, sampling.nd)}; % initial guess

for iteri = 1:nr_rep

    mu = nobs*3e-4/lambda;

    zkite  = znut+dist(iteri)*v;

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda, mu, kmax, tol);
    
    ALRC{iteri} = A{1};
    Akite{iteri} = A{2};

    clear A

    relerrkite(iteri) = norm(Fkite{iteri}-Akite{iteri},'fro')/norm(Fkite{iteri},'fro');
    relerrLRC(iteri) = norm(Fall{iteri}-Fkite{iteri}-ALRC{iteri},'fro')/norm(Fall{iteri}-Fkite{iteri},'fro');

end

%% Plot relative reconstruction errors:

f = figure();
f.Position = [100 200 950 450];

semilogy(dist, relerrkite, '--*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(dist, relerrLRC, '--o','Color', 'blue','LineWidth', 1.5)
hold on

xlabel('$|\mathbf{c}_1-\mathbf{c}_2|$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

lgd = legend({'$\varepsilon_{\mathrm{rel}}^1$' , '$\varepsilon_{\mathrm{rel}}^{L}$'},'Interpreter','latex', 'NumColumns', 3);
lgd.Location = 'southeast';

xlim([dist(1) dist(end)])
ylim([10^(-3) 1])

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_N.eps -depsc

%% Plot geometry of scatterers:

figure()

d = (0:100)/100*2*pi;

% nut-shaped scatterer:
[x_nut,~,~,~] = kurve(100,'nut', paramsnut);
x_nut = [x_nut;x_nut(1,:)];
plot(x_nut(:,1),x_nut(:,2),'LineWidth', 2, 'Color', 'blue')
hold on

% kite-shaped scatterers:
for iteri = 1:nr_rep

    zkite  = znut+dist(iteri)*v;
    paramskite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4];
    
    [x_kite,~,~,~] = kurve(100, 'kite', paramskite);
    x_kite = [x_kite;x_kite(1,:)];
    plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', [1 0.7 0.7])
    hold on
    scatter(zkite(1),zkite(2),100,'Marker','+','LineWidth', 2,'MarkerEdgeColor', [.7 .7 .7])
    hold on

    clear zkite paramskite

end

zkite  = znut+dist(7)*v;
paramskite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4];

[x_kite,~,~,~] = kurve(100, 'kite', paramskite);
x_kite = [x_kite;x_kite(1,:)];
plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', 'red')
hold on
scatter(zkite(1),zkite(2),100,'Marker','+','LineWidth', 2,'MarkerEdgeColor', 'black')
hold on

grid on
axis equal

text(17,9,'$D_2$','Color','blue','FontSize',32, 'Interpreter','latex')
text(-37,-18,'$D_1$','Color','red','FontSize',32, 'Interpreter','latex')

xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')

set(gca,'XTick',[-60 -40 -20 0 20 40])
set(gca,'YTick',[-40 -20 0 20])

xlim([-70 50])
ylim([-50 30])

ax = gca;
ax.FontSize = 32;

print geometry_vary_N.eps -depsc