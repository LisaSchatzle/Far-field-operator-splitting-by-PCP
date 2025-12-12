%% Example 4.3:
%% Tests for varying sizes of scatterer's components

clear all
close all
clc

%% Tests for varying size of nut-shaped scatterer:
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
Rnut = 1:13;
Rkite = 5;

% Positions of scatterers:
zkite  = [-28;-30];
znut = [26;-3];
z = [zkite, znut];

% s = [-0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7; -0.6,0.7];

paramskite = [1 .5 [zkite.'+[-1.1,1.1]*.7] 135 3.65]; % parameters defining kite-shaped scatterer

nr_rep = length(Rnut);

%% Simulate far field data for varying size of nut-shaped scatterer:

[Fkite, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, q(1), k, nobs, 0);
Fkite = Fkite*2*pi/nobs;

Fall = cell(1,nr_rep);

for iteri = 1:nr_rep
    
    K = 4.7/5*Rnut(iteri); % scaling factors for size of nut-shaped scatterer
    paramsnut = [1 .5 znut.' 45 K];
    params = [paramskite; paramsnut];
    
    [Fall{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite','nut'}, params, q, k, nobs, 0);

    clear params K paramsnut
    
    Fall{iteri} = Fall{iteri}*2*pi/nobs;

end

%% Initializations for numerical reconstructions for varying size of nut-shaped scatterer:

Akite = cell(nr_rep,1);
ALRC = cell(nr_rep,1);

relerrkite = zeros(nr_rep,1);
relerrLRC = zeros(nr_rep,1);

%% RPCP reconstruction for varying size of nut-shaped scatterer:

kmax = 200;
lambda = nobs * [7.5 7.4 7.3 7.2 7.4 7.7 7.7 7.5 7.8 8.2 8.5 8.5 7.9]*1e-4;
tol = 1e-4;
X = {zeros(sampling.nxhat, sampling.nd), zeros(sampling.nxhat, sampling.nd)}; % initial guess

for iteri = 1:nr_rep

    mu = nobs * 3e-4/lambda(iteri);

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda(iteri), mu, kmax, tol);
    
    ALRC{iteri} = A{1};
    Akite{iteri} = A{2};

    clear A

    relerrkite(iteri) = norm(Fkite-Akite{iteri},'fro')/norm(Fkite,'fro');
    relerrLRC(iteri) = norm(Fall{iteri}-Fkite-ALRC{iteri},'fro')/norm(Fall{iteri}-Fkite,'fro');

end

%% Plot relative reconstruction errors for varying size of nut-shaped scatterer:

f = figure();
f.Position = [100 200 950 450];

semilogy(Rnut, relerrkite, '--*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(Rnut, relerrLRC, '--o','Color', 'blue','LineWidth', 1.5)
hold on

xlabel('$R_2$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

lgd = legend({'$\varepsilon_{\mathrm{rel}}^1$' , '$\varepsilon_{\mathrm{rel}}^{L}$'},'Interpreter','latex', 'NumColumns', 3);
lgd.Location = 'southeast';

xlim([Rnut(1) Rnut(end)])
ylim([10^(-3) 1])

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_N.eps -depsc

%% Geometry of scatterers for varying size of nut-shaped scatterer:

figure()

% Nut-shaped scatterers:

Rnut = [1 3 5 7 9 11 13];

for iteri = 1:7

    K = 4.7/5*Rnut(iteri); % scaling factors for size of nut-shaped scatterer
    paramsnut = [1 .5 znut.' 45 K];
    
    [x_nut,~,~,~] = kurve(100,'nut', paramsnut);
    x_nut = [x_nut;x_nut(1,:)];
    plot(x_nut(:,1),x_nut(:,2),'LineWidth', 2, 'Color', [.7 .7 1])
    hold on

end

[x_nut,~,~,~] = kurve(100,'nut', [1 .5 znut.' 45 4.7]);
x_nut = [x_nut;x_nut(1,:)];
plot(x_nut(:,1),x_nut(:,2),'LineWidth', 2, 'Color', [0 0 1])
hold on

% Kite-shaped scatterer:

[x_kite,~,~,~] = kurve(100, 'kite', paramskite);
x_kite = [x_kite;x_kite(1,:)];
plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', [1 0 0])
hold on

scatter(zkite(1),zkite(2),100,'Marker','+','LineWidth', 2,'MarkerEdgeColor', 'black')
hold on

grid on
axis equal

text(17,15,'$D_2$','Color','blue','FontSize',32, 'Interpreter','latex')
text(-37,-19,'$D_1$','Color','red','FontSize',32, 'Interpreter','latex')

xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')

set(gca,'XTick',[-60 -40 -20 0 20 40])
set(gca,'YTick',[-40 -20 0 20])

xlim([-70 50])
ylim([-50 30])

ax = gca;
ax.FontSize = 32;

print geometry_vary_size_nut.eps

%% Tests for varying size of nut-shaped scatterer:

clear all
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
Rkite = 1:13;

% Positions of scatterers:
zkite  = [-28;-30];
znut = [26;-3];
z = [zkite, znut];

paramsnut = [1 .5 znut.' 45 4.7]; % parameters defining nut-shaped scatterer

nr_rep = length(Rkite);

%% Simulate far field data for varying size of kite-shaped scatterer:

Fkite = cell(1,nr_rep);
Fall = cell(1,nr_rep);

for iteri = 1:nr_rep
    
    K = 3.65/5*Rkite(iteri); % scaling factors for size of kite-shaped scatterer
    paramskite = [1 .5 [zkite.'+[-1.1,1.1]*.7/5*Rkite(iteri)] 135 K];
    params = [paramskite; paramsnut];
    
    [Fkite{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite'}, paramskite, q(1), k, nobs, 0);
    [Fall{iteri}, ~, ~] = evaluateFarfieldNystrom({'kite','nut'}, params, q, k, nobs, 0);

    clear params paramskite K
    
    Fkite{iteri} = Fkite{iteri}*2*pi/nobs;
    Fall{iteri} = Fall{iteri}*2*pi/nobs;

end

%% Initializations for numerical reconstructions for varying size of kite-shaped scatterer:

Akite = cell(nr_rep,1);
ALRC = cell(nr_rep,1);

relerrkite = zeros(nr_rep,1);
relerrLRC = zeros(nr_rep,1);

%% RPCP reconstruction for varying size of kite-shaped scatterer:

kmax = 200;
lambda = nobs * [12.9 9.8 8.9 8.1 7.5 6.9 6.8 6.8 6.3 6.2 5.9 6 5.6] * 1e-4;
tol = 1e-4;
X = {zeros(sampling.nxhat, sampling.nd), zeros(sampling.nxhat, sampling.nd)}; % initial guess

for iteri = 1:nr_rep

    mu = nobs * 3e-4/lambda(iteri);

    [A,~,~] = RPCP(Fall{iteri}, X, sampling, k, ones(sampling.nxhat, sampling.nd), zkite, lambda(iteri), mu, kmax, tol);
    
    ALRC{iteri} = A{1};
    Akite{iteri} = A{2};

    clear A

    relerrkite(iteri) = norm(Fkite{iteri}-Akite{iteri},'fro')/norm(Fkite{iteri},'fro');
    relerrLRC(iteri) = norm(Fall{iteri}-Fkite{iteri}-ALRC{iteri},'fro')/norm(Fall{iteri}-Fkite{iteri},'fro');
end

%% Plot relative reconstruction errors for varying size of kite-shaped scatterer:

f = figure();
f.Position = [100 200 950 450];

semilogy(Rkite, relerrkite, '--*','Color', 'red','LineWidth', 1.5)
hold on
semilogy(Rkite, relerrLRC, '--o','Color', 'blue','LineWidth', 1.5)
hold on

xlabel('$R_1$', 'Interpreter', 'latex')
ylabel('$\varepsilon_{\mathrm{rel}}$', 'Interpreter', 'latex')

lgd = legend({'$\varepsilon_{\mathrm{rel}}^1$' , '$\varepsilon_{\mathrm{rel}}^{L}$'},'Interpreter','latex', 'NumColumns', 3);
lgd.Location = 'southeast';

xlim([Rkite(1) Rkite(end)])
ylim([10^(-3) 1])

grid on

ax = gca;
ax.FontSize = 30;

print errors_vary_size_kite.eps -depsc

%% Geometry for varying size of kite-shaped scatterer:

figure()

% Nut-shaped scatterer:

[x_nut,~,~,~] = kurve(100,'nut', paramsnut);
x_nut = [x_nut;x_nut(1,:)];
plot(x_nut(:,1),x_nut(:,2),'LineWidth', 2, 'Color', [0 0 1])
hold on

% Kite-shaped scatterer:

Rkite = [1 3 5 7 9 11 13];

for iteri = 1:7

    K = 3.65/5*Rkite(iteri);
    paramskite = [1 .5 [zkite.'+[-1.1,1.1]*.7/5*Rkite(iteri)] 135 K];
    [x_kite,~,~,~] = kurve(100, 'kite', paramskite);
    x_kite = [x_kite;x_kite(1,:)];
    plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', [1 .7 .7])
    hold on

end

K = 3.65/5*Rkite(3);
paramskite = [1 .5 [zkite.'+[-1.1,1.1]*.7/5*Rkite(3)] 135 K];
[x_kite,~,~,~] = kurve(100, 'kite', paramskite);
x_kite = [x_kite;x_kite(1,:)];
plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', [1 0 0])
hold on

scatter(zkite(1),zkite(2),100,'Marker','+','LineWidth', 2,'MarkerEdgeColor', 'black')
hold on

grid on
axis equal

text(17,9,'$D_2$','Color','blue','FontSize',32, 'Interpreter','latex')
text(-37,-14,'$D_1$','Color','red','FontSize',32, 'Interpreter','latex')

xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')

set(gca,'XTick',[-60 -40 -20 0 20 40])
set(gca,'YTick',[-40 -20 0 20])

xlim([-70 50])
ylim([-50 30])

ax = gca;
ax.FontSize = 32;

print geometry_vary_size_kite.eps -depsc