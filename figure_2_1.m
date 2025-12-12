clear all
close all 
clc

%% Initializations:

k = 1;  % wave number

M = 170; % number of illumination and observation directions

sampling.nxhat = M;  % number of observation directions, i.e. number of rows 
sampling.nd = M; % number of illumination directions, i.e. number of columns
sampling.xhat = (0:sampling.nxhat-1)'/sampling.nxhat*2*pi; % observation directions
sampling.d = (0:sampling.nd-1)/sampling.nd*2*pi; % illumination directions

% Values of piecewise constant constrast inside scatterer's components:
q = [1,-0.5];

% Positions of scatterer's components:
zkite  = [-28;-30];
znut = [26;-3];

% Parameters defining scatterer's components
parkite = [1 .5 [zkite.'+[-0.6,0.7]] 135 3.4];
parnut = [1 .5 znut.' 45 4.5];

% Radii of scatterer's components:
Rkite = 5;
Rnut = 5;

%% Plot geometry of scatterers:

s = [50 50 600 600];
figure('Renderer', 'painters', 'Position', s)

d = (0:100)/100*2*pi;

[x_kite,~,~,~] = kurve(100, 'kite', parkite);
x_kite = [x_kite;x_kite(1,:)];
plot(x_kite(:,1),x_kite(:,2),'LineWidth', 2, 'Color', 'blue')
hold on

[x_nut,~,~,~] = kurve(100, 'nut', parnut);
x_nut = [x_nut;x_nut(1,:)];
plot(x_nut(:,1),x_nut(:,2),'LineWidth', 2, 'Color', 'blue')
hold on

plot(Rkite*cos(d)+zkite(1) ,Rkite*sin(d)+zkite(2),'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'black')
hold on
scatter(zkite(1),zkite(2),100,'Marker','+','LineWidth', 1.5,'MarkerEdgeColor', 'black')
hold on
plot(Rnut*cos(d)+znut(1) ,Rnut*sin(d)+znut(2),'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'black')
hold on
scatter(znut(1),znut(2),100,'Marker','+','LineWidth', 1.5,'MarkerEdgeColor', 'black')

text(-20,-30,'$D_1$','Color','blue','FontSize',32, 'Interpreter','latex')
text(10,3,'$D_2$','Color','blue','FontSize',32, 'Interpreter','latex')

set(gca,'XTick',[-30 -15 0 15 30])
set(gca,'YTick',[-45 -30 -15 0 15])

xlim([-40 40])
ylim([-55 25])

grid on
ax = gca;
ax.FontSize = 32;

xlabel('$x_1$', 'Interpreter', 'LaTex', 'Fontsize', 32)
ylabel('$x_2$', 'Interpreter', 'LaTex', 'Fontsize', 32)

print figure_geometry.eps -depsc

%% Generate far field data:

[F] = evaluateFarfieldNystrom({'kite', 'nut'}, [parkite; parnut], q, k, sampling.nxhat, 0);
[F1, ~, ~] = evaluateFarfieldNystrom({'kite'}, parkite, q(1), k, sampling.nxhat, 0);
F21 = F - F1;

%% Plot modulated Fourier coefficients of F1:

s = [50 50 600 600];
figure('Renderer', 'painters', 'Position', s)

N = ceil(exp(1)/2*k*Rkite); % truncation index

Ft1 = translOp(F1, sampling, zkite, zkite, k);

FF1 = 1/M^2*fftshift(fft2(Ft1));
FF1 = [FF1, FF1(:,1)];
FF1 = [FF1; FF1(1,:)];
FF1 = flipud(FF1);

L = 80;

imagesc((-L:L),(-L:L),abs(FF1(M/2+1-L:M/2+1+L,M/2+1-L:M/2+1+L)))
hold on

plot([-N,N],[N,N],'--k', 'LineWidth', 1.5)
hold on
plot([-N,N],[-N,-N],'--k', 'LineWidth', 1.5)
hold on
plot([N,N],[-N,N],'--k', 'LineWidth', 1.5)
hold on
plot([-N,-N],[-N,N],'--k', 'LineWidth', 1.5)

load('MyCMap.mat')
cmp = CMap;
colormap(cmp);

set(gca,'ColorScale','log')
cx = clim;
clim([10^(-4) cx(2)])

set(gca,'XTick',[-60 -30 0 30 60]);
set(gca,'YTick',[-60 -30 0 30 60]);

grid on

xlabel('$m$', 'Interpreter', 'LaTex', 'Fontsize', 27)
ylabel('$n$', 'Interpreter', 'LaTex', 'Fontsize', 27)

axis square

colorbar

set(gca,'Fontsize',27)
set(gca,'YDir','normal')

print figure_exp_coeff_F1.eps -depsc

%% Plot modulated Fourier coefficients of F21:

s = [50 50 600 600];
figure('Renderer', 'painters', 'Position', s)

N = ceil(exp(1)/2*k*Rnut); % truncation index

Ft21 = translOp(F21, sampling, znut, znut, k);

FF21 = 1/M^2*fftshift(fft2(Ft21));
FF21 = [FF21, FF21(:,1)];
FF21 = [FF21; FF21(1,:)];
FF21 = flipud(FF21);

L = 80;

imagesc((-L:L),(-L:L),abs(FF21(M/2+1-L:M/2+1+L,M/2+1-L:M/2+1+L)))
hold on

plot([-M/2,-N],[N,N],'--k', 'LineWidth', 1.5)
hold on
plot([N,M/2],[N,N],'--k', 'LineWidth', 1.5)
hold on
plot([-M/2,-N],[-N,-N],'--k', 'LineWidth', 1.5)
hold on
plot([N,M/2],[-N,-N],'--k', 'LineWidth', 1.5)
hold on
plot([N,N],[-M/2,-N],'--k', 'LineWidth', 1.5)
hold on
plot([N,N],[N,M/2],'--k', 'LineWidth', 1.5)
hold on
plot([-N,-N],[-M/2,-N],'--k', 'LineWidth', 1.5)
hold on
plot([-N,-N],[N,M/2],'--k', 'LineWidth', 1.5)

set(gca,'ColorScale','log')
cx = clim;
clim([10^(-4) cx(2)])

load('MyCMap.mat')
cmp = CMap;
colormap(cmp);

set(gca,'XTick',[-60 -30 0 30 60]);
set(gca,'YTick',[-60 -30 0 30 60]);

grid on

xlabel('$m$', 'Interpreter', 'LaTex', 'Fontsize', 27)
ylabel('$n$', 'Interpreter', 'LaTex', 'Fontsize', 27)

axis square

colorbar

set(gca,'Fontsize',27)
set(gca,'YDir','normal')

print figure_exp_coeff_F21.eps -depsc