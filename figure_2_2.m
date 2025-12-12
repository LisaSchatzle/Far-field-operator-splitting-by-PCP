clear all
close all
clc

%% Plot of beta and upper bound for N=5, 5<=t<=1500:

N = 5;
t = 5:5:1500;
b = 0.7595; % constant in upper bound of Bessel function

beta = zeros(1,length(t));
bound = zeros(1,length(t));

A = zeros(1,length(t));
B = zeros(1,length(t));

for itert = 1:length(t)
    for itern0 = -10*t(itert):10*t(itert)
        term = besselj(itern0,t(itert))^2;
        for itern = 1:N
            term = term + besselj(itern0+itern,t(itert))^2 + besselj(itern0-itern,t(itert))^2;
        end
        if term >= beta(itert)
            beta(itert) = max(term, beta(itert));
            n0 = itern0;
        end
    end
    bound(itert) = b*sqrt(2*N+1)/((t(itert))^(1/2));
end

f = figure();
f.Position = [100 200 950 450];

loglog([5, 1500],[0.1,0.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
loglog([192, 192],[0.03,1.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
loglog([636, 636],[0.03,1.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
plot(t,beta,'b:','linewidth',2)
hold on
plot(t,bound,'r--','linewidth',2)
hold on

title('$N_2=5$', 'Interpreter', 'LaTex')
xlabel('$t$', 'Interpreter', 'LaTex')
legend({'','','','$\beta_{N_1}^{\mathbf{c}_1-\mathbf{c}_2}$','$b\sqrt{\frac{2N_2+1}t}$'},'Interpreter', 'LaTex','Location','southwest')

set(gcf, 'Color', 'w')
set(gca, 'Color', 'w') 

grid on

set(gca,'Fontsize',26)
set(gca, 'Linewidth', 1.5);

set(gca,'YTick',[0.05, 0.1, 0.2, 0.5, 1])
set(gca,'XTick',[10, 100, 1000])
ax = gca;
ax.XLim = [5, 1500];
ax.YLim = [0.03,1.1];

print beta_upperbound_N5.eps -depsc

%% Plot of beta and upper bound for N=10, 5<=t<=1500:

N = 10;
t = 5:5:1500;

b = 0.7595; % constant in upper bound of Bessel function

beta = zeros(1,length(t));
bound = zeros(1,length(t));

for itert = 1:length(t)
    for itern0 = -10*t(itert):10*t(itert)
        term = besselj(itern0,t(itert))^2;
        for itern = 1:N
            term = term + besselj(itern0+itern,t(itert))^2 + besselj(itern0-itern,t(itert))^2;
        end
        if term >= beta(itert)
            beta(itert) = max(term, beta(itert));
            n0 = itern0;
        end
    end
    bound(itert) = b*sqrt(2*N+1)/((t(itert))^(1/2));
end

f = figure();
f.Position = [100 200 950 450];


loglog([5, 1500],[0.1,0.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
loglog([373, 373],[0.03,1.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
loglog([1213, 1213],[0.03,1.1],'-','color',[.3,.3,.3],'linewidth',2)
hold on
plot(t,beta,'b:','linewidth',2)
hold on
plot(t,bound,'r--','linewidth',2)

title('$N_2=10$', 'Interpreter', 'LaTex')
xlabel('$t$', 'Interpreter', 'LaTex')
legend({'','','','$\beta_{N_1}^{\mathbf{c}_1-\mathbf{c}_2}$','$b\sqrt{\frac{2N_2+1}{t}}$'},'Interpreter', 'LaTex','Location','southwest')

set(gcf, 'Color', 'w')
set(gca, 'Color', 'w') 

grid on

set(gca,'Fontsize',26)
set(gca, 'Linewidth', 1.5);

set(gca,'XTick',[10, 100, 1000])
set(gca,'YTick',[0.05, 0.1, 0.2, 0.5, 1])
ax = gca;
ax.XLim = [5, 1500];
ax.YLim = [0.03,1.1];
ax.Toolbar.Visible = 'off';

print beta_upperbound_N10.eps -depsc
