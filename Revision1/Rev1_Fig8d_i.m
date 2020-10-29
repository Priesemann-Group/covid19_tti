%%% Rtcrit(eta) %%%

clear all
close all
clc

%% Definici?n de parametros "Default"

xi = 0.15; xi0 = xi;
varphi = 0.2;
xi = xi + (1-xi)*varphi; xim = 1-xi;
lambda_s = 0.1;
lambda_r = 0;
etaesc = 0.66;
Phi = 15;
lmax = 1;
etamax = 1;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
eta = 0:etamax/100:etamax;

%% Definici?n par?metros barrido
kmax = length(eta);
Rt = zeros(size(eta));

Rtcmax = fzero(@(Rtk) maxVp(Pools(xi,nu,Rtk,Gamma,lambda_s,lambda_r,etaesc,epsilon)),2); %eta
eta = etaesc;

%% calculo 10%

Rcrit10 = 1.1*Rtcmax;
Rcrit20 = 1.2*Rtcmax;
%eta
eta10 = fzero(@(eta) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
eta20 = fzero(@(eta) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
%lambda_s
lambda_s10 = fzero(@(lambda_s) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
lambda_s20 = fzero(@(lambda_s) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
%eta
lambda_r10 = fzero(@(lambda_r) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
lambda_r20 = fzero(@(lambda_r) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
%eta
xi10 = fzero(@(xi) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
xi20 = fzero(@(xi) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
varphi10 = (xi10-xi0)/(1-xi0);
varphi20 = (xi20-xi0)/(1-xi0);
%eta
ep10 = fzero(@(epsilon) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
ep20 = fzero(@(epsilon) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
%eta
nu10 = fzero(@(nu) maxVp(Pools(xi,nu,Rcrit10,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
nu20 = fzero(@(nu) maxVp(Pools(xi,nu,Rcrit20,Gamma,lambda_s,lambda_r,eta,epsilon)),2);

ETA = [eta eta10 eta20];
Lambda_s = [lambda_s lambda_s10 lambda_s20];
Lambda_r = [lambda_r lambda_r10 lambda_r20];
Varphi = [varphi varphi10 varphi20];
Epsilon = [epsilon ep10 ep20];
Nu = [nu nu10 nu20];

%% Ploteo

fact_axis = 2;
fact_label = 2.5;
fact_curva = 2;
siz = 15;

load('DefColors.mat')
%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = Lambda_s;
C = Default(1,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
xlabel('Increase in $R_{\mbox{crit}}^{H}$','interpreter','latex','FontSize',15*fact_label)
ylabel('Sympt. test. $\lambda_s$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';

%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = Lambda_r;
C = Default(2,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
hold on
plot([0 8],[0.002 0.002],'r--','HandleVisibility','off','LineWidth',2*fact_curva)
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
xlabel('Increase in $R_{\mbox{crit}}^{H}$','interpreter','latex','FontSize',15*fact_label)
ylabel('Rand. test. $\lambda_r$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';

%%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = ETA;
C = Default(3,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
hold on
plot([0 8],[1 1],'r--','HandleVisibility','off','LineWidth',2*fact_curva)
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
ylim([0 1.2])
xlabel('Increase in $R_{\mbox{crit}}^{H}$','interpreter','latex','FontSize',15*fact_label)
ylabel('Tracing eff. $\eta$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';

%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = Epsilon;
C = Default(4,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
ylim([0 0.15])
xlabel('Increase in $R_{\mbox{crit}}^{H}$','interpreter','latex','FontSize',15*fact_label)

ylabel('Miss. cont. (tr.) $\epsilon$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';

%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = Nu;
C = Default(5,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
ylim([0 0.15])
xlabel('Increase in $R_{\mbox{crit}}^{H}$','interpreter','latex','FontSize',15*fact_label)
ylabel('Isol. factor (tr.) $\nu$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';


%%%%%%%%
figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

v = Varphi;
C = Default(6,:);
for i = 1:length(v)
    hold on
    bar(1+2*(i-1),v(i),1.5,'FaceColor',C,'FaceAlpha',i/length(v))
end
set(gca, 'XTick', [1 3 5])
set(gca, 'XTickLabel', {'$+0\%$' '$+10\%$' '$+20\%$'})
set(gca,'FontSize',15*fact_axis)
hold on
xlim([0 6])
ylim([0 0.25])
xlabel('Increase in $R_{\mbox{\mbox{crit}}}^{H}$','interpreter','latex','FontSize',15*fact_label)

ylabel('Test avoidance $\varphi$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';
