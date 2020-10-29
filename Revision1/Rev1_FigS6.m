
clear all
close all
clc

%% Definici?n de parametros "Default"

xi = 0.15; xi0 = xi;
varphi = 0.2;
% xi = xi + (1-xi)*varphi; xim = 1-xi;
xi = 0.32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_s = 0.1;
lambda_r = 0;
etaesc = 0.66;
Phi = 1;
lmax = 1;
etamax = 1;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
eta = 0:etamax/100:etamax;

%% Definici?n par?metros barrido
kmax = length(eta);
lsmax = 1;
etamax = 1;
red = [1 0 0]; blue = [0 0 1]; green = [0 1 0];
fact_axis = 2;
fact_label = 2.5;
fact_curva = 2;
siz = 15;
eta = etaesc;
Rtcmax = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,1)),2); %eta


%% calculo 10%

Nmax = 100;  RT1 = 1; RT2 = 3;
Chi = linspace(0,2,Nmax);
Rtcrit = NaN(size(Chi));
for i = 1:Nmax
    try
    Rtcrit(i) = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Chi(i))),Rtcmax);
    catch err
        Rtcrit(i) = NaN;
    end
end


load('DefColors.mat')

%% Ploteo

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.15 0.25 0.65 0.65/2];
ax.ActivePositionProperty = 'position';

plot(Chi,Rtcrit,'Color',Default(1,:),'LineWidth',3*fact_curva)
hold on
plot(Chi,Rtcmax*ones(size(Chi)),'r--','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
plot(Chi,3.3*ones(size(Chi)),'r--','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
set(gca,'FontSize',15*fact_axis)
ylabel('$R^{H}_{\rm crit}$','interpreter','latex','FontSize',15*fact_label)
xlim([Chi(1) Chi(end)])
ylim([0.5 3.5])
%extractData(Chi(1),Chi(end),0.5,3.5,'Rev1_FigS6.csv')
xlabel('$\chi$','interpreter','latex','FontSize',15*fact_label)
title(strcat('$\xi^{\rm ap}=$',num2str(xi),', $R^{H}_{\rm crit}=$',num2str(round(Rtcmax,2))),'interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';
