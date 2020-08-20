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

Rtcmax = fzero(@(Rtk) maxVp(Pools(xi,nu,Rtk,Gamma,lambda_s,lambda_r,etaesc,epsilon)),2); %eta
eta = etaesc;

%% calculo 10%

Nmax = 100;  RT1 = 1; RT2 = 3;
RT = linspace(RT1,RT2,Nmax);
ETA = zeros(Nmax,1)'; 
Lambda_s = zeros(Nmax,1)'; 
Lambda_r = zeros(Nmax,1)'; 
Epsilon = zeros(Nmax,1)'; 
Nu = zeros(Nmax,1)';
Varphi = zeros(Nmax,1)'; 
for i = 1:Nmax
    if ETA(i) <= etamax
        try
            ETA(i) = fzero(@(eta) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
        catch err
            ETA(i) = NaN;
        end
    else
        ETA(i) = NaN;
    end
    try
    Lambda_s(i) = fzero(@(lambda_s) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
    catch err
        Lambda_s(i) = NaN;
    end
    try
    Lambda_r(i) = fzero(@(lambda_r) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
    catch err
        Lambda_r(i) = NaN;
    end
    try
    Epsilon(i) = fzero(@(epsilon) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
    catch err
        Epsilon(i) = NaN;
    end
    try
    Nu(i) = fzero(@(nu) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
    catch err
        Nu(i) = NaN;
    end
    try
        xi10 = fzero(@(xi) maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s,lambda_r,eta,epsilon)),2);
        Varphi(i) = (xi10-xi0)/(1-xi0);
    catch err
        Varphi(i) = NaN;
    end
end


load('DefColors.mat')

%% Ploteo

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.15 0.25 0.65 0.65/2];
ax.ActivePositionProperty = 'position';

plot(RT,Lambda_s,'Color',Default(1,:),'LineWidth',3*fact_curva)
hold on
ylabel('$\lambda_s$','interpreter','latex','FontSize',15*fact_label)
yyaxis right
plot(RT,Lambda_r,'Color',Default(2,:),'LineWidth',3*fact_curva)
hold on
plot([RT1 RT2],[0.002 0.002],'r:','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
plot([RT1 RT2],[0.02 0.02],'r--','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
ylim([0 0.025])

set(gca,'FontSize',15*fact_axis)
ylabel('$\lambda_r$','interpreter','latex','FontSize',15*fact_label)
xlim([RT1 RT2])
set(gca, 'XTick', [1.0 1.5 2 2.5 3])
set(gca, 'XTickLabel', {'1.0' '1.5' '2.0' '2.5' '3.0'})
xlabel('$R^H$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';


figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.15 0.25 0.65 0.65/2];
ax.ActivePositionProperty = 'position';

plot(RT,ETA,'Color',Default(3,:),'LineWidth',3*fact_curva)
hold on
ylabel('$\eta$','interpreter','latex','FontSize',15*fact_label)
ylim([0 1])
set(gca,'FontSize',15*fact_axis)
xlim([RT1 RT2])

set(gca, 'XTick', [1.0 1.5 2 2.5 3])
set(gca, 'XTickLabel', {'1.0' '1.5' '2.0' '2.5' '3.0'})
xlabel('$R^H$','interpreter','latex','FontSize',15*fact_label)
ax.TickLabelInterpreter='latex';


figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

plot(RT,Epsilon,'Color',Default(4,:),'LineWidth',3*fact_curva)
hold on
plot(RT,Nu,'Color',Default(5,:),'LineWidth',3*fact_curva)
hold on
plot(RT,Varphi,'Color',Default(6,:),'LineWidth',3*fact_curva)
hold on
ylabel('External factors','interpreter','latex','FontSize',15*fact_label)
ylim([0 1])
set(gca,'FontSize',15*fact_axis)
xlim([RT1 RT2])

set(gca, 'XTick', [1.0 1.5 2 2.5 3])
set(gca, 'XTickLabel', {'1.0' '1.5' '2.0' '2.5' '3.0'})
xlabel('$R^H$','interpreter','latex','FontSize',15*fact_label)
legend({'$\epsilon$ ``leak'''' fract.','$\nu$ isol. factor','$\varphi$ test avoid.'},'interpreter','latex','FontSize',15*fact_label);%,'Orientation','horizontal')

ax.TickLabelInterpreter='latex';

