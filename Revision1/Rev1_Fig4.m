%%% Resoluci?n del problema directo con din?mica saturable %%%

clear all
close all
clc

%% Choose tact

tact = [0 45 90 135 180];
% tact = [0 90 180];
str1 = 'Rev1_Fig4a.csv';
str2 = 'Rev1_Fig4b.csv';
str3 = 'Rev1_Fig4c.csv';
%% Par?metros solver

tmin = -80; tminplot = -20;
tmax = 240;
dt = 0.005;
t = tmin:dt:tmax;
ylimsupI = 3e4;
ylimsupN = 4e3;
Rtlim = 2.5;
factRtcrit = 0.95;
loc = 'northeast';
Phimax = 3331;          % PhiNmax
sigma = 2;
Phi = 15;
nmax = 300;
eta = 0.66;

%% Definici?n de variables

xi = 0.15; 
varphi = 0.2;
xi = xi + (1-xi)*varphi; xim = 1-xi;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
lambda_s = 0.1;
lambda_r = 0;

Rtcrit = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
Rtcrit0 = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s,lambda_r,0,epsilon)),2);
Rt = factRtcrit*Rtcrit;
%% Visualizaci?n

fact_axis = 2;
fact_label = 3;
fact_curva = 2;
siz = 15;

%% Construcci?n de la adivinanza inicial en equilibrio

[Te,He,Hse,~,~,~] = equilibrio(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax);
x0 = [Te;He;Hse];
[~,~,~,Neq,ne,Neqcrit] = equilibrio_nmax(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax);

%% Problema directo

[x] = ode4(@(t,x) Pools_solver_phit(t,x,xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi,Phimax,sigma,tact),t,x0);
x = x';
N_hat = newInfections(t,x,Gamma,nu,Rt,lambda_s,lambda_r,nmax,eta);
N_sum = newInfections_Total(t,x,Gamma,nu,epsilon,Rt,PHI(t,tact,sigma,Phi,Phimax));

%% Discretizacion

idx = t==floor(t);
T = t(idx);
N_hat = N_hat(idx);
N_sum = N_sum(idx);

%% Calculo de los casos nuevos por d?a

N_hat_obs = EstimDelay(N_hat,4,1,0.95);
tau = 4;
Rt_hat_obs = N_hat_obs./[NaN(tau,1) ; N_hat_obs(1:end-tau)];
tau = 4;
N_T = EstimDelay(N_sum,4,1,0.95);
Rt_hat = N_T./[NaN(tau,1) ; N_T(1:end-tau)];

%% Visualizaci?n

load('Colores.mat')
lt = '-';
ls = '-.';
lh = ':';

%% Fig a

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

plot(t,x(1,:),lt,'Color',bl1,'LineWidth',3*fact_curva)
hold on
plot(t,x(1,:) + x(2,:),ls,'Color',bl2,'LineWidth',3*fact_curva)
hold on
plot(t,x(2,:),lh,'Color',bl3,'LineWidth',3*fact_curva)
hold on
legend({'Traced $T$','Total $T+H$','Hidden $H$'},'interpreter','latex','FontSize',15*fact_axis,'Location',loc);
set(gca,'FontSize',15*fact_axis)
hold on
xlabel('days','interpreter','latex','FontSize',15*fact_label)
ylabel('Active cases','interpreter','latex','FontSize',15*fact_label)
xlim([tminplot tmax])
ax.TickLabelInterpreter='latex';
ylim([0 ylimsupI])
%extractData(tminplot,tmax,0,ylimsupI,str1)
hold on

%% Fig b

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

plot(T,N_hat_obs,lt,'Color',or1,'LineWidth',3*fact_curva)
hold on
plot(T,N_sum,ls,'Color',or2,'LineWidth',3*fact_curva)
hold on

Phi = PHI(t,tact,sigma,Phi,Phimax);
p1 = plot(T,Phi(idx),'k--','LineWidth',3*fact_curva);
p1.Color(4) = 0.25;
hold on
plot(T,Neqcrit*ones(size(T)),'--','Color','k','LineWidth',2*fact_curva);
hold on
legend({'$\hat{N}^{\mbox{obs}}$','$N$','$\Phi(t)$'},'interpreter','latex','FontSize',15*fact_axis);
set(gca,'FontSize',15*fact_axis)
hold on
xlabel('days','interpreter','latex','FontSize',15*fact_label)
ylabel('New infections','interpreter','latex','FontSize',15*fact_label)
xlim([tminplot tmax])
hold on
ylim([0 ylimsupN])
%extractData(tminplot,tmax,0,ylimsupN,str2)
ax.TickLabelInterpreter='latex';

%% Calculo del Rt

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';


plot(T,Rt_hat_obs,lt,'Color',red1,'LineWidth',3*fact_curva)
hold on
plot(T,Rt_hat,ls,'Color',red2,'LineWidth',3*fact_curva)
hold on
plot(T,ones(size(T)),'k--','LineWidth',1*fact_curva);
hold on
legend({'$\hat{R}_t^{\mbox{obs}}$','$\hat{R}_t$'},'interpreter','latex','FontSize',15*fact_axis);
set(gca,'FontSize',15*fact_axis)
hold on
xlabel('days','interpreter','latex','FontSize',15*fact_label)
ylabel('Rep. Number','interpreter','latex','FontSize',15*fact_label)
xlim([tminplot tmax])
hold on
ylim([0.5 Rtlim])
hold on
set(gca,'FontSize',15*fact_axis)
%extractData(tminplot,tmax,0.5,Rtlim,str3)
ax.TickLabelInterpreter='latex';
