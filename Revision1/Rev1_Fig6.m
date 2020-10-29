%%% Figure 6 revised manuscript %%%

clear all
close all
clc

%% Data generation

str1 = 'Rev1_Fig6a.csv';
str2 = 'Rev1_Fig6b.csv';
str3 = 'Rev1_Fig6c.csv';

%% Par?metros solver

ti = [-100 0];
tf = [0 240]; tminplot = -80;
tmax = tf(end);
tmin = -80;
dt = 0.005;
factRtcrit = 0.95;
loc = 'northwest';
Phi = 15;
Rt2 = 2;

%% Definici?n de variables

xi = 0.15; 
varphi = 0.2;
xi = xi + (1-xi)*varphi; xim = 1-xi;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
eta = 0.66;
lambda_s = 0.1;
lambda_r = 0;
nmax = 300;

Rtcrit = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon)),2);
Rtcrit0 = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s,lambda_r,0,epsilon)),2);
Rt = factRtcrit*Rtcrit;

%% Calculo de as?ntotas

mode1 = maxVp(Pools(xi,nu,Rt2,Gamma,lambda_s,lambda_r,eta,epsilon));
mode2 = maxVp(Pools(xi,nu,Rt2,Gamma,lambda_s,lambda_r,0,epsilon));

%% Construcci?n de la adivinanza inicial en equilibrio

[Te,He,Hse,~,~,~] = equilibrio(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax);
x0 = [Te;He;Hse];
[~,~,~,Neq,ne,Neqcrit] = equilibrio_nmax(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax);
    
%% Problema directo
N_hat = [];
N_sum = [];
t = ti(1):dt:tf(end);
X = [];
rt = [Rt 2];
for j = 1:length(ti)
    if j>1
        x0 = X(:,end)'; X(:,end)=[];
    end
    x = solver_por_tramos(ti(j),tf(j),dt,xi,nu,rt(j),Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi,x0);
    if j>1
        N_hat = [N_hat(1:end-1) ; newInfections(ti(j):dt:tf(j),x,Gamma,nu,rt(j),lambda_s,lambda_r,nmax,eta)];
        N_sum = [N_sum(1:end-1) ; newInfections_Total(ti(j):dt:tf(j),x,Gamma,nu,epsilon,rt(j),Phi)];
    else
        N_hat = newInfections(ti(j):dt:tf(j),x,Gamma,nu,rt(j),lambda_s,lambda_r,nmax,eta);
        N_sum = newInfections_Total(ti(j):dt:tf(j),x,Gamma,nu,epsilon,rt(j),Phi);
    end
    X = [X x];
end
x = X; 

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
fact_axis = 2;
fact_label = 3;
fact_curva = 2;
siz = 15;
ylimsupI = 5e4;
ylimsupN = 5e3;
factRtcrit = 0.95;
Rtlim0 = 0.9;
Rtlim = 1.3;

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
xlim([tmin tmax])
ax.TickLabelInterpreter='latex';
hold on
ylim([0 ylimsupI])
%extractData(tmin,tmax,0,ylimsupI,str1)


%% Fig b

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

plot(T,N_hat_obs,lt,'Color',or1,'LineWidth',3*fact_curva)
hold on
plot(T,N_sum,ls,'Color',or2,'LineWidth',3*fact_curva)
hold on
% plot(T,Nsumobs,lt,'Color',or3,'LineWidth',3*fact_curva)
% hold on

plot([tmin tmax],[Neqcrit Neqcrit],'k--','LineWidth',3*fact_curva,'HandleVisibility','off')
legend({'$\hat{N}^{\mbox{obs}}$','$N$'},'interpreter','latex','FontSize',15*fact_axis);
set(gca,'FontSize',15*fact_axis)
hold on
xlabel('days','interpreter','latex','FontSize',15*fact_label)
ylabel('New infections','interpreter','latex','FontSize',15*fact_label)
xlim([tmin tmax])
hold on
ylim([0 ylimsupN])
ax.TickLabelInterpreter='latex';
%extractData(tmin,tmax,0,ylimsupN,str2)

%% Calculo del Rt

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';


plot(T,Rt_hat_obs,lt,'Color',red1,'LineWidth',3*fact_curva)
hold on
plot(T,Rt_hat,ls,'Color',red2,'LineWidth',3*fact_curva)
hold on
hold on
plot(T,ones(size(T)),'r--','LineWidth',1*fact_curva);
hold on
plot(T,exp(tau*mode1)*ones(size(T)),'k:','LineWidth',2*fact_curva);
hold on
plot(T,exp(tau*mode2)*ones(size(T)),'k--','LineWidth',2*fact_curva);
hold on
% legend({'$\hat{R}_t^{\mbox{obs}}$','$\hat{R}_t$'},'interpreter','latex','FontSize',15*fact_axis);
set(gca,'FontSize',15*fact_axis)
hold on
xlabel('days','interpreter','latex','FontSize',15*fact_label)
ylabel('Rep. Number','interpreter','latex','FontSize',15*fact_label)
xlim([tmin tmax])
hold on
ylim([0.9 Rtlim])
%extractData(tmin,tmax,0.9,Rtlim,str3)

ax.TickLabelInterpreter='latex';
