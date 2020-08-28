%%% script problema directo correcciones por tramos %%%

clear all
close all
clc

%% Temporal

ti = [-35 0];
tf = [0 360];
tmax = tf(end);
tmin = -20;
nmax = 300; 
dt = 0.1;

%% Variables iniciales (default)
% 
xi= 0.15; 
varphi = 0.2;
eta_sc = 0.66;
eta = [0 eta_sc];
lambda_s = [0 0.1];
xi = xi + (1-xi)*varphi; xim = 1-xi;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
lambda_r = 0;
Phi = 1;

nrt = 1500;
RT = linspace(0.5,4,nrt);
Rtcrit = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s(end),lambda_r,eta(end),epsilon)),2);
Rtcrit0 = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s(end),lambda_r,0,epsilon)),2);

%% Definici?n par?metros barrido

t = ti(1):dt:tf(end);

kmax = length(RT);
Xs = cell(kmax,1);
Rts = cell(kmax,1);
Rhatobs = cell(kmax,1);
Rhat = cell(kmax,1);
idx = t==floor(t);
T = t(idx);
larg_eigenv = zeros(size(RT));
larg_eigenv2 = zeros(size(RT));

for i = 1:kmax
    larg_eigenv(i) = maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s(end),lambda_r,eta(end),epsilon));
    larg_eigenv2(i) = maxVp(Pools(xi,nu,RT(i),Gamma,lambda_s(end),lambda_r,0,epsilon));
end


for k = 1:kmax
    rt = [RT(k) RT(k)];
    X = [];
    for j = 1:length(ti)
        if j>1
            x0 = X(:,end)'; X(:,end)=[];
        else
            x0 = [0 ; 1 ;0];
        end
        x = solver_por_tramos(ti(j),tf(j),dt,xi,nu,rt(j),Gamma,lambda_s(j),lambda_r,eta(j),epsilon,nmax,Phi,x0);
        if j>1
            N_hat = [N_hat(1:end-1) ; newInfections(ti(j):dt:tf(j),x,Gamma,nu,rt(j),lambda_s(j),lambda_r,nmax,eta(j))];
            N_sum = [N_sum(1:end-1) ; newInfections_Total(ti(j):dt:tf(j),x,Gamma,nu,epsilon,rt(j),Phi)];
        else
            N_hat = newInfections(ti(j):dt:tf(j),x,Gamma,nu,rt(j),lambda_s(j),lambda_r,nmax,eta(j));
            N_sum = newInfections_Total(ti(j):dt:tf(j),x,Gamma,nu,epsilon,rt(j),Phi);
        end
        X = [X x];
    end
    x = X;
    N_hat_obs = EstimDelay(N_hat,4,1,0.9);
    tau = 4;
    Rt_hat_obs = N_hat_obs./[NaN(tau,1) ; N_hat_obs(1:end-tau)];
    N_T = EstimDelay(N_sum,4,1,0.9);
    Rt_hat = N_T./[NaN(tau,1) ; N_T(1:end-tau)];
    Rhatobs{k} = Rt_hat_obs;
end

%% Visualizaci?n

fact_axis = 2;
fact_label = 3;
fact_curva = 2;
siz = 15;

green = [0 1 0]; 
blue = [0 0 1];

ylimsupI = 250;
ylimsupN = 30;
loc = 'northeast';


figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';

hold on
plot([1 1],[0 2],'-.','LineWidth',3*fact_curva,'Color',[0 0.5 0],'HandleVisibility','off');
hold on
plot([Rtcrit0 Rtcrit0],[0 2],'-.','LineWidth',3*fact_curva,'Color','y','HandleVisibility','off');
hold on
plot([Rtcrit Rtcrit],[0 2],'-.','LineWidth',3*fact_curva,'Color',[0.6 0 0],'HandleVisibility','off');
hold on
plot(RT,max(1,exp(tau*larg_eigenv)),'k:','LineWidth',2*fact_curva);%,'HandleVisibility','off');
hold on
plot(RT,max(1,exp(tau*larg_eigenv2)),'k--','LineWidth',2*fact_curva);%,'HandleVisibility','off');
hold on
plot([0 Rtcrit0],[1 1],'k-','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
legend({'$\eta\!=\!0.66$','$\eta\!=\!0$'},'interpreter','latex','FontSize',15*fact_axis,'Location','southeast');

set(gca,'FontSize',15*fact_axis)
hold on
xlabel('$R_t^{H}$','interpreter','latex','FontSize',15*fact_label)
ylabel('Rep. number','interpreter','latex','FontSize',15*fact_label)
hold on
ylim([0.5 1.7])
xlim([0.5 4])
ax.TickLabelInterpreter='latex';
