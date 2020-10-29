%%% Figure 3a %%%

clear all
close all
clc

%% Temporal

ti = [-35 0];
tf = [0 360];
tmax = tf(end);
tmin = -20;
nmax = 300; 
dt = 0.005;

%% Variables iniciales (default)
% 
xi= 0.15; 
varphi = 0.2;
eta_sc = 0.66;
eta = [0 eta_sc];
lambda_s = [0 0.1];
ylimsupI = 250;
ylimsupN = 30;
loc = 'northeast';
xi = xi + (1-xi)*varphi; xim = 1-xi;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
lambda_r = 0;
Phi = 15;

nrt = 250;
Rtcrit = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s(end),lambda_r,eta(end),epsilon)),2);
Rtcrit0 = fzero(@(Rt) maxVp(Pools(xi,nu,Rt,Gamma,lambda_s(end),lambda_r,0,epsilon)),2);
RT = linspace(0,Rtcrit,nrt);

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

NeqRTH = NaN(size(RT));
NeqRT0 = NaN(size(RT));

for k = 1:kmax
    rt = [RT(k) RT(k)];
    X = [];
    if RT(k)<Rtcrit
        [Te,He,Hse,Neq,ne,Neqcrit] = equilibrio(xi,nu,RT(k),Gamma,lambda_s(end),lambda_r,eta(end),epsilon,Phi,nmax);
        NeqRTH(k) = Neq;
    end
    if RT(k) <Rtcrit0
        [Te,He,Hse,Neq,ne,Neqcrit] = equilibrio(xi,nu,RT(k),Gamma,lambda_s(end),lambda_r,0,epsilon,Phi,nmax);
        NeqRT0(k) = Neq;
    end
end
rtcrit = linspace(Rtcrit0,Rtcrit,100);
neqcrit = NaN(size(rtcrit));
for k = 1:100
    rt = [rtcrit(k) rtcrit(k)];
    X = [];
    if rtcrit(k)<Rtcrit
        [Te,He,Hse,Neq,ne,Neqcrit] = equilibrio(xi,nu,rtcrit(k),Gamma,lambda_s(end),lambda_r,eta(end),epsilon,Phi,nmax);
        neqcrit(k) = Neqcrit;
    end
end

%% Visualization

fact_axis = 2;
fact_label = 3;
fact_curva = 2;
siz = 15;

green = [0 1 0]; 
blue = [0 0 1];

figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';


maxim = max(NeqRTH);
hold on
plot([1 1],[0 maxim],'-.','LineWidth',3*fact_curva,'Color',[0 0.5 0],'HandleVisibility','off');
hold on
plot([Rtcrit0 Rtcrit0],[0 maxim],'-.','LineWidth',3*fact_curva,'Color','y','HandleVisibility','off');
hold on
plot([Rtcrit Rtcrit],[0 maxim],'-.','LineWidth',3*fact_curva,'Color',[0.6 0 0],'HandleVisibility','off');
hold on
plot(RT,NeqRTH,'LineWidth',3*fact_curva);
hold on
plot(RT,NeqRT0,'LineWidth',3*fact_curva);
hold on
plot(rtcrit,neqcrit,'LineWidth',3*fact_curva,'HandleVisibility','off');
hold on

legend({'$\hat{N}_{\mbox{eq}}^{\mbox{obs}}$($\eta\!=\!0.66$)',...
    '$\hat{N}_{\mbox{eq}}^{\mbox{obs}}$($\eta\!=\!0$)',...
    '$N_{\mbox{max}}$'},'interpreter','latex','FontSize',15*fact_axis,'Location','southeast');

set(gca,'FontSize',15*fact_axis)
hold on
xlabel('$R^{H}_t$','interpreter','latex','FontSize',15*fact_label)
ylabel('$\hat{N}_{\mbox{eq}}^{\mbox{obs}}$','interpreter','latex','FontSize',15*fact_label)
hold on
ylim([0 1000])
xlim([0.5 3.5])
%extractData(0.5,3.5,0,1000,'Rev1_FigS5b.csv')
ax.TickLabelInterpreter='latex';
