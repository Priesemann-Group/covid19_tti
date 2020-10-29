%%% Fig S2 a %%%

clear all
close all
clc

%% Definici?n de parametros "Default"

xi = 0.15; 
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
eta = etaesc;
Rtcmax = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,1)),2); %eta


%% Definicion de rangos para los parametros
Nmax = 100;
% Chi = linspace(0,2,Nmax);
nu0 = 0; nu1 = 1;
epsilon0 = 0; epsilon1=1;
ls0 = 0; ls1 = 1;
eta0 = 0 ; eta1 = 1;
xi0 = 0; xi1=1;
Nu = linspace(nu0,nu1,Nmax);
Epsilon = linspace(epsilon0,epsilon1,Nmax);
Lambda_s = linspace(ls0,ls1,Nmax);
Eta = linspace(eta0,eta1,Nmax);
Xi = linspace(xi0,xi1,Nmax);
Rtcrit = NaN(5,Nmax);
for i = 1:Nmax
    try
    Rtcrit(1,i) = fzero(@(Rt) maxVp(Pools_forReview(xi,Nu(i),Rt,Gamma,lambda_s,lambda_r,eta,epsilon,1)),Rtcmax);
    catch err
        Rtcrit(1,i) = NaN;
    end
    try
    Rtcrit(2,i) = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,Epsilon(i),1)),Rtcmax);
    catch err
        Rtcrit(2,i) = NaN;
    end
    try
    Rtcrit(3,i) = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,Lambda_s(i),lambda_r,eta,epsilon,1)),Rtcmax);
    catch err
        Rtcrit(3,i) = NaN;
    end
    try
    Rtcrit(4,i) = fzero(@(Rt) maxVp(Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,Eta(i),epsilon,1)),Rtcmax);
    catch err
        Rtcrit(4,i) = NaN;
    end
    try
    Rtcrit(5,i) = fzero(@(Rt) maxVp(Pools_forReview(Xi(i),nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,1)),Rtcmax);
    catch err
        Rtcrit(5,i) = NaN;
    end
end


load('DefColors.mat')

%% Ploteo

M = [Nu ; Epsilon ; Lambda_s ; Eta ; Xi];
xl = {'$\nu$'; '$\epsilon$'; '$\lambda_s$' ; '$\eta$' ; '$\xi^{\rm ap}$'};
num = {'1';'2';'3';'5';'4'};
W = 8; H = 6;

for i = 1:5
    figure('units','centimeters','position',[5 5 W H]);
    ax = subplot(1,1,1);
    ax.Position = [0.28 0.32 0.65 0.55];
    ax.ActivePositionProperty = 'position';
    
    plot(M(i,:),Rtcrit(i,:),'Color',Default(1,:),'LineWidth',3*fact_curva)
    hold on
    plot(M(i,:),Rtcmax*ones(size(M(i,:))),'r--','LineWidth',2*fact_curva,'HandleVisibility','off');
    hold on
    set(gca,'FontSize',15*fact_axis)
    ylabel('$R^{H}_{\rm crit}$','interpreter','latex','FontSize',15*fact_label)
    xlim([0 1])
    ylim([0.5 3.5])
    xlabel(xl{i},'interpreter','latex','FontSize',15*fact_label)
    ax.TickLabelInterpreter='latex';
    %extractData(0,1,0.5,3.5,strcat('Rev1_FigS2a_',num{i},'.csv'))
end
