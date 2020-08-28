%%% An?lisis de sensibilidad (pairwise) par?metros no-controlados sobre RHcrit %%%

clear all
close all
clc

%% Definici?n de parametros "Default"

lambda_s = 0.1;
lambda_r = 0;
Gamma = 0.1;
N_sampling = 1e5;

%% Definici?n par?metros para ser barridos

[axi,bxi] = findBparam(0.15,0.0025);           % mean = 0.15, variance pm 0.0025;
[avarxi,bvarxi] = findBparam(0.2,0.005);       % mean = 0.2, variance = 0.005;
[aeta,beta] = findBparam(0.66,0.005/4);        % mean = 0.66, variance = 0.0025;
[anu,bnu] = findBparam(0.1,0.005/2);             % mean = 0.1, variance = 0.005;
[aepsilon,bepsilon] = findBparam(0.1,0.005/2);   % mean = 0.1, variance = 0.005;
[als,bls] = findBparam(0.1,0.005/6);   % mean = 0.1, variance = 0.005;

%% Definici?n par?metros barrido

etadef = 0.66;
nudef = 0.1;
epsilondef = 0.1;
xi = 0.15; varphi = 0.2;
xidef = xi + (1-xi)*varphi; 


Rtcmax = fzero(@(Rtk) maxVp(Pools(xidef,nudef,Rtk,Gamma,lambda_s,lambda_r,etadef,epsilondef)),2); %eta

%% calculo

RHcrit = zeros(N_sampling,1);
Xireal =  betarnd(axi,bxi,N_sampling,1);
Varphi =  betarnd(avarxi,bvarxi,N_sampling,1);
Eta =  betarnd(aeta,beta,N_sampling,1);
Lambda_s = betarnd(als,bls,N_sampling,1);
Nu =  betarnd(anu,bnu,N_sampling,1);
Epsilon =  betarnd(aepsilon,bepsilon,N_sampling,1);
Xi = Xireal + (1-Xireal).*Varphi; 
M = [Lambda_s Xi Eta Epsilon Nu];
Mdef = [lambda_s*ones(size(Lambda_s)) xidef*ones(size(Lambda_s)) etadef*ones(size(Lambda_s)) epsilondef*ones(size(Lambda_s)) nudef*ones(size(Lambda_s))];


load('DefColors.mat')
red = [1 0 0]; blue = [0 0 1]; green = [0 1 0];

str = {'$\lambda_s$','$\xi^{\mbox{ap}}$','$\eta$','$\epsilon$','$\nu$','$\varphi$'};
str2 = {'Sympt-based testing $\lambda_s$','Ap. asympt. ratio $\xi^{\mbox{ap}}$','Tracing eff. $\eta$','Leak ratio $\epsilon$','Isol. ratio $\nu$','$\varphi$'};

RHcrit = zeros(N_sampling,5);
RHcrit0 = zeros(N_sampling,5);
RHpercentiles = zeros(5,3);
for j = [1 2 3 4 5]
    Param = Mdef; Param(:,j) = M(:,j);
    for i = 1:N_sampling
        try
            RHcrit(i,j) = fzero(@(Rt) maxVp(Pools(Param(i,2),Param(i,5),Rt,Gamma,Param(i,1),lambda_r,Param(i,3),Param(i,4))),Rtcmax);
            RHcrit0(i,j) = fzero(@(Rt) maxVp(Pools(Param(i,2),Param(i,5),Rt,Gamma,Param(i,1),lambda_r,0,Param(i,4))),Rtcmax);
        catch err
            RHcrit(i,j) = NaN;
            RHcrit0(i,j) = NaN;
        end
    end
    
    RHpercentiles(j,:) = [mean(RHcrit(:,j)) prctile(RHcrit(:,j),2.5) prctile(RHcrit(:,j),97.5)];
end

BW = 0.025;
BWparam = 0.025;
fact_axis = 1.2;
fact_label = 1.3;
fact_curva = 3;
siz = 15;
W = 8; H = 6;
  name = {'Rcrit_ls','Rcrit_xiap','Rcrit_eta','Rcrit_ep','Rcrit_nu'};

for j = [1 2 3 4 5]
    figure('units','centimeters','position',[5 5 W H]);
    ax = subplot(1,1,1);
    ax.Position = [0.28 0.32 0.65 0.55];
    ax.ActivePositionProperty = 'position';
    C = Default(j,:);
    if not(j==3)
        idx = RHcrit0(:,j)>0;
        set(gca,'FontSize',siz*fact_axis)
        hold on
        pdi = fitdist(RHcrit0(idx,j),'gev');
        xi = 1:0.01:4;
        yi = pdf(pdi,xi); yi = yi/100;
        p1 = plot(xi,yi,'LineWidth',fact_curva,'Color',C);
        hold on
        patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
        hold on
        idx = RHcrit(:,j)>0;
        pdi = fitdist(RHcrit(idx,j),'gev');
        xi = 1:0.01:5;
        yi = pdf(pdi,xi);  yi = yi/100;
        hold on
        p2 = plot(xi,yi,'LineWidth',fact_curva,'Color',C);
        hold on
        patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.75,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
        hold on
        p1.Color(4) = 0.25;
        p2.Color(4) = 0.75;
        ax.TickLabelInterpreter='latex';
    else
        idx = RHcrit(:,j)>0;
        pdi = fitdist(RHcrit(idx,j),'gev');
        xi = 1:0.01:5;
        yi = pdf(pdi,xi);  yi = yi/100;
        hold on
        p2 = plot(xi,yi,'LineWidth',fact_curva,'Color',C);
        hold on
        patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.75,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
        hold on
        p1.Color(4) = 0.25;
        p2.Color(4) = 0.75;
        ax.TickLabelInterpreter='latex';
    end
    set(gca,'FontSize',15*fact_axis)
    hold on
    xlim([0.75 3])
    ylim([0 0.15])
            set(gca, 'XTick', [1.0 1.5 2 2.5 3])
        set(gca, 'YTick', [0 0.1])
    xlabel('$R^H_{\mbox{crit}}$','interpreter','latex','FontSize',15*fact_label)
    ax.TickLabelInterpreter='latex';
end