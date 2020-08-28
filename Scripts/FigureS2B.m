%%% An?lisis de sensibilidad par?metros no-controlados sobre RHcrit %%%

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

%% calculo 10%
RHcrit = zeros(N_sampling,1);
RHcrit0 = RHcrit;
Xireal =  betarnd(axi,bxi,N_sampling,1);
Varphi =  betarnd(avarxi,bvarxi,N_sampling,1);
Eta =  betarnd(aeta,beta,N_sampling,1);
Lambda_s = betarnd(als,bls,N_sampling,1);
Nu =  betarnd(anu,bnu,N_sampling,1);
Epsilon =  betarnd(aepsilon,bepsilon,N_sampling,1);
Xi = Xireal + (1-Xireal).*Varphi; 
    
for i = 1:N_sampling
    try
    RHcrit(i) = fzero(@(Rt) maxVp(Pools(Xi(i),Nu(i),Rt,Gamma,Lambda_s(i),lambda_r,Eta(i),Epsilon(i))),Rtcmax);
    catch err
        RHcrit(i) = NaN;
    end
    try
    RHcrit0(i) = fzero(@(Rt) maxVp(Pools(Xi(i),Nu(i),Rt,Gamma,Lambda_s(i),lambda_r,0,Epsilon(i))),Rtcmax);
    catch err
        RHcrit0(i) = NaN;
    end
end


load('DefColors.mat')

fact_axis = 1.2;
fact_label = 1.3;
fact_curva = 3;
siz = 15;
W = 8; H = 6;
BW = 0.01;
BWparam = 0.001;

%% Dist conjunta total
% histogram(RHcrit0(idx),'Normalization','pdf','FaceColor',C,'EdgeColor','k','FaceAlpha',0.25,'BinWidth',BW)
% set(gca,'FontSize',siz*fact_axis)
% histogram(RHcrit(idx),'Normalization','pdf','FaceColor',C,'EdgeColor','k','FaceAlpha',0.5,'BinWidth',BW)
figure('units','centimeters','position',[5 5 2*W 2*H]);
    ax = subplot(1,1,1);
    ax.Position = [0.28 0.3 0.65 0.55];
    ax.ActivePositionProperty = 'position';
C = Default(7,:);

idx = RHcrit0>0;
set(gca,'FontSize',siz*fact_axis)
hold on
pdi = fitdist(RHcrit0(idx),'gev');
xi = 1:0.01:4;
yi = pdf(pdi,xi); yi = yi/100;
p1 = plot(xi,yi,'LineWidth',fact_curva,'Color',C);
hold on
patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
hold on
idx = RHcrit>0;
pdi = fitdist(RHcrit(idx),'gev');
xi = 1:0.01:5;
yi = pdf(pdi,xi);  yi = yi/100;
hold on
p2 = plot(xi,yi,'LineWidth',fact_curva,'Color',C);
hold on
patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.75,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
hold on
p1.Color(4) = 0.25;
p2.Color(4) = 0.75;
set(gca, 'XTick', [1.0 1.5 2 2.5 3])
set(gca, 'YTick', [0 0.04])
xlim([0.75 3])
ylim([0 0.05])
xlabel('$R^H_{\mbox{crit}}$','interpreter','latex','FontSize',siz*fact_label)
ylabel('probability','interpreter','latex','FontSize',siz*fact_label)
legend({'$\eta=0$','$\eta=0.66$'},'interpreter','latex','FontSize',siz*fact_label);
ax.TickLabelInterpreter='latex';
%% Distribuciones individuales

%%%%%%%%

M = [Lambda_s Xi Eta Epsilon Nu Varphi];
str = {'$\lambda_s$','$\xi^{\mbox{ap}}$','$\eta$','$\epsilon$','$\nu$','$\varphi$'};
name = {'ls','xiap','eta','ep','nu'};

for i = [1 2 3 4 5]
    %     figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
    figure('units','centimeters','position',[5 5 W H]);
    ax = subplot(1,1,1);
    ax.Position = [0.28 0.3 0.65 0.55];
    ax.ActivePositionProperty = 'position';
    C = Default(i,:);
    v = M(:,i);
%     histogram(v,'Normalization','pdf','FaceColor',C,'EdgeColor','none','FaceAlpha',0.5,'BinWidth',BWparam,'LineWidth',0.1,'HandleVisibility','off')
    pdi = fitdist(v,'beta');
    xi = 0:0.01:1;
    yi = pdf(pdi,xi); yi = yi/100;
    hold on
    plot(xi,yi,'LineWidth',fact_curva,'Color',C)
    hold on
    patch([xi fliplr(xi)],[yi zeros(size(yi))],C,'FaceAlpha',0.5,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    set(gca,'FontSize',siz*fact_axis)
    hold on
    legend({str{i}},'interpreter','latex','FontSize',siz*fact_label);%,'Orientation','horizontal')
    
    xlim([0 1])
    ylim([0 0.15])
    set(gca, 'YTick', [0 0.1])

    xlabel(str{i},'interpreter','latex','FontSize',siz*fact_label)
    ylabel('probability','interpreter','latex','FontSize',siz*fact_label)
    ax.TickLabelInterpreter='latex';
end