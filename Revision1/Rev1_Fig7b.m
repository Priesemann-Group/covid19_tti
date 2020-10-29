%%% main level plot paper %%%

clear all
close all
clc

%% Definici?n de parametros "fijos"

xi = 0.15;
varphi = 0.20;
Rt = [3 2.5 2 1.5];
lmax = 0.08;
etamax = 1;
xi = xi + (1-xi)*varphi; xim = 1-xi;
Gamma = 0.1;
nu = 0.1;
epsilon = 0.1;
lambda_r = 0:lmax/100:lmax;
eta = 0;
%% Parametros ploteo

red = [1 0 0]; blue = [0 0 1];
fact_axis = 2;
fact_label = 3;
fact_curva = 2;
siz = 15;

%% Definici?n par?metros barrido
kmax = length(Rt);
lambda_s = NaN(kmax,length(lambda_r));
str = cell(1,kmax);

for k = 1:kmax
    for i = 1:length(lambda_r)
        try
            lambda_s(k,i) = fzero(@(lambda_s) maxVp(Pools(xi,nu,Rt(k),Gamma,lambda_s,lambda_r(i),eta,epsilon)),0);
        catch err
            disp('oops')
        end
    end
end


figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
ax = subplot(1,1,1);
ax.Position = [0.25 0.25 0.65 0.65];
ax.ActivePositionProperty = 'position';



plot([0 1],[0.002 0.002],'r:','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
plot([0 1],[0.02 0.02],'r--','LineWidth',2*fact_curva,'HandleVisibility','off');
hold on
for i = 1:kmax
    a = 0.5*(i-1)/(kmax-1);
    plot(lambda_s(i,:),lambda_r,'Color',(1-a)*red+a*blue,'LineWidth',3*fact_curva);
    str{1,i} = char(strcat("$R^H\!=",num2str(Rt(i)),"$"));
    hold on
end
ylim([0 lmax])
xlim([0 1])
%extractData(0,lmax,0,1,'Rev1_Fig7b.csv')
annotation('textbox',[.75116 .79004 .13923 .11438],'String','Controlled','interpreter','latex','FontSize',10*fact_axis,'Linewidth',2,'FitBoxToText','on','BackgroundColor','w');
annotation('textbox',[.2172 .27628 .15648 .081522],'String','Uncontrolled','interpreter','latex','FontSize',10*fact_axis,'Linewidth',2,'FitBoxToText','on','BackgroundColor','w');
legend(str,'interpreter','latex','FontSize',siz*fact_axis,'Location','northeast');


set(gca,'FontSize',siz*fact_axis)
hold on
ylabel('$\lambda_r$','interpreter','latex','FontSize',siz*fact_label)
xlabel('$\lambda_s$','interpreter','latex','FontSize',siz*fact_label)
ax.TickLabelInterpreter='latex';