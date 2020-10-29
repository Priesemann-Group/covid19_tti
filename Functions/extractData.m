function extractData(xi,xf,yi,yf,name)
h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
[n,~] = size(x);
if n == 1
    x = {x}; y = {y};
end
N = length(x{1});
Labels = cell(1,2*n);
Data = NaN(N,2*n);
for i = 1:n
    Labels{2*i-1} = strcat('xline',num2str(i));
    Labels{2*i} = strcat('yline',num2str(i));
    X = x{i}'; Y = y{i}';
    idx1 = X<xi | X>xf; X(idx1) = NaN;
    idx2 = Y<yi | Y>yf; Y(max(idx1,idx2)) = NaN;
    Data(:,2*i-1) = X;
    Data(:,2*i) = Y;
end
idxf = not(isnan(Data(:,1))); Data = Data(idxf,:);
T = array2table(Data,'VariableNames',Labels);
writetable(T,name)