function Ys = mediamovil(Y,k)
n = length(Y);
Ys = zeros(size(Y));
for i = 1:n
    if i<k+1
        Ys(i) = mean(Y(1:i+k));
    elseif i<n-k
        Ys(i) = mean(Y(i-k:i+k));
    else
        Ys(i) = mean(Y(i-k:end));
    end
end