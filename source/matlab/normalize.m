function x = normalize(x)
x = x ./ sqrt(sum(x.^2,1));
