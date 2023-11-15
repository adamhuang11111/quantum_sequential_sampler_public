function norms = norm_beta(N,betap,x_all)


norms = zeros(N+1,1);

for j = 0:N
    alpha = j + betap; 
    betad = (N - j) + betap;
    norms(j+1,1) = sum(betapdf(x_all,alpha,betad));
end
