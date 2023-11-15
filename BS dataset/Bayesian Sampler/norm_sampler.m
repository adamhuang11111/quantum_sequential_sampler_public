function norms = norm_sampler(N,betap,x)


norms = zeros(N+1,1);

for j = 0:N
    alpha = j + betap; 
    betad = (N - j) + betap;
    norms(j+1,1) = sum(betapdf(x,alpha,betad));
end
