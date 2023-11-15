function x_prob = beta_sampler(N,x,betap,p,norms)

j = 0:N;
alpha = j + betap; 
betad = (N - j) + betap;
x_prob = binopdf(j,N,p)*(betapdf(x,alpha,betad)'./norms);


