function [nLL, Pred, Ms] = FitIndBSBeta3(parm,Sdat,cs)

beta = parm(4);
N2 = 1+round(parm(6));
N1 = N2+round(parm(5));
x = [ .005 (.01:.01:.99) .995];

%build bayesian probability
pA = parm(1);
pnA = 1 - pA;
pBgA = parm(2);
pnBgA = 1 - pBgA;
pBgnA = parm(3);
pnBgnA = 1 - pBgnA;
pB = pA*pBgA + pnA*pBgnA;
pnB = 1 - pB;

pAB = pA * pBgA;
pnAB = pnA * pBgnA;
pAnB = pA * pnBgA;
pnAnB = pnA * pnBgnA;

pAgB = pAB/pB;
pnAgB = pnAB/pB;
pAgnB = pAnB/pnB;
pnAgnB = pnAnB/pnB;

Pred = [pA, pB, pnA, pnB, pAB, pAnB, pnAB, pnAnB, (1 - pnAnB), (1 - pnAB), (1 - pAnB), ...
    (1 - pAB), pAgB, pAgnB, pnAgB, pnAgnB, pBgA, pBgnA, pnBgA, pnBgnA]';

%build sampler
ns = size(Pred,1);
nLL = 0;
Ms = zeros(ns,1);
for k = 1:ns
    Pk = Pred(k,1);
    if 4 < k && k < 13
        N = N2;
    else
        N = N1;
    end
    for i = 1:3
        Rk = Sdat(k,i);
        norms = norm_sampler(N,beta,x);
        PR = eps;
        %cs = 5 refers to the rounding mechanism into 5s and 10s.
        if cs == 5
            for j = 0:(cs-1)
                Rkj = (Rk + j)/100 +.005;
                PR = PR +  beta_sampler(N, Rkj, beta, Pk, norms);
            end
        else 
            Rkj = (Rk==0).*(.005) + (Rk==100).*(.995) + (Rk>0).*(Rk<100).*(Rk/100);
            PR = beta_sampler(N, Rkj, beta, Pk, norms);
        end
        %log-likelihood
        nLL = nLL + log(PR);
    end
    ms = 0;
    for ij = 1:101
        ms = ms + beta_sampler(N, x(ij), beta, Pk, norms)* (ij-1);
    end
    Ms(k) = ms;
end
nLL = -2 * nLL;














