function [nLL, Pred,Psyfs,Ms] = FitIndMarkov3_classical(parm,Sdat,cs)

%Each repeated measure drift in the same rate but different time
alphakk = parm(4);
ckk = parm(5);

%Build Initial State
m = 101;
x = [ .005 (.01:.01:.99) .995];
alpha = parm(6);
beta = alpha;
x0 = betapdf(x,alpha,beta)';
Psy0 = x0./sum(x0);

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
nq = size(Sdat,2);

nLL = 0;
Ms = zeros(ns,1);
Psyfs = zeros(ns,101,3);

for k = 1:ns
    Pk = Pred(k,1);
    K = drift_add([alphakk,ckk,Pk],m);
    T1 = expm(K) ;
    Psyf = T1*Psy0;
    Ms(k) = (0:100)*Psyf;
%     Psyfs(k,:,:) = Psyf;
    for i = 1:nq
        Rk = Sdat(k,i);
        PR = eps;
        cc = cs;
        if (cs == 5) && (Rk >= 95)
            cc = 6;
        end
        for j = 0:(cc-1)
            Rkj = Rk + j;
            PR = PR + Psyf(Rkj+1);
        end
        %log-likelihood
        nLL = nLL + log(PR);
    end
end
nLL = -2 * nLL;