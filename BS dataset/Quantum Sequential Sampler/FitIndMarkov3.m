function [nLL, Pred,Psyfs,Ms] = FitIndMarkov3(parm,Sdat,cs)

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
o1 = 0;
%build probability model
pA = parm(1); pnA = 1 - pA;
pB = parm(2); pnB = 1 - pB;
pBgA = parm(3); pnBgA = 1 - pBgA;

o2 = (-(pA>=pB) + (pA<pB)).*o1;
o3 = o1;

%% A & B
%%%%%%%%%%%%%%%%
% known conjunctions
pAB =   pA*pBgA ;
pAnB =  pA*pnBgA;

% order effects
pBA = pAB - o1;
pBA = (pBA >= pB)*pB + (pBA <= 0)*0 + (pBA > 0)*(pBA < pB)*pBA;
pnBA = pAnB - o2;
pnBA = (pnBA >= pnB)*pnB + (pnBA <= 0)*0 + (pnBA > 0)*(pnBA < pnB)*pnBA;
pnBnA = pnB - pnBA;
pBnA = pB - pBA;

% reversed order effects
pnAnB = pnBnA - o3;
pnAnB = (pnAnB >= pnA)*pnA + (pnAnB <= 0)* 0 + (pnAnB > 0)*(pnAnB < pnA)*pnAnB;
pnAB = pnA - pnAnB; 

% conditionals
pAgB = pBA/pB;      pnAgB = 1-pAgB;
pAgnB = pnBA/pnB ; pnAgnB = 1-pAgnB ;

% reversed conditionals
pnBgnA = pnAnB/pnA; pBgnA = 1 - pnBgnA;
%% AB pairs
if pA > pB
    pAandB = pAB;
    pnAandnB = pnBnA;    
else % pB > pA
    pAandB = pBA;
    pnAandnB = pnAnB;   
end

if pA > pnB
    pAandnB = pAnB;
    pnAandB = pBnA;
else
    % pnB > pA
    pAandnB = pnBA;
    pnAandB = pnAB;
end

pAorB = 1 - pnAandnB;
pnAorB = 1 - pAandnB;
pAornB = 1 - pnAorB;
pnAornB = 1 - pAandB;

%%
Pred = [pA, pB, pnA, pnB, pAandB, pAandnB, pnAandB, pnAandnB, pAorB, pAornB, pnAorB, ...
    pnAornB, pAgB, pAgnB, pnAgB, pnAgnB, pBgA, pBgnA, pnBgA, pnBgnA]';

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