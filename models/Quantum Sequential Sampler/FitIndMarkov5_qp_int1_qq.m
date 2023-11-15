function [nLL, nLLc, Pred, Psy0,Ks,Ms,Ords] = FitIndMarkov5_qp_int1_qq(parm,Sdat,cs)
nd = size(Sdat,1);
alphakk = parm(7);
ckk = parm(8);
o1 = parm(10);
o3 = o1;
o4 = o1;
o6 = o1;
o7 = o1;
o9 = o1;

%% Build Initial State
m = 101;

%Gaussian initial state
% mid = (m+1)/2;
% wd = parm(9);
%
%
% x0 = ((1:m) - mid)';
% x0 = (x0./wd).^2;       % normal dist around mid with std = wd
% x0 = exp(-x0);
% Psy0 = x0./sum(x0);

% Beta Initial State
x = [ .005 (.01:.01:.99) .995];
alpha = parm(9);
beta = alpha;
x0 = betapdf(x,alpha,beta)';
Psy0 = x0./sum(x0);
%%
% true probabilities
% A = A1,  B = A2,  C = A3
pA = parm(1); pnA = 1 - pA;
pB = parm(2); pnB = 1 - pB;
pC = parm(3); pnC = 1 - pC;

%Classical models have these bounds in order to make sure that pAB is not
%greater than pB as a result of pA*pBgA can be fitted to be greater than
%pB.
o2 = (-(pA>=pB) + (pA<pB)).*o1;
o5 = (-(pA>=pC) + (pA<pC)).*o1;
o8 = (-(pB>=pC) + (pB<pC)).*o1;

pBgA = parm(4); pnBgA = 1 - pBgA;
pCgA = parm(5); pnCgA = 1 - pCgA;
pCgB = parm(6); pnCgB = 1 - pCgB;

%% A & B
%%%%%%%%%%%%%%%%
% known conjunctions
pAB =   pA*pBgA ;
pAnB =  pA*pnBgA;

% order effects
pBA = pAB - o1;
%Quantum constrains
pBA = (pBA >= pB)*pB + (pBA <= 0)*0 + (pBA > 0)*(pBA < pB)*pBA;
Ord1 = pAB - pBA;
pnBA = pAnB - o2;
pnBA = (pnBA >= pnB)*pnB + (pnBA <= 0)*0 + (pnBA > 0)*(pnBA < pnB)*pnBA;
Ord2 = pAnB - pnBA;
pnBnA = pnB - pnBA;
pBnA = pB - pBA;

% reversed order effects
pnAnB = pnBnA - o3;
pnAnB = (pnAnB >= pnA)*pnA + (pnAnB <= 0)* 0 + (pnAnB > 0)*(pnAnB < pnA)*pnAnB;
Ord3 = pnAnB - pnBnA;
pnAB = pnA - pnAnB;

% conditionals
pAgB = pBA/pB;      pnAgB = 1-pAgB;
pAgnB = pnBA/pnB ; pnAgnB = 1-pAgnB ;

% reversed conditionals
pnBgnA = pnAnB/pnA; pBgnA = 1 - pnBgnA;
%% A & C

pAC =   pA*pCgA ;
pAnC =  pA*pnCgA;

% order effects
pCA = pAC - o4;
pCA = (pCA >= pC)*pC + (pCA <= 0)*0 + (pCA > 0)*(pCA < pC)*pCA;
Ord4 = pAC - pCA;
%Since this is A first, we need interference on C
pnCA = pAnC - o5;
pnCA = (pnCA >= pnC)*pnC + (pnCA <= 0)*0 + (pnCA > 0)*(pnCA < pnC)*pnCA;
Ord5 = pAnC - pnCA;
pnCnA = pnC - pnCA;
pCnA = pC - pCA;

% reversed order effects
pnAnC = pnCnA - o6;
pnAnC = (pnAnC >= pnA)*pnA + (pnAnC <= 0)* 0 + (pnAnC > 0)*(pnAnC < pnA)*pnAnC;
Ord6 = pnAnC - pnCnA;
pnAC = pnA - pnAnC;

% conditionals
pAgC = pCA/pC;      pnAgC = 1-pAgC;
pAgnC = pnCA/pnC ; pnAgnC = 1-pAgnC ;

% reversed conditionals
pnCgnA = pnAnC/pnA; pCgnA = 1 - pnCgnA;
%% B & C

pBC =   pB*pCgB;
pBnC =  pB*pnCgB;

% order effects
pCB = pBC - o7;
pCB = (pCB >= pC)*pC + (pCB <= 0)*0 + (pCB > 0)*(pCB < pC)*pCB;
Ord7 = pBC - pCB;
pnCB = pBnC - o8;
pnCB = (pnCB >= pnC)*pnC + (pnCB <= 0)*0 + (pnCB > 0)*(pnCB < pnC)*pnCB;
Ord8 = pBnC - pnCB;
pnCnB = pnC - pnCB;
pCnB = pC - pCB;

% reversed order effects
pnBnC = pnCnB - o9;
pnBnC = (pnBnC >= pnB)*pnB + (pnBnC <= 0)* 0 + (pnBnC > 0)*(pnBnC < pnB)*pnBnC;
Ord9 = pnBnC - pnCnB;
pnBC = pnB - pnBnC;

% conditionals
pBgC = pCB/pC;      pnBgC = 1-pBgC;
pBgnC = pnCB/pnC ; pnBgnC = 1-pBgnC ;

% reversed conditionals
pnCgnB = pnBnC/pnB; pCgnB = 1 - pnCgnB;
%%
% MoreLikelyFirst goes here
%Order effects -- for conjunctions more likely first, for disjunctions less
%likely first

%AB pairs
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


%AC pairs
if pA > pC
    pAandC = pAC;
    pnAandnC = pnCnA;
else % pC > pA
    pAandC = pCA;
    pnAandnC = pnAnC;
end

if pA > pnC
    pAandnC = pAnC;
    pnAandC = pCnA;
else
    % pnC > pA
    pAandnC = pnCA;
    pnAandC = pnAC;
end

pAorC = 1 - pnAandnC;
pnAorC = 1 - pAandnC;
pAornC = 1 - pnAorC;
pnAornC = 1 - pAandC;



%BC pairs
if pB > pC
    pBandC = pBC;
    pnBandnC = pnCnB;
else % pC > pB
    pBandC = pCB;
    pnBandnC = pnBnC;
end

if pB > pnC
    pBandnC = pBnC;
    pnBandC = pCnB;
else
    % pnC > pB
    pBandnC = pnCB;
    pnBandC = pnBC;
end

pBorC = 1 - pnBandnC;
pnBorC = 1 - pBandnC;
pBornC = 1 - pnBorC;
pnBornC = 1 - pBandC;
%%
% A = A1,  B = A2,  C = A3

%       1       2      3       4         5       6
%      P(A1)   P(A2)  P(A3)  P(¬A1)    P(¬A2) P(¬A3)
Pred = [ pA ;   pB  ;   pC ;    pnA  ;  pnB ;   pnC ; ...

%    7      8            9       10
% P(A2|A1) P(A2|¬A1) P(¬A2|A1) P(¬A2|¬A1)      4

pBgA  ;    pBgnA  ;    pnBgA  ;   pnBgnA ; ...

%   11         12           13          14
% P(A1 ∩ A2) P(A1 ∩ ¬A2) P(¬A1 ∩ A2) P(¬A1 ∩ ¬A2)   4

pAandB     ;     pAandnB  ;      pnAandB ;   pnAandnB ; ...

%   15          16           17          18
% P(A1 ∪ A2) P(A1 ∪ ¬A2)  P(¬A1 ∪ A2) P(¬A1 ∪ ¬A2)  4

pAorB  ;    pAornB   ;    pnAorB   ;    pnAornB ; ...

%    19        20         21        22
%  P(A3|A2) P(A3|¬A2)  P(¬A3|A2) P(¬A3|¬A2)   4

pCgB   ;     pCgnB  ;    pnCgB ;    pnCgnB ; ...

%  23          24              25            26
% P(A2 ∩ A3) P(A2 ∩ ¬A3)    P(¬A2 ∩ A3)    P(¬A2 ∩ ¬A3) 4

pBandC  ;        pBandnC  ;       pnBandC ;      pnBandnC ; ...

%  27               28         29          30
% P(A2 ∪ A3) P(A2 ∪ ¬A3)  P(¬A2 ∪ A3)   P(¬A2 ∪ ¬A3)  4

pBorC  ;        pBornC  ;       pnBorC ;      pnBornC ; ...

%    31         32         33             34
%  P(A1|A3) P(A1|¬A3)   P(¬A1|A3)      P(¬A1|¬A3)   4

pAgC ;     pAgnC  ;    pnAgC   ;    pnAgnC ;  ...

%    35          36            37         38
% P(A3 ∩ A1) P(A3 ∩ ¬A1)    P(¬A3 ∩ A1) P(¬A3 ∩ ¬A1)  4

pAandC    ;     pnAandC    ;    pAandnC ;       pnAandnC  ;  ...

%    39         40          41          42
% P(A3 ∪ A1) P(A3 ∪ ¬A1) P(¬A3 ∪ A1) P(¬A3 ∪ ¬A1)  4

pAorC    ;     pnAorC    ;    pAornC ;       pnAornC  ;  ...

%  43          44          45          46
% P(A1|A2)   P(¬A1|A2)    P(A1|¬A2)  P(¬A1|¬A2)     4

pAgB   ;    pnAgB  ;    pAgnB ;   pnAgnB ;   ...

%   47          48          49             50
% P(A2 ∩ A1) P(¬A2 ∩ A1) P(A2 ∩ ¬A1)    P(¬A2 ∩ ¬A1)    4

pAandB     ;     pAandnB  ;      pnAandB ;   pnAandnB ; ...

%   51           52         53           54
% P(A2 ∪ A1) P(¬A2 ∪ A1) P(A2 ∪ ¬A1)  P(¬A2 ∪ ¬A1)

pAorB  ;    pAornB   ;    pnAorB   ;    pnAornB ;  ...

%    55          56          57         58
% P(A2|A3)   P(¬A2|A3)    P(A2|¬A3)   P(¬A2|¬A3)

pBgC   ;    pnBgC  ;      pBgnC  ;    pnBgnC  ;  ...

%   59          60           61            62
% P(A3 ∩ A2) P(¬A3 ∩ A2)   P(A3 ∩ ¬A2)   P(¬A3 ∩ ¬A2)

pBandC  ;        pBandnC  ;       pnBandC ;      pnBandnC ;   ...

%   63           64               65        66
% P(A3 ∪ A2) P(¬A3 ∪ A2)    P(A3 ∪ ¬A2) P(¬A3 ∪ ¬A2)

pBorC  ;        pBornC  ;       pnBorC ;      pnBornC ; ...

%   67       68            69          70
% P(A3|A1) P(¬A3|A1)     P(A3|¬A1)  P(¬A3|¬A1)

pCgA   ;     pnCgA ;    pCgnA ;    pnCgnA  ; ...

%    71         72           73         74
% P(A1 ∩ A3) P(¬A1 ∩ A3) P(A1 ∩ ¬A3) P(¬A1 ∩ ¬A3)

pAandC    ;     pnAandC    ;    pAandnC ;       pnAandnC  ; ...

%    75          76            77            78
% P(A1 ∪ A3) P(¬A1 ∪ A3)   P(A1 ∪ ¬A3)   P(¬A1 ∪ ¬A3)

pAorC    ;     pnAorC    ;    pAornC ;       pnAornC  ; ];


Ords = [Ord1,Ord2,Ord3,Ord4,Ord5,Ord6,Ord7,Ord8,Ord9];

%% Cross validation
% conj_index = [11 12 13 14 23 24 25 26 35 36 37 38 47 48 49 50 59 60 ...
%     61 62 71 72 73 74];
% disj_index = [15 16 17 18 27 28 29 30 39 40 41 42 51 52 53 54 ... 
%     63 64 65 66 75 76 77 78];
%% Build Measurement Model

LL = eps;
LLc = eps;

Ks = zeros(nd,101);
Ms = zeros(nd,1);

for k = 1:nd
    Rk = Sdat(k,1);   % subj's rating 0 to 100
    Pk = Pred(k,1);   % subjective probability for question
    %Cross validation
%     if ismember(k,disj_index) == 0
        %         K = BuildInt5([alphau,alphad,Pk],m);
    K = drift_add([alphakk,ckk,Pk],m);
    T1 = expm(K) ;
    Psyf = T1*Psy0;
    Ks(k,:) = Psyf;
    Ms(k) = (0:100)*Psyf;
    % We fit the model using cs = 1 for our dataset, but for Zhu et al.
    % (2020) dataset, since there is biased towards 5s and 10s, we put
    % judgments into case of 5s and 10s when fitting.
    cc = cs;
    if (cs == 5) && (Rk >= 95)
        cc = 6;
    end

    PR = eps;

    for j = 1:cc
        Rkj = Rk + j ;   % shifted up one because lowest index =1
        PR = PR + Psyf(Rkj,1);
    end
    LL = LL + log(PR);
    
    %Cross validation test result
%     if ismember(k,disj_index) == 1
%         LLc = LLc + log(PR);
%     end
    % end
end  % states

nLL = -2*LL;
nLLc = -2*LLc;

end   % function


