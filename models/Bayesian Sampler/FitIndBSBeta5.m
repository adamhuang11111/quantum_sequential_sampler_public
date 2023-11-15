function [nLL, nLLc, Pred, Ms, Ks] = FitIndBSBeta5(parm,Cdat,cs)

beta = parm(7);
N2 = 1 + round(parm(9));
% N3 = 1 + round(parm(10));
% N4 = 1 + round(parm(11));
N3 = N2;
N4 = N2;
N1 = max([N2 N3 N4]) + round(parm(8));

% x_all is used to compute the beta likelihood (use 0.005 to approximate
% 0, and 0.995 to approximate 1)
x_all = [ 0.005 (.01:.01:.99) 0.995];
nd = size(Cdat,1);
nLLc = 0;
Pred = zeros(78,1);
Ks = zeros(nd,101);
Ms = zeros(nd,1);


%%
% true probabilities
% A = A1,  B = A2,  C = A3
pA = parm(1); pnA = 1 - pA;
pB = parm(2); pnB = 1 - pB;
pC = parm(3); pnC = 1 - pC;
pBgA = parm(4); pnBgA = 1 - pBgA;
pCgA = parm(5); pnCgA = 1 - pCgA;
pCgB = parm(6); pnCgB = 1 - pCgB;

%% A & B
%%%%%%%%%%%%%%%%
pAB =   pA*pBgA ;
pAnB =  pA*pnBgA;
pBA = pAB;
%Classical models have these bounds in order to make sure that pAB is not
%greater than pB as a result of pA*pBgA can be fitted to be greater than
%pB.
if (pBA < 0) || (pBA > pB)
    nLL = 10000;
    return
end
pnBA = pAnB;
if (pnBA < 0) || (pnBA > pnB)
    nLL = 10000;
    return
end
pnBnA = pnB - pnBA;
pBnA = pB - pBA;

pnAnB = pnBnA;
if (pnAnB < 0) || (pnAnB > pnA)
    nLL = 10000;
    return
end
pnAB = pnA - pnAnB; 
pAgB = pBA/pB;      pnAgB = 1-pAgB;
pAgnB = pnBA/pnB ; pnAgnB = 1-pAgnB ;
pnBgnA = pnAnB/pnA; pBgnA = 1 - pnBgnA;
%% A & C
pAC =   pA*pCgA ;
pAnC =  pA*pnCgA;
pCA = pAC;
if (pCA < 0) || (pCA > pC)
    nLL = 10000;
    return
end
pnCA = pAnC;
if (pnCA < 0) || (pnCA > pnC)
    nLL = 10000;
    return
end
pnCnA = pnC - pnCA;
pCnA = pC - pCA;
pnAnC = pnCnA;
if (pnAnC < 0) || (pnAnC > pnA)
    nLL = 10000;
    return
end
pnAC = pnA - pnAnC; 
pAgC = pCA/pC;      pnAgC = 1-pAgC;
pAgnC = pnCA/pnC ; pnAgnC = 1-pAgnC ;
pnCgnA = pnAnC/pnA; pCgnA = 1 - pnCgnA;
%% B & C
pBC =   pB*pCgB;
pBnC =  pB*pnCgB;
pCB = pBC;
if (pCB < 0) || (pCB > pC)
    nLL = 10000;
    return
end
pnCB = pBnC;
if (pnCB < 0) || (pnCB > pnC)
    nLL = 10000;
    return
end
pnCnB = pnC - pnCB;
pCnB = pC - pCB;
pnBnC = pnCnB;
if (pnBnC < 0) || (pnBnC > pnB)
    nLL = 10000;
    return
end
pnBC = pnB - pnBnC; 
pBgC = pCB/pC;      pnBgC = 1-pBgC;
pBgnC = pnCB/pnC ; pnBgnC = 1-pBgnC ;
pnCgnB = pnBnC/pnB; pCgnB = 1 - pnCgnB;
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

pAB     ;     pAnB  ;      pnAB ;   pnAnB ; ...

%   15          16           17          18
% P(A1 ∪ A2) P(A1 ∪ ¬A2)  P(¬A1 ∪ A2) P(¬A1 ∪ ¬A2)  4

1-pnAnB  ;    1-pnAB   ;    1-pAnB   ;    1-pAB ; ...

%    19        20         21        22
%  P(A3|A2) P(A3|¬A2)  P(¬A3|A2) P(¬A3|¬A2)   4

pCgB   ;     pCgnB  ;    pnCgB ;    pnCgnB ; ...

%  23          24              25            26
% P(A2 ∩ A3) P(A2 ∩ ¬A3)    P(¬A2 ∩ A3)    P(¬A2 ∩ ¬A3) 4

pBC  ;        pBnC  ;       pnBC ;      pnBnC ; ...

%  27               28         29          30
% P(A2 ∪ A3) P(A2 ∪ ¬A3)  P(¬A2 ∪ A3)   P(¬A2 ∪ ¬A3)  4

1-pnBnC   ;  1-pnBC   ;    1-pBnC  ;    1-pBC ; ...

%    31         32         33             34
%  P(A1|A3) P(A1|¬A3)   P(¬A1|A3)      P(¬A1|¬A3)   4

pAgC ;     pAgnC  ;    pnAgC   ;    pnAgnC ;  ...

%    35          36            37         38
% P(A3 ∩ A1) P(A3 ∩ ¬A1)    P(¬A3 ∩ A1) P(¬A3 ∩ ¬A1)  4

pCA  ;       pCnA   ;      pnCA  ;     pnCnA ;  ...

%    39         40          41          42
% P(A3 ∪ A1) P(A3 ∪ ¬A1) P(¬A3 ∪ A1) P(¬A3 ∪ ¬A1)  4

1-pnCnA   ;   1-pnCA  ;    1-pCnA  ;    1-pCA ;  ...

%  43          44          45          46
% P(A1|A2)   P(¬A1|A2)    P(A1|¬A2)  P(¬A1|¬A2)     4

pAgB   ;    pnAgB  ;    pAgnB ;   pnAgnB ;   ...

%   47          48          49             50
% P(A2 ∩ A1) P(¬A2 ∩ A1) P(A2 ∩ ¬A1)    P(¬A2 ∩ ¬A1)    4

pBA     ;    pnBA   ;     pBnA   ;     pnBnA ; ...

%   51           52         53           54
% P(A2 ∪ A1) P(¬A2 ∪ A1) P(A2 ∪ ¬A1)  P(¬A2 ∪ ¬A1)

1-pnBnA  ;    1-pBnA  ;   1-pnBA   ;    1-pBA ;  ...

%    55          56          57         58
% P(A2|A3)   P(¬A2|A3)    P(A2|¬A3)   P(¬A2|¬A3)

pBgC   ;    pnBgC  ;      pBgnC  ;    pnBgnC  ;  ...

%   59          60           61            62
% P(A3 ∩ A2) P(¬A3 ∩ A2)   P(A3 ∩ ¬A2)   P(¬A3 ∩ ¬A2)

pCB   ;       pnCB   ;      pCnB ;      pnCnB ;   ...

%   63           64               65        66
% P(A3 ∪ A2) P(¬A3 ∪ A2)    P(A3 ∪ ¬A2) P(¬A3 ∪ ¬A2)

1-pnCnB  ;   1-pCnB   ;    1-pnCB ;    1-pCB  ; ...

%   67       68            69          70
% P(A3|A1) P(¬A3|A1)     P(A3|¬A1)  P(¬A3|¬A1)

pCgA   ;     pnCgA ;    pCgnA ;    pnCgnA  ; ...

%    71         72           73         74
% P(A1 ∩ A3) P(¬A1 ∩ A3) P(A1 ∩ ¬A3) P(¬A1 ∩ ¬A3)

pAC    ;     pnAC    ;    pAnC ;       pnAnC  ; ...

%    75          76            77            78
% P(A1 ∪ A3) P(¬A1 ∪ A3)   P(A1 ∪ ¬A3)   P(¬A1 ∪ ¬A3)

1-pnAnC  ;    1-pAnC  ;   1-pnAC  ;     1-pAC ];

%%  Partition different part of the predictions into different pairs (not needed for the base model)
% ListAB = [ 11:18 47:54 ]' ;
% ListBC = [ 23:30 59:66 ]' ;
% ListAC = [ 35:42 71:78]' ;
ListCFDF = [ 11:18 23:30 35:42 47:54 59:66 71:78]' ;
% sL = size(ListAC,1);
sL = size(ListCFDF,1);
%% Cross validation
% conj_index = [11 12 13 14 23 24 25 26 35 36 37 38 47 48 49 50 59 60 ...
%     61 62 71 72 73 74];
% disj_index = [15 16 17 18 27 28 29 30 39 40 41 42 51 52 53 54 ... 
%     63 64 65 66 75 76 77 78];
%% Fitting
LL = eps;
LLc = eps;
norms_N1 = norm_beta(N1,beta,x_all);
norms_N2 = norm_beta(N2,beta,x_all);
% norms_N3 = norm_sampler(N3,beta,x_all);
% norms_N4 = norm_sampler(N4,beta,x_all);

for k = 1:nd     % loop thru 78 questions
    %cross validation
%     if ismember(k,disj_index) == 0
        Rk = Cdat(k,1);   % subj's rating 0 to 100
        Pk = Pred(k,1);   % subjective probability for question
        EL = sum(k*ones(sL,1) == ListCFDF);
        N = EL.*N2 + (1-EL).*N1;
    %     EAB = sum(k*ones(sL,1) == ListAB);
    %     EBC = sum(k*ones(sL,1) == ListBC);
    %     EAC = sum(k*ones(sL,1) == ListAC);
    %     N = EAB.*N2 + EBC.*N3 + EAC.*N4 + (1-EAB-EBC-EAC).*N1;
        if N == N1
            norms = norms_N1;
        else
            norms = norms_N2;
        end
    %     elseif N == N3
    %         norms = norms_N3;
    %     else
    %         norms = norms_N4;
    %     end
        PR = eps;
        % We fit the model using cs = 1 for our dataset, but for Zhu et al.
        % (2020) dataset, since there is biased towards 5s and 10s, we put
        % judgments into case of 5s and 10s when fitting.
        if cs == 5
            for j = 0:(cs-1)
                Rkj = (Rk + j)/100 +.005;
                PR = PR +  beta_sampler(N, Rkj, beta, Pk, norms);
            end
        else 
            Rkj = (Rk==0).*(.005) + (Rk==100).*(.995) + (Rk>0).*(Rk<100).*(Rk/100);
            PR = PR + beta_sampler(N, Rkj, beta, Pk, norms);
%         end
        ms = 0;
        % Compute mean and distribution for each participant.
%         for ij = 1:101
%             ms = ms + beta_sampler(N, x_all(ij), beta, Pk, norms)* (ij-1);
%             Ks(k,ij) = beta_sampler(N, x_all(ij), beta, Pk, norms);
%         end
%         Ms(k) = ms;
        LL = LL + log(PR);
%       compute cross validation test
%         if ismember(k,disj_index) == 1
%             LLc = LLc + log(PR);
%         end
    %end
end  % states

nLL = -2*LL;
nLLc = -2*LLc;

end   % function


