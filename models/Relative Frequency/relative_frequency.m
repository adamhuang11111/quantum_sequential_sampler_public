function [SSE, Ms, Pred]= relative_frequency(parm,Cdat,cs)

nd = size(Cdat,1);
%%
% true probabilities
% A = A1,  B = A2,  C = A3

% probabilities: note that pA,pB,pC can not be exactly 0, we can set it to
% be close to 0 like 0.0001, but never exactly 0.
pA = parm(1); pnA = 1 - pA;
pB = parm(2); pnB = 1 - pB;
pC = parm(3); pnC = 1 - pC;
pBgA = parm(4); pnBgA = 1 - pBgA;
pCgA = parm(5); pnCgA = 1 - pCgA;
pCgB = parm(6); pnCgB = 1 - pCgB;
%% A & B
%%%%%%%%%%%%%%%%
% known conjunctions
pAB =   pA*pBgA ;
pAnB =  pA*pnBgA;
pBA = pAB;
%Classical models have these bounds in order to make sure that pAB is not
%greater than pB as a result of pA*pBgA can be fitted to be greater than
%pB.
if (pBA < 0) || (pBA > pB)
    SSE = 10000;
    return
end
pnBA = pAnB;
if (pnBA < 0) || (pnBA > pnB)
    SSE = 10000;
    return
end
pnBnA = pnB - pnBA;
pBnA = pB - pBA;
pnAnB = pnBnA;
if (pnAnB < 0) || (pnAnB > pnA)
    SSE = 10000;
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
    SSE = 10000;
    return
end
pnCA = pAnC;
if (pnCA < 0) || (pnCA > pnC)
    SSE = 10000;
    return
end
pnCnA = pnC - pnCA;
pCnA = pC - pCA;
pnAnC = pnCnA;
if (pnAnC < 0) || (pnAnC > pnA)
    SSE = 10000;
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
    SSE = 10000;
    return
end
pnCB = pBnC;
if (pnCB < 0) || (pnCB > pnC)
    SSE = 10000;
    return
end
pnCnB = pnC - pnCB;
pCnB = pC - pCB;
pnBnC = pnCnB;
if (pnBnC < 0) || (pnBnC > pnB)
    SSE = 10000;
    return
end
pnBC = pnB - pnBnC; 
pBgC = pCB/pC;      pnBgC = 1-pBgC;
pBgnC = pnCB/pnC ; pnBgnC = 1-pBgnC ;
pnCgnB = pnBnC/pnB; pCgnB = 1 - pnCgnB;
%%

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

SSE = eps;
Ms = zeros(nd,1);

for k = 1:nd
    Rk = Cdat(k,1);
    Pk = Pred(k,1);
    SE = eps;
    if cs == 5
        for j = 0:(cs-1)
            Rkj = (Rk + j)/100;
            SE = SE + (Pk - Rkj).^2;
        end
    else 
        Rkj = Rk/100;
        % The relative frequency model is fitted by square error
        SE = SE + (Pk - Rkj).^2;
    end
    SSE = SSE + SE;
    Ms(k,1) = Pk*100;
end
