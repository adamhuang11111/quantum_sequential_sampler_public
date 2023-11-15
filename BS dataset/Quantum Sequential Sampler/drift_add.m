function K = drift_add(P,m)

 alpha = P(1); % drift rate 
 c = P(2);
 prob = P(3); % subjective probability

if c >= 0
    c1 = c;
    c2 = 0;
else
    c1 = 0;
    c2 = -c;
end
    
 
up = prob*alpha*ones(m-2,1) + 1 + c1;    %rate up = (prob dn)/h
dn = (1 - prob)*alpha*ones(m-2,1) + 1 + c2;    % rate dn = (prob dn)/h
cn = -(up + dn);

A = diag(up);
B = diag(dn);
C = diag(cn);
K = zeros(m,m);


K(2:(m-1),2:(m-1)) = K(2:(m-1),2:(m-1)) + C;
K(1:(m-2),2:(m-1)) = K(1:(m-2),2:(m-1)) + B;     
K(3:m,2:(m-1)) = K(3:m,2:(m-1)) + A;
K(2,1) = up(1,1);
K(1,1) = -K(2,1);
K(100,101) = dn(1,1);
K(101,101) = -K(100,101);

% in col, out row
% col's sum to zero
