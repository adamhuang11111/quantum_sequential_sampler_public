% Pred_all = zeros(size(Rdat));
Ns = size(data,2);
nLLS_all = zeros(Ns/3,1);
Ms_all = zeros(Ns/3,20);
%Average data
Ms_data_all = zeros(Ns/3,20);
subjn = 1;
for subj = linspace(1,(Ns - 2),Ns/3)
    Sdat = double(data(:,subj:subj+2));
    Ms_data_all(subjn,:) = mean(Sdat,2);
    parm = ParmS(subjn,:);
    %Round data to nearest 5 and tens for the rounding mechanism
    cs = 5;
    if cs == 5
    Cdat = floor(Sdat/cs) * cs;
    Cdat = (Cdat == 100).*(100-cs) + (Cdat < 100).*Cdat;
    else
        Cdat = Sdat;
    end
    [nLL, Pred, Ms] = FitIndBSBeta3(parm,Cdat,cs);
%     Pred_all(subjn,:)= Pred;
    nLLS_all(subjn) = nLL;
    Ms_all(subjn,:) = Ms;
    disp(subjn)
    subjn = subjn + 1;
end