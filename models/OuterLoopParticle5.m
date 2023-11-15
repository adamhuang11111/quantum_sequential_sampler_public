% Main Loop  across subjects

clear
clc
close('all')

 

% Contains Triplet Comp Rdat code
load IndDat

Ns = size(Rdat,1);    % no subj
Nd = size(Rdat,2);    % no conditions

np = 9; % no parameters



nLLS =  zeros(Ns,1);
ParmS = zeros(Ns,np);
MeanS = zeros(Ns,Nd);


options = optimoptions('particleswarm','SwarmSize',50,'UseParallel',true,'Display','off','MaxIter',1000);


ns = 2;    % number of subjects to run , set equal to Ns to run all
reps = 3;  % number of replicatins per subject
cs = 1;   % categorization, if changed , re-run BuildFitIndQuant3
count = 1;

%     for subj = [1,4]
for subj =  601:1162
    %% Quantum upper and lower bounds
    lb = .0001*zeros(np,1);
    ub = .9999*ones(np,1);
    ub(7) = 200; %upper bound drift 
    ub(8) = 200; %upper bound additive bias
    ub(9) = 200; %upper bound symmetric beta
    lb(8) = -200; %lower bound additive bias
%     lb(10) = -0.9999; %interference
    %% Bayesian Sampler upper and lower bounds
%     lb = .0001*zeros(np,1);
%     ub = .9999*ones(np,1);
%     ub(7) = 200; %upper bound symmetric beta
%     ub(8) = 300; %upper bound sample size 1
%     ub(9) = 300; %upper bound sample size 2 (marginals = sample size 1 + 2)
    %% Fitting
    Sdat = double(Rdat(subj, :))';

    if cs == 5
        Cdat = floor(Sdat/cs) * cs;
        Cdat = (Cdat == 100).*(100-cs) + (Cdat < 100).*Cdat;
    else
        Cdat = Sdat;
    end

    nLLV = zeros(reps,1);
    ParmM = zeros(reps,np);
    wgt = .2;    % amount of jitter
    for n = 1:reps
        %change it to the model you want to fit
        BSM = @(parm) FitIndMarkov5_qp_int1_qq_mex(parm,Cdat,cs);
        [parm, nLL] = particleswarm(BSM,np,lb,ub,options);
        nLLV(n) =  nLL;
        ParmM(n,:) =  parm';    
    end  % reps
    [nLL, Ind] = min(nLLV);    % pick best fit Index
    parm = ParmM(Ind,:);
    nLLS(subj) = nLL;
    ParmS(subj,:) = parm';
    disp(string(subj - 1))
    disp(nLL)
end
%% Save data
% newnLL = "IndDat_int1_qq_nLLS.mat";
% newParm = "IndDat_int1_qq_ParmS.mat";
% save(newnLL,"nLLS");
% save(newParm,"ParmS");