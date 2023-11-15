% Main Loop  across subjects

clear
clc
close('all')
% 
file_list = ["exp1_frosty","exp1_normal", "exp2_cloudy", "exp2_rainy", "exp2_snowy"];

for file = file_list

% Contains Triplet Comp Rdat code
    filename = strcat(file, ".mat");
    load(filename);
    
    Ns = size(data,2);    % no subj
    Nd = size(data,1);    % no conditions
    
    np = 6;
%     np = 7; % no parameters
    nLLS =  zeros(Ns/3,1);
    ParmS = zeros(Ns/3,np);
    MeanS = zeros(Ns/3,Nd);

    options = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'Display','off','MaxIter',1000);


    ns = 2;    % number of subjects to run , set equal to Ns to run all
    reps = 3;  % number of replicatins per subject
    cs = 5;   % categorization, if changed , re-run BuildFitIndQuant3
    count = 1;
    
    for subj = linspace(1,(Ns - 2),Ns/3)
        lb = .001*zeros(np,1);
        ub = .999*ones(np,1);
        ub(4) = 200;  %bound of drift rate
        ub(5) = 200;  %bound of additive bias
        lb(5) = -200; %bound of additive bias
        ub(6) = 200; %bound of beta parameter
%         lb(7) = -0.999; %bound of interference parameter
        Sdat = double(data(:,subj:subj+2));
        %Round data to nearest 5 and tens for the rounding mechanism
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
            BSM = @(parm) FitIndMarkov3_classical_mex(parm,Cdat,cs);
            [parm, nLL] = particleswarm(BSM,np,lb,ub,options);
            nLLV(n) =  nLL;
            ParmM(n,:) =  parm';       
        end  % reps

        [nLL, Ind] = min(nLLV);    % pick best fit Index
        parm = ParmM(Ind,:);
        % parm comes out of optimization as a row vector
        nLLS(count) = nLL;
        ParmS(count,:) = parm';
        count = count + 1;
        disp(strcat(strcat(file," "), string(count - 1)))
        disp(nLL)
    end
    
    newnLL = strcat(file,"_classical_nLLS.mat");
    newParm = strcat(file,"_classical_ParmS.mat");
    save(newnLL,"nLLS");
    save(newParm,"ParmS");

end