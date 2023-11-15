% Parallel compute across subjects
clear
clc
close('all')

file_list = ["exp1_frosty", "exp1_normal", "exp2_cloudy", "exp2_rainy", "exp2_snowy"];

for file = file_list
    
    filename = strcat(file, ".mat");
    
    load(filename);
    % Contains Triplet Comp Rdat code
    Ns = size(data,2);   % no conditions
    Nd = size(data,1);
    options = optimoptions('particleswarm','SwarmSize',160,'UseParallel',true,'Display','off','MaxIter',300);

    parm0 = [ .5*ones(3,1); 1.1; 50; 50] ;
    np = size(parm0,1);
    
    lb = .001*zeros(np,1);
    ub = .999*ones(np,1);
    ub(5:6,1) = 200; % Bounds for the two sample sizes
    ub(4,1) = 200; %Bound for beta parameter

    nLLS =  zeros(Ns/3,1);
    ParmS = zeros(Ns/3,np);
    MeanS = zeros(Ns/3,Nd);

    ns = 2;    % number of subjects to run , set equal to Ns to run all
    reps = 5;  % number of replicatins per subject
    cs = 5;   % categorization, if changed , re-run BuildFitIndQuant3
    cs_bol = 1; % whether the rounding mechanism is needed (1 means turn on
    ... rounding and 0 means turn off)
    count = 1;

    for subj = linspace(1,(Ns - 2),Ns/3)
%     for subj = 1:N
        Sdat = double(data(:,subj:subj+2));
        %Round data to nearest 5 and tens for the rounding mechanism
        if cs_bol == 1
            Cdat = floor(Sdat/cs) * cs;
            %100 round to 95 and else just keep it.
            Cdat = (Cdat == 100).*(100-cs) + (Cdat < 100).*Cdat;
        else
            Cdat = Sdat;
            cs = 1;
        end

        [nLL,parm] = MainFitIndBS3(Cdat,parm0,lb,ub,options,reps,cs);
        nLLS(count) = nLL;
        ParmS(count,:) = parm';
        count = count + 1;
        disp(strcat(strcat(file," "), string(count - 1)))
        disp(nLL);
    end
    newnLL = strcat(file,"_bs_nLLS.mat");
    newParm = strcat(file,"_bs_ParmS.mat");
    save(newnLL,"nLLS");
    save(newParm,"ParmS");
end