% Generalization test across subjects
% Fits the model to 54 of the 78 questions and reports the G square on the
% 24 held out questions.
%   heldout = 1  train without the disjunctions, test on the disjunctions
%                -> IndDat_<model>_test_disj_*
%   heldout = 2  train without the conjunctions, test on the conjunctions
%                -> IndDat_<model>_test_conj_*
% Results are saved in model_generalization_test. Progress is saved after
% every subject; rerunning resumes where it stopped.

clear
clc
close('all')

model = "classical";   % "int1_qq" (quantum) or "classical"
heldouts = [1 2];      % splits to run in turn: 1 = hold out disjunctions, 2 = hold out conjunctions

here = fileparts(mfilename('fullpath'));   % models folder
root = fileparts(here);

switch model
    case "int1_qq"
        np = 10; % no parameters (int1 quantum: 9 + interference)
        fitfun = @FitIndMarkov5_qp_int1_qq;
        addpath(fullfile(here, 'Quantum Sequential Sampler'))
    case "classical"
        np = 9;  % no parameters (classical: no interference term)
        fitfun = @FitIndMarkov5_qp_classical;
        addpath(fullfile(here, 'Classical Sequential Sampler'))
end

% Contains Triplet Comp Rdat code
load(fullfile(root, 'data', 'IndDat.mat'))

Ns = size(Rdat,1);    % no subj

for heldout = heldouts
    if heldout == 1
        tag = "test_disj";
    else
        tag = "test_conj";
    end
    stem = fullfile(root, 'model_generalization_test', "IndDat_" + model + "_" + tag);
    progress = stem + "_progress.mat";

    if isfile(progress)
        load(progress, "nLLStrain", "nLLStest", "ParmS", "done")
    else
        nLLStrain = zeros(Ns,1);
        nLLStest = zeros(Ns,1);
        ParmS = zeros(Ns,np);
        done = false(Ns,1);
    end

    options = optimoptions('particleswarm','SwarmSize',50,'UseParallel',true,'Display','off','MaxIter',1000);

    reps = 3;  % number of replicatins per subject
    cs = 1;    % categorization

    for subj = 1:Ns
        if done(subj)
            continue
        end
        %% Quantum upper and lower bounds
        lb = .0001*zeros(np,1);
        ub = .9999*ones(np,1);
        ub(7) = 200; %upper bound drift
        ub(8) = 200; %upper bound additive bias
        ub(9) = 200; %upper bound symmetric beta
        lb(8) = -200; %lower bound additive bias
        if np == 10
            lb(10) = -0.9999; %interference
        end
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
        for n = 1:reps
            BSM = @(parm) fitfun(parm,Cdat,cs,heldout);   % training G square only
            [parm, nLL] = particleswarm(BSM,np,lb,ub,options);
            nLLV(n) =  nLL;
            ParmM(n,:) =  parm;
        end  % reps

        [nLLtrain, Ind] = min(nLLV);    % pick best fit Index
        parm = ParmM(Ind,:);
        [~, nLLtest] = fitfun(parm,Cdat,cs,heldout);   % G square on the held out questions
        nLLStrain(subj) = nLLtrain;
        nLLStest(subj) = nLLtest;
        ParmS(subj,:) = parm;
        done(subj) = true;
        save(progress, "nLLStrain", "nLLStest", "ParmS", "done", "model", "heldout", "cs", "reps")
        disp(string(subj - 1))
        disp([nLLtrain nLLtest])
    end

    %% Save data
    if all(done)
        nLLS = nLLStrain;    % same variable names as the published files
        nLLSc = nLLStest;
        save(stem + "_nLLStrain.mat", "nLLS")
        save(stem + "_nLLStest.mat", "nLLSc")
        save(stem + "_ParmS.mat", "ParmS")
    end
end  % heldouts
