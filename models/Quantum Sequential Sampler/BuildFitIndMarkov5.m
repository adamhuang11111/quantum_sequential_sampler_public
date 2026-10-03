% build mex file for quant model
clear
clc


here = fileparts(mfilename('fullpath'));   % Quantum Sequential Sampler folder
root = fileparts(fileparts(here));
cd(here)   % codegen puts the mex file in the current folder
load(fullfile(root, 'data', 'IndDat.mat'))
% Contains Triplet Comp Rdat code

subj = 1;
Sdat = double(Rdat(1,:))';


% if you change cs, you need to rerun this program for the mex file
cs = 1;

if cs == 5
    Cdat = floor(Sdat/cs) * cs;
    Cdat = (Cdat == 100).*95 + (Cdat < 100).*Cdat;
else
    Cdat = Sdat;
end


parm = [0.5;.5*ones(8,1);0]';   % 10 parms: parm(10) is the interference term

codegen FitIndMarkov5_qp_int1_qq -args {parm,Cdat,cs}



