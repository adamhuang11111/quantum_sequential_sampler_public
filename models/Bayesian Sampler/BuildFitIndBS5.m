% build mex file for quant model
clear
clc


load IndDat.mat
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


parm = [0.5;.5*ones(8,1)]';

codegen FitIndBSBeta5 -args {parm,Cdat,cs}



