% build mex file for quant model
clear
clc


load exp1_frosty.mat
% Contains Triplet Comp Rdat code

cs = 5;
Ns = size(data,2);

subj = 1;

Sdat = double(data(:,subj:subj+2));

if cs == 5
    Cdat = floor(Sdat/cs) * cs;
    Cdat = (Cdat == 100).*(100-cs) + (Cdat < 100).*Cdat;
else
    Cdat = Sdat;
end


parm = [ .5*ones(3,1) ; 1.1; 50; 5]' ;

codegen FitIndBSBeta3 -args {parm,Cdat,cs}

