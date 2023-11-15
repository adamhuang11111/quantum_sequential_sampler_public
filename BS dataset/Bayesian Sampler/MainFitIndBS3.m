function [nLL,parm] = MainFitIndBS3(Cdat,parm0,lb,ub,options,reps,cs)

np = size(parm0,1);
nLLV = zeros(reps,1);
ParmM = zeros(reps,np);

for n = 1:reps
    BSM = @(parm) FitBSBeta3_mex(parm,Cdat,cs);
    [parm,nLL] = particleswarm(BSM,np,lb,ub,options);
    nLLV(n) =  nLL;
    ParmM(n,:) =  parm';
end  % reps

[nLL, Ind] = min(nLLV);    % pick best fit Index
parm = ParmM(Ind,:);


