nLLS_all = zeros(1162,1);
% nLLSc_all = zeros(1162,1);
% Pred_all = zeros(size(Rdat));
% Psyfs_all = zeros(size(Rdat,1),101,size(Rdat,2));
Ms_all = zeros(size(Rdat));
% Ords_all = zeros(size(Rdat,1),9);
for i = 1:1162
    Sdat = double(Rdat(i,:))';
    parm = ParmS(i,:);
    cs = 1;
    [nLL, nLLc, Pred, Psy0,Ks,Ms,Ords] = FitIndMarkov5_qp_classical(parm,Sdat,cs);
%     [nLL, Ms] = relative_frequency(parm,Sdat,cs);
%     [nLL, nLLc, Pred, Ms, Ks] = FitIndBSBeta5(parm,Sdat,cs);
    nLLS_all(i) = nLL;
%     nLLSc_all(i) = nLLc;
%     Pred_all(i,:)= Pred;
%     Psyfs_all(i,:,:) = Psyfs;
    Ms_all(i,:) = Ms;
%     Ords_all(i,:) = Ords;
    disp(i)
end